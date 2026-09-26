#include "ramcore/SamToNTuple.h"
#include "ramcore/QualityBlocks.h"
#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TStopwatch.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>

#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <TROOT.h>

#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

void HandleHeaderLine(TList &headers, const std::string &tag, const std::string &content)
{
   headers.Add(new TNamed(tag.c_str(), content.c_str()));

   if (tag == "@SQ") {
      size_t sn_pos = content.find("SN:");
      if (sn_pos != std::string::npos) {
         sn_pos += 3;
         size_t tab_pos = content.find('\t', sn_pos);
         std::string ref_name =
            content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
         RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name);
      }
   }
}

// Everything but the two reference ids, which the caller resolves.
void FillRecordFields(const ramcore::SamRecord &sam_record, RAMNTupleRecord &rec, uint32_t quality_policy)
{
   rec.SetBit(quality_policy);
   rec.SetQNAME(sam_record.qname);
   rec.SetFLAG(sam_record.flag);
   rec.SetPOS(sam_record.pos);
   rec.SetMAPQ(sam_record.mapq);
   rec.SetCIGAR(sam_record.cigar);
   rec.SetPNEXT(sam_record.pnext);
   rec.SetTLEN(sam_record.tlen);
   rec.SetSEQ(sam_record.seq);
   rec.SetQUAL(sam_record.qual);

   rec.ResetNOPT();
   for (const auto &opt : sam_record.optional_fields)
      rec.SetOPT(opt);
}

} // namespace

namespace {

// One output file, held open while the input streams past it.
struct ChromosomeWriter {
   std::unique_ptr<TFile> file{};
   std::unique_ptr<ROOT::RNTupleWriter> writer{};
   std::unique_ptr<QualityBlockWriter> out;
   int64_t rows = 0;
   int32_t last_pos = -1;
   bool sorted = true;
};

} // namespace

void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy)
{
   RAMNTupleRecord::InitializeRefs();

   std::map<std::string, ChromosomeWriter> writers;
   TList headers;
   headers.SetName("headers");
   headers.SetOwner(true);

   auto header_callback = [&](const std::string &tag, const std::string &content) {
      HandleHeaderLine(headers, tag, content);
   };

   auto open_writer = [&](const std::string &chr) -> ChromosomeWriter & {
      auto it = writers.find(chr);
      if (it != writers.end())
         return it->second;

      ChromosomeWriter &cw = writers[chr];
      std::string filename{output_prefix};
      filename += "_";
      filename += chr;
      filename += ".root";
      cw.file.reset(TFile::Open(filename.c_str(), "RECREATE"));
      if (!cw.file || !cw.file->IsOpen())
         throw std::runtime_error("cannot create " + filename);

      ROOT::RNTupleWriteOptions writeOptions;
      writeOptions.SetCompression(compression_algorithm);
      writeOptions.SetMaxUnzippedPageSize(64000);
      // Every chromosome's buffers are open at once, so keep the clusters small.
      writeOptions.SetApproxZippedClusterSize(8 * 1024 * 1024);

      cw.writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *cw.file, writeOptions);
      auto *writer = cw.writer.get();
      cw.out = std::make_unique<QualityBlockWriter>(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                                                    [writer](ROOT::REntry &e) { writer->Fill(e); });
      return cw;
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
      // A record with no reference has no chromosome file to go to.
      if (sam_record.rname == "*")
         return;

      ChromosomeWriter &cw = open_writer(sam_record.rname);
      RAMNTupleRecord &rec = cw.out->Record();

      FillRecordFields(sam_record, rec, quality_policy);
      rec.SetREFID(sam_record.rname);
      rec.SetREFNEXT(sam_record.rnext);

      RAMNTupleRecord::NoteRefSpan(rec.GetRefSpan());
      cw.out->Add();
      cw.rows++;

      // One reference per file, so the order check is on the position alone.
      const int32_t pos = rec.GetPOS() - 1;
      if (pos < cw.last_pos)
         cw.sorted = false;
      cw.last_pos = pos;
   };

   ramcore::SamParser parser;
   if (!parser.ParseFile(datafile, header_callback, record_callback)) {
      printf("Failed to parse SAM file %s\n", datafile);
      return;
   }

   // The reference table and the longest span are only complete once the whole
   // input has been read, so every file is finished here.
   for (auto &[chr, cw] : writers) {
      cw.out->Finish();
      RAMNTupleRecord::SetQualBlockEnds(cw.out->TakeBlockEnds());
      cw.out.reset();
      cw.writer.reset();

      RAMNTupleRecord::SetCoordinateSorted(cw.sorted);
      if (!cw.sorted)
         fprintf(stderr, "%s: %s is not in coordinate order; region queries will read it in full.\n", datafile,
                 chr.c_str());
      RAMNTupleRecord::WriteAllRefs(*cw.file);

      cw.file->cd();
      headers.Write("headers", TObject::kSingleKey);
      cw.file->Close();

      printf("%s_%s.root: %lld records\n", output_prefix, chr.c_str(), static_cast<long long>(cw.rows));
   }
}

// Single-file conversion: the main thread reads blocks of whole lines, worker
// threads fill them into their own RNTupleFillContext and commit the staged
// clusters in block order.

namespace {

struct BlockOrder {
   RAMCoordinateOrder check;
   bool has_placed = false;
   int32_t first_refid = -1;
   int32_t first_pos = -1;
   uint32_t max_span = 0;

   void Note(int32_t refid, int32_t pos)
   {
      if (refid >= 0 && !has_placed) {
         has_placed = true;
         first_refid = refid;
         first_pos = pos;
      }
      check.Note(refid, pos);
   }
};

struct FileOrder {
   RAMCoordinateOrder check;
   uint32_t max_span = 0;

   void Add(const BlockOrder &b)
   {
      if (!b.check.sorted)
         check.sorted = false;
      if (b.has_placed) {
         check.Note(b.first_refid, b.first_pos);
         check.Note(b.check.last_refid, b.check.last_pos);
      }
      if (b.check.seen_unplaced)
         check.seen_unplaced = true;
      max_span = std::max(max_span, b.max_span);
   }
};

struct Block {
   size_t seq = 0;        ///< Position in the input, from 0.
   size_t first_line = 0; ///< 1-based number of the block's first line, for warnings.
   std::vector<char> data;
};

// Reads blocks of whole lines; every block ends with '\n'.
class BlockReader {
   FILE *m_file;
   size_t m_blockBytes;
   std::vector<char> m_carry;
   bool m_eof = false;
   bool m_error = false;

public:
   BlockReader(FILE *file, size_t block_bytes) : m_file(file), m_blockBytes(std::max<size_t>(block_bytes, 1)) {}

   /// A read error ends the input early; Next() then returns false.
   [[nodiscard]] bool Failed() const { return m_error; }

   bool Next(std::vector<char> &out)
   {
      out.swap(m_carry);
      m_carry.clear();
      while (!m_eof) {
         const size_t old = out.size();
         out.resize(old + m_blockBytes);
         const size_t n = fread(&out[old], 1, m_blockBytes, m_file);
         out.resize(old + n);
         if (n < m_blockBytes) {
            if (ferror(m_file)) {
               m_error = true;
               return false;
            }
            // A short read from a pipe is not the end.
            if (feof(m_file))
               m_eof = true;
         }
         const auto nl = std::find(out.rbegin(), out.rend(), '\n');
         if (nl != out.rend()) {
            m_carry.assign(nl.base(), out.end());
            out.erase(nl.base(), out.end());
            return true;
         }
      }
      if (out.empty())
         return false;
      if (out.back() != '\n')
         out.push_back('\n');
      return true;
   }
};

// Bounded FIFO. Not TBB (TTaskGroup, TThreadExecutor): it does not start tasks
// in order, and a task waiting for an earlier block could then deadlock.
class BlockQueue {
   std::mutex m_mutex;
   std::condition_variable m_notEmpty;
   std::condition_variable m_notFull;
   std::deque<Block> m_blocks;
   size_t m_capacity;
   bool m_closed = false;

public:
   explicit BlockQueue(size_t capacity) : m_capacity(std::max<size_t>(capacity, 1)) {}

   void Push(Block &&block)
   {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_notFull.wait(lock, [&] { return m_blocks.size() < m_capacity || m_closed; });
      if (m_closed)
         return;
      m_blocks.push_back(std::move(block));
      m_notEmpty.notify_one();
   }

   bool Pop(Block &block)
   {
      std::unique_lock<std::mutex> lock(m_mutex);
      m_notEmpty.wait(lock, [&] { return !m_blocks.empty() || m_closed; });
      if (m_blocks.empty())
         return false;
      block = std::move(m_blocks.front());
      m_blocks.pop_front();
      m_notFull.notify_one();
      return true;
   }

   void Close()
   {
      const std::lock_guard<std::mutex> lock(m_mutex);
      m_closed = true;
      m_notEmpty.notify_all();
      m_notFull.notify_all();
   }

   void Abort()
   {
      const std::lock_guard<std::mutex> lock(m_mutex);
      m_closed = true;
      m_blocks.clear();
      m_notEmpty.notify_all();
      m_notFull.notify_all();
   }
};

struct RefCache {
   std::string name;
   int id = -1;
   bool valid = false;

   int Lookup(RAMNTupleRefs &refs, const std::string &rname)
   {
      if (valid && rname == name)
         return id;
      id = refs.GetRefId(rname);
      name = rname;
      valid = true;
      return id;
   }
};

struct Progress {
   std::mutex mutex;
   std::condition_variable next_turn;
   size_t next_seq = 0; ///< The block whose clusters may be committed next.
   bool failed = false; ///< A worker or the reader gave up; everyone stops.
   FileOrder order;
   size_t records = 0;
   std::vector<uint64_t> qual_block_ends;
   std::vector<std::pair<std::string, std::string>> late_headers;

   void Fail()
   {
      const std::lock_guard<std::mutex> lock(mutex);
      failed = true;
      next_turn.notify_all();
   }

   bool Failed()
   {
      const std::lock_guard<std::mutex> lock(mutex);
      return failed;
   }
};

// A worker leaving by exception stops the others; the exception itself
// travels in the worker's future.
class StopOthersOnException {
   Progress &m_progress;
   BlockQueue &m_queue;
   const int m_exceptions = std::uncaught_exceptions();

public:
   StopOthersOnException(Progress &progress, BlockQueue &queue) : m_progress(progress), m_queue(queue) {}
   StopOthersOnException(const StopOthersOnException &) = delete;
   StopOthersOnException &operator=(const StopOthersOnException &) = delete;
   StopOthersOnException(StopOthersOnException &&) = delete;
   StopOthersOnException &operator=(StopOthersOnException &&) = delete;
   ~StopOthersOnException()
   {
      if (std::uncaught_exceptions() > m_exceptions) {
         m_progress.Fail();
         m_queue.Abort();
      }
   }
};

class CloseQueueOnExit {
   BlockQueue &m_queue;

public:
   explicit CloseQueueOnExit(BlockQueue &queue) : m_queue(queue) {}
   CloseQueueOnExit(const CloseQueueOnExit &) = delete;
   CloseQueueOnExit &operator=(const CloseQueueOnExit &) = delete;
   CloseQueueOnExit(CloseQueueOnExit &&) = delete;
   CloseQueueOnExit &operator=(CloseQueueOnExit &&) = delete;
   ~CloseQueueOnExit() { m_queue.Abort(); }
};

BlockOrder ProcessBlock(Block &block, QualityBlockWriter &out, uint32_t quality_policy, RefCache &rname_cache,
                        RefCache &rnext_cache, ramcore::SamRecord &sam_record, size_t &records,
                        std::vector<std::pair<std::string, std::string>> &late_headers)
{
   BlockOrder order;
   RAMNTupleRefs &rname_refs = *RAMNTupleRecord::GetRnameRefs();
   RAMNTupleRefs &rnext_refs = *RAMNTupleRecord::GetRnextRefs();

   std::vector<char> &data = block.data;
   size_t start = 0;
   size_t line_number = block.first_line;
   while (start < data.size()) {
      const auto nl = std::find(data.begin() + static_cast<std::ptrdiff_t>(start), data.end(), '\n');
      if (nl == data.end())
         break;
      const auto line_end = static_cast<size_t>(nl - data.begin());
      data[line_end] = '\0';
      char *line = &data[start];
      const size_t line_start = start;
      start = line_end + 1;
      ramcore::StripCRLF(line);
      const size_t this_line = line_number++;

      if (data[line_start] == '\0')
         continue;

      if (data[line_start] == '@') {
         const std::string_view header(line);
         const size_t tab = header.find('\t');
         if (tab != std::string_view::npos)
            late_headers.emplace_back(header.substr(0, tab), header.substr(tab + 1));
         else
            late_headers.emplace_back(header, "");
         continue;
      }

      sam_record.Clear();
      if (!ramcore::SamParser::ParseRecord(line, sam_record, this_line))
         continue;

      RAMNTupleRecord &rec = out.Record();
      FillRecordFields(sam_record, rec, quality_policy);
      rec.refid = rname_cache.Lookup(rname_refs, sam_record.rname);
      rec.refnext = rnext_cache.Lookup(rnext_refs, sam_record.rnext);

      order.max_span = std::max(order.max_span, rec.GetRefSpan());
      order.Note(rec.refid, rec.pos);
      out.Add();
      records++;
   }
   return order;
}

void WorkerMain(BlockQueue &queue, Progress &progress, const std::shared_ptr<ROOT::RNTupleFillContext> &ctx,
                uint32_t quality_policy)
{
   const StopOthersOnException guard(progress, queue);

   QualityBlockWriter out(ctx->CreateEntry(), ctx->CreateEntry(), [&ctx](ROOT::REntry &e) { ctx->Fill(e); });
   RefCache rname_cache;
   RefCache rnext_cache;
   ramcore::SamRecord sam_record;

   Block block;
   while (queue.Pop(block)) {
      size_t records = 0;
      std::vector<std::pair<std::string, std::string>> late_headers;
      const BlockOrder order =
         ProcessBlock(block, out, quality_policy, rname_cache, rnext_cache, sam_record, records, late_headers);
      // Quality blocks end with the input block, whose clusters are committed as a unit.
      out.Finish();
      const std::vector<uint64_t> block_ends = out.TakeBlockEnds();
      ctx->FlushCluster();

      std::unique_lock<std::mutex> lock(progress.mutex);
      progress.next_turn.wait(lock, [&] { return progress.next_seq == block.seq || progress.failed; });
      if (progress.failed)
         return;
      ctx->CommitStagedClusters();
      progress.order.Add(order);
      for (const uint64_t end : block_ends)
         progress.qual_block_ends.push_back(progress.records + end);
      progress.records += records;
      progress.late_headers.insert(progress.late_headers.end(), late_headers.begin(), late_headers.end());
      progress.next_seq++;
      progress.next_turn.notify_all();
   }
}

// Returns the offset of the first record line, npos if the block is all header.
size_t ConsumeHeader(std::vector<char> &data, TList &headers, size_t &lines)
{
   size_t offset = 0;
   while (offset < data.size()) {
      const auto nl = std::find(data.begin() + static_cast<std::ptrdiff_t>(offset), data.end(), '\n');
      size_t line_end = static_cast<size_t>(nl - data.begin());
      const size_t next = nl != data.end() ? line_end + 1 : data.size();
      while (line_end > offset && (data[line_end - 1] == '\r' || data[line_end - 1] == '\n'))
         --line_end;
      if (line_end == offset) {
         offset = next;
         lines++;
         continue;
      }
      if (data[offset] != '@')
         return offset;
      std::string line(data.begin() + static_cast<std::ptrdiff_t>(offset),
                       data.begin() + static_cast<std::ptrdiff_t>(line_end));
      const size_t tab = line.find('\t');
      if (tab != std::string::npos)
         HandleHeaderLine(headers, line.substr(0, tab), line.substr(tab + 1));
      else
         HandleHeaderLine(headers, line, "");
      offset = next;
      lines++;
   }
   return std::string::npos;
}

} // namespace

bool samtoramntuple(const char *datafile, const char *treefile, int compression_algorithm, uint32_t quality_policy,
                    int threads, size_t block_bytes)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   threads = std::max(threads, 1);

   std::unique_ptr<FILE, int (*)(FILE *)> input(fopen(datafile, "r"), fclose);
   if (!input) {
      printf("Failed to parse SAM file %s\n", datafile);
      return false;
   }

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (!rootFile || !rootFile->IsOpen()) {
      printf("Failed to create RAM file %s\n", treefile);
      return false;
   }

   ROOT::EnableThreadSafety();
   RAMNTupleRecord::InitializeRefs();

   TList headers;
   headers.SetName("headers");

   BlockReader reader(input.get(), block_bytes);
   size_t lines = 0;
   Block first;
   bool have_records = false;
   {
      std::vector<char> data;
      while (reader.Next(data)) {
         const size_t offset = ConsumeHeader(data, headers, lines);
         if (offset == std::string::npos) {
            data.clear();
            continue;
         }
         first.data.assign(data.begin() + static_cast<std::ptrdiff_t>(offset), data.end());
         have_records = true;
         break;
      }
   }

   // RNEXT ids in header order, so they do not depend on worker timing.
   {
      RAMNTupleRefs &rnext = *RAMNTupleRecord::GetRnextRefs();
      rnext.GetRefId("=");
      for (const auto &name : RAMNTupleRecord::GetRnameRefs()->GetRefs())
         rnext.GetRefId(name);
   }

   auto model = ROOT::RNTupleModel::CreateBare();
   model->MakeField<RAMNTupleRecord>("record");
   model->MakeField<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(compression_algorithm);
   writeOptions.SetMaxUnzippedPageSize(64000);
   // Required by RNTupleParallelWriter.
   writeOptions.SetUseBufferedWrite(true);

   Progress progress;
   {
      auto writer = ROOT::RNTupleParallelWriter::Append(std::move(model), "RAM", *rootFile, writeOptions);

      BlockQueue queue(static_cast<size_t>(threads));
      // Destruction order on an exception: close the queue, join the futures,
      // drop the contexts, then the writer.
      std::vector<std::shared_ptr<ROOT::RNTupleFillContext>> contexts;
      std::vector<std::future<void>> workers;
      const CloseQueueOnExit closer(queue);
      for (int i = 0; i < threads; i++) {
         auto ctx = writer->CreateFillContext();
         ctx->EnableStagedClusterCommitting();
         contexts.push_back(ctx);
         workers.push_back(
            std::async(std::launch::async, WorkerMain, std::ref(queue), std::ref(progress), ctx, quality_policy));
      }

      if (have_records) {
         size_t seq = 0;
         Block block = std::move(first);
         bool more = true;
         while (more) {
            block.seq = seq++;
            block.first_line = lines + 1;
            lines += static_cast<size_t>(std::count(block.data.begin(), block.data.end(), '\n'));
            queue.Push(std::move(block));
            block = Block{};
            more = !progress.Failed() && reader.Next(block.data);
         }
      }
      queue.Close();
      for (auto &w : workers)
         w.get();

      contexts.clear();
      writer.reset();
   }
   if (reader.Failed()) {
      printf("Failed to read SAM file %s\n", datafile);
      return false;
   }

   for (const auto &[tag, content] : progress.late_headers)
      HandleHeaderLine(headers, tag, content);

   RAMNTupleRecord::SetCoordinateSorted(progress.order.check.sorted);
   RAMNTupleRecord::NoteRefSpan(progress.order.max_span);
   RAMNTupleRecord::SetQualBlockEnds(std::move(progress.qual_block_ends));
   if (!progress.order.check.sorted)
      fprintf(stderr, "%s is not in coordinate order; region queries will read it in full.\n", datafile);
   RAMNTupleRecord::WriteAllRefs(*rootFile);

   headers.Write("headers", TObject::kSingleKey);
   rootFile->Close();

   printf("\nRAM file created: %s\n", treefile);
   printf("Number of entries: %zu\n", progress.records);

   RAMNTupleRecord::GetRnameRefs()->Print();
   RAMNTupleRecord::GetRnextRefs()->Print();

   printf("\nProcessed %d SAM headers\n", headers.GetSize());
   printf("Processed %zu SAM records with %d threads\n\n", progress.records, threads);

   stopwatch.Print();
   return true;
}
