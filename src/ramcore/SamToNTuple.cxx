#include "ramcore/SamToNTuple.h"
#include "BlockPipeline.h"
#include "ramcore/RecordBlocks.h"
#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"
#include "rntuple/RAMTags.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TStopwatch.h>
#include <TFile.h>

#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <TROOT.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

void HandleHeaderLine(std::vector<std::string> &headers, RAMNTupleRefs &rnames, const std::string &line)
{
   headers.push_back(line);

   if (line.rfind("@SQ\t", 0) == 0) {
      size_t sn_pos = line.find("SN:", 4);
      if (sn_pos != std::string::npos) {
         sn_pos += 3;
         size_t tab_pos = line.find('\t', sn_pos);
         std::string ref_name =
            line.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
         rnames.GetRefId(ref_name);
      }
   }
}

// Everything but the two reference ids, which the caller resolves.
void FillRecordFields(const ramcore::SamRecord &sam_record, RAMNTupleRecord &rec, uint32_t quality_policy,
                      RAMTagWriter &tags)
{
   rec.SetBit(quality_policy);
   rec.SetQNAME(sam_record.qname);
   rec.SetFLAG(sam_record.flag);
   rec.SetPOS(sam_record.pos);
   rec.SetMAPQ(sam_record.mapq);
   rec.SetCIGAR(sam_record.cigar);
   rec.SetMate(sam_record.pnext, sam_record.tlen, RAMMateOnSameRef(sam_record.rname, sam_record.rnext));
   rec.SetSEQ(sam_record.seq);
   rec.SetQUAL(sam_record.qual);
   tags.Split(sam_record.optional_fields, rec);
}

// The mate's reference as htslib's mtid has it, for flagstat: the record's own
// for "=", -1 for "*" or a name the header lacks.
int32_t MateRefId(const RAMNTupleRefs &rnames, const std::string &rnext, int32_t refid)
{
   return rnext == "=" ? refid : rnames.FindRefId(rnext);
}

} // namespace

namespace {

// One output file, held open while the input streams past it.
struct ChromosomeWriter {
   std::unique_ptr<TFile> file{};
   std::unique_ptr<ROOT::RNTupleWriter> writer{};
   std::unique_ptr<BlockWriter> out;
   int64_t rows = 0;
   int32_t last_pos = -1;
   bool sorted = true;
   RAMRefCounts counts;
   RAMFlagCounts flags;
};

} // namespace

void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy)
{
   std::map<std::string, ChromosomeWriter> writers;
   std::vector<std::string> headers;
   // Shared by every output file: the ids are those of the whole input.
   RAMNTupleRefs rnames;
   RAMNTupleRefs rnexts;
   RAMTagTables tagTables;
   RAMTagWriter tags(tagTables);
   uint32_t max_span = 0;

   // The parser splits a header line at its first tab.
   auto header_callback = [&](const std::string &tag, const std::string &content) {
      HandleHeaderLine(headers, rnames, content.empty() ? tag : tag + "\t" + content);
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
      writeOptions.SetMaxUnzippedPageSize(RAMNTupleRecord::kMaxUnzippedPageSize);
      // Every chromosome's buffers are open at once, so keep the clusters small.
      writeOptions.SetApproxZippedClusterSize(8 * 1024 * 1024);

      cw.writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *cw.file, writeOptions);
      auto *writer = cw.writer.get();
      cw.out = std::make_unique<BlockWriter>(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                                             [writer](ROOT::REntry &e) { writer->Fill(e); });
      return cw;
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
      // A record with no reference has no chromosome file to go to.
      if (sam_record.rname == "*")
         return;

      ChromosomeWriter &cw = open_writer(sam_record.rname);
      RAMNTupleRecord &rec = cw.out->Record();

      FillRecordFields(sam_record, rec, quality_policy, tags);
      rec.refid = rnames.GetRefId(sam_record.rname);
      rec.refnext = rnexts.GetRefId(sam_record.rnext);

      max_span = std::max(max_span, rec.GetRefSpan());
      cw.counts.Note(rec.refid, rec.flag);
      cw.flags.Note(rec.flag, rec.refid, MateRefId(rnames, sam_record.rnext, rec.refid), rec.mapq);
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
   RAMFileMetadata meta;
   meta.header_lines = std::move(headers);
   meta.rname_refs = rnames.GetRefs();
   meta.rnext_refs = rnexts.GetRefs();
   tagTables.Store(meta);
   meta.max_ref_span = max_span;
   for (auto &[chr, cw] : writers) {
      cw.out->Finish();
      meta.block_ends = cw.out->TakeBlockEnds();
      cw.out.reset();
      cw.writer.reset();

      meta.coordinate_sorted = cw.sorted;
      meta.counts = std::move(cw.counts);
      meta.flags = cw.flags;
      if (!cw.sorted)
         fprintf(stderr, "%s: %s is not in coordinate order; region queries will read it in full.\n", datafile,
                 chr.c_str());
      RAMNTupleRecord::WriteMetadata(*cw.file, meta);
      cw.file->Close();

      printf("%s_%s.root: %lld records\n", output_prefix, chr.c_str(), static_cast<long long>(cw.rows));
   }
}

// Single-file conversion: the main thread reads blocks of whole lines, worker
// threads fill them into their own RNTupleFillContext and commit the staged
// clusters in block order (BlockPipeline.h).

namespace {

using ramcore::blocks::Block;
using ramcore::blocks::BlockOrder;
using ramcore::blocks::Progress;
using ramcore::blocks::RefCache;
using ramcore::blocks::SharedRefs;

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

BlockOrder ProcessBlock(Block &block, BlockWriter &out, uint32_t quality_policy, SharedRefs &refs,
                        RefCache &rname_cache, RefCache &rnext_cache, RAMTagWriter &tags,
                        ramcore::SamRecord &sam_record, size_t &records, std::vector<std::string> &late_headers)
{
   BlockOrder order;

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
         late_headers.emplace_back(line);
         continue;
      }

      sam_record.Clear();
      if (!ramcore::SamParser::ParseRecord(line, sam_record, this_line))
         continue;

      RAMNTupleRecord &rec = out.Record();
      FillRecordFields(sam_record, rec, quality_policy, tags);
      rec.refid = rname_cache.Lookup(refs.rnames, sam_record.rname);
      rec.refnext = rnext_cache.Lookup(refs.rnexts, sam_record.rnext);

      order.max_span = std::max(order.max_span, rec.GetRefSpan());
      order.Note(rec.refid, rec.pos);
      order.counts.Note(rec.refid, rec.flag);
      order.flags.Note(rec.flag, rec.refid, MateRefId(refs.rnames, sam_record.rnext, rec.refid), rec.mapq);
      out.Add();
      records++;
   }
   return order;
}

void WorkerMain(ramcore::blocks::BlockQueue &queue, Progress &progress, ROOT::RNTupleFillContext &ctx, SharedRefs &refs,
                uint32_t quality_policy)
{
   RefCache rname_cache;
   RefCache rnext_cache;
   RAMTagWriter tags(refs.tags);
   ramcore::SamRecord sam_record;
   ramcore::blocks::FillBlocks(
      queue, progress, ctx,
      [&](Block &block, BlockWriter &out, size_t &records, std::vector<std::string> &late_headers) {
         return ProcessBlock(block, out, quality_policy, refs, rname_cache, rnext_cache, tags, sam_record, records,
                             late_headers);
      });
}

// Returns the offset of the first record line, npos if the block is all header.
size_t ConsumeHeader(std::vector<char> &data, std::vector<std::string> &headers, RAMNTupleRefs &rnames, size_t &lines)
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
      HandleHeaderLine(headers, rnames,
                       std::string(data.begin() + static_cast<std::ptrdiff_t>(offset),
                                   data.begin() + static_cast<std::ptrdiff_t>(line_end)));
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

   std::vector<std::string> headers;
   SharedRefs refs;

   BlockReader reader(input.get(), block_bytes);
   size_t lines = 0;
   Block first;
   bool have_records = false;
   {
      std::vector<char> data;
      while (reader.Next(data)) {
         const size_t offset = ConsumeHeader(data, headers, refs.rnames, lines);
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
   refs.rnexts.GetRefId("=");
   for (const auto &name : refs.rnames.GetRefs())
      refs.rnexts.GetRefId(name);

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(compression_algorithm);
   writeOptions.SetMaxUnzippedPageSize(RAMNTupleRecord::kMaxUnzippedPageSize);
   // Required by RNTupleParallelWriter.
   writeOptions.SetUseBufferedWrite(true);

   Progress progress;
   {
      auto writer =
         ROOT::RNTupleParallelWriter::Append(ramcore::blocks::MakeBareModel(), "RAM", *rootFile, writeOptions);
      bool first_pending = have_records;
      ramcore::blocks::FillInBlocks(
         *writer, threads, progress,
         [&refs, quality_policy](ramcore::blocks::BlockQueue &queue, Progress &p, ROOT::RNTupleFillContext &ctx) {
            WorkerMain(queue, p, ctx, refs, quality_policy);
         },
         [&](Block &block) {
            if (first_pending) {
               block.data = std::move(first.data);
               first_pending = false;
            } else if (!have_records || !reader.Next(block.data)) {
               return false;
            }
            block.first_line = lines + 1;
            lines += static_cast<size_t>(std::count(block.data.begin(), block.data.end(), '\n'));
            return true;
         });
      writer.reset();
   }
   if (reader.Failed()) {
      printf("Failed to read SAM file %s\n", datafile);
      return false;
   }

   for (const auto &line : progress.late_headers)
      HandleHeaderLine(headers, refs.rnames, line);
   const size_t header_lines = headers.size();

   RAMFileMetadata meta;
   meta.header_lines = std::move(headers);
   meta.rname_refs = refs.rnames.GetRefs();
   meta.rnext_refs = refs.rnexts.GetRefs();
   refs.tags.Store(meta);
   meta.max_ref_span = progress.order.max_span;
   meta.coordinate_sorted = progress.order.check.sorted;
   meta.block_ends = std::move(progress.block_ends);
   meta.counts = std::move(progress.order.counts);
   meta.flags = progress.order.flags;
   if (!meta.coordinate_sorted)
      fprintf(stderr, "%s is not in coordinate order; region queries will read it in full.\n", datafile);
   RAMNTupleRecord::WriteMetadata(*rootFile, meta);
   rootFile->Close();

   printf("\nRAM file created: %s\n", treefile);
   printf("Number of entries: %zu\n", progress.records);

   refs.rnames.Print();
   refs.rnexts.Print();

   printf("\nProcessed %zu SAM headers\n", header_lines);
   printf("Processed %zu SAM records with %d threads\n\n", progress.records, threads);

   stopwatch.Print();
   return true;
}
