#include "ramcore/BamtoNTuple.h"
#include "BlockPipeline.h"
#include "ramcore/RecordBlocks.h"

#include "rntuple/RAMNTupleRecord.h"
#include "rntuple/RAMTags.h"

#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TFile.h>
#include <TROOT.h>
#include <TStopwatch.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Conversion on threads, as from SAM (BlockPipeline.h): the calling thread
// reads the records into blocks, workers fill them into fill contexts of their
// own, and the clusters are committed in block order.

namespace {

using ramcore::blocks::Block;
using ramcore::blocks::BlockOrder;
using ramcore::blocks::Progress;
using ramcore::blocks::RefCache;
using ramcore::blocks::SharedRefs;

std::string GetSeq(const bam1_t *b)
{
   const int len = b->core.l_qseq;
   // A BAM record with no sequence is SEQ "*" in SAM, not an empty column.
   if (len <= 0)
      return "*";
   const uint8_t *enc = bam_get_seq(b);
   std::string s(static_cast<std::size_t>(len), 'N');
   for (int i = 0; i < len; ++i)
      s[static_cast<std::size_t>(i)] = seq_nt16_str[bam_seqi(enc, i)];
   return s;
}

std::string GetQual(const bam1_t *b)
{
   const int len = b->core.l_qseq;
   if (len <= 0)
      return "*";
   const uint8_t *q = bam_get_qual(b);
   if (q == nullptr || q[0] == 0xff)
      return "*";
   std::string s(static_cast<std::size_t>(len), '!');
   for (int i = 0; i < len; ++i)
      s[static_cast<std::size_t>(i)] = static_cast<char>(q[i] + 33);
   return s;
}

std::string GetCigar(const bam1_t *b)
{
   const uint32_t n = b->core.n_cigar;
   if (n == 0)
      return "*";
   const uint32_t *cigar = bam_get_cigar(b);
   static constexpr std::array<char, 9> ops = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
   std::string s;
   s.reserve(n * 4U);
   for (uint32_t i = 0; i < n; ++i) {
      s += std::to_string(bam_cigar_oplen(cigar[i]));
      s += ops.at(static_cast<std::size_t>(bam_cigar_op(cigar[i])));
   }
   return s;
}

std::string FormatTag(const char *name, char type_char, const std::string &value)
{
   std::string result;
   result.reserve(strlen(name) + 3 + value.size());
   result += name;
   result += ':';
   result += type_char;
   result += ':';
   result += value;
   return result;
}

std::vector<std::string> GetTags(const bam1_t *b)
{
   std::vector<std::string> tags;
   const uint8_t *aux = bam_get_aux(b);
   const uint8_t *end = b->data + b->l_data;

   while (aux < end) {
      const std::array<char, 3> tag = {static_cast<char>(aux[0]), static_cast<char>(aux[1]), '\0'};
      aux += 2;
      const char type = static_cast<char>(*aux++);

      switch (type) {
      case 'A':
         tags.emplace_back(FormatTag(tag.data(), 'A', std::string(1, static_cast<char>(*aux))));
         aux += 1;
         break;
      case 'c': {
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(static_cast<int8_t>(*aux))));
         aux += 1;
         break;
      }
      case 'C': {
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(*aux)));
         aux += 1;
         break;
      }
      case 's': {
         int16_t v{};
         memcpy(&v, aux, sizeof(v));
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(v)));
         aux += sizeof(v);
         break;
      }
      case 'S': {
         uint16_t v{};
         memcpy(&v, aux, sizeof(v));
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(v)));
         aux += sizeof(v);
         break;
      }
      case 'i': {
         int32_t v{};
         memcpy(&v, aux, sizeof(v));
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(v)));
         aux += sizeof(v);
         break;
      }
      case 'I': {
         uint32_t v{};
         memcpy(&v, aux, sizeof(v));
         tags.emplace_back(FormatTag(tag.data(), 'i', std::to_string(v)));
         aux += sizeof(v);
         break;
      }
      case 'f': {
         float v{};
         memcpy(&v, aux, sizeof(v));
         // Must match htslib's rendering exactly, or a converted record no
         // longer compares equal to the BAM it came from. htslib formats
         // float aux values with "%g" (6 significant digits, trailing zeros
         // stripped, scientific notation where shorter); std::to_string uses
         // "%f" and always emits 6 decimal places, turning htslib's
         // "pa:f:0.862" into "pa:f:0.862000".
         char buf[32];
         std::snprintf(buf, sizeof(buf), "%g", static_cast<double>(v));
         tags.emplace_back(FormatTag(tag.data(), 'f', buf));
         aux += sizeof(v);
         break;
      }
      case 'Z':
      case 'H': {
         const char *str = static_cast<const char *>(static_cast<const void *>(aux));
         tags.emplace_back(FormatTag(tag.data(), type, str));
         aux += strlen(str) + 1;
         break;
      }
      default: return tags;
      }
   }
   return tags;
}

void FillRecord(RAMNTupleRecord *rec, const bam1_t *b, const sam_hdr_t *hdr, uint32_t quality_policy, SharedRefs &refs,
                RefCache &rname_cache, RefCache &rnext_cache, RAMTagWriter &tags)
{
   const char *rname = (b->core.tid < 0) ? "*" : sam_hdr_tid2name(hdr, b->core.tid);

   const char *rnext = "*";
   if (b->core.mtid >= 0)
      rnext = (b->core.mtid == b->core.tid) ? "=" : sam_hdr_tid2name(hdr, b->core.mtid);

   rec->SetBit(quality_policy);
   rec->SetQNAME(bam_get_qname(b));
   rec->SetFLAG(b->core.flag);
   rec->refid = rname_cache.Lookup(refs.rnames, rname);
   rec->SetPOS(static_cast<int32_t>(b->core.pos + 1));
   rec->SetMAPQ(b->core.qual);
   rec->SetCIGAR(GetCigar(b));
   rec->refnext = rnext_cache.Lookup(refs.rnexts, rnext);
   rec->SetMate(static_cast<int32_t>(b->core.mpos + 1), static_cast<int32_t>(b->core.isize),
                RAMMateOnSameRef(rname, rnext));
   rec->SetSEQ(GetSeq(b));
   rec->SetQUAL(GetQual(b));
   tags.Split(GetTags(b), *rec);
}

// A block holds whole records, each as its core fields, the length of its
// data and the data, as sam_read1 left them.
void AppendRecord(std::vector<char> &data, const bam1_t *b)
{
   const std::size_t at = data.size();
   const auto length = static_cast<uint32_t>(b->l_data);
   data.resize(at + sizeof(b->core) + sizeof(length) + length);
   std::memcpy(&data[at], &b->core, sizeof(b->core));
   std::memcpy(&data[at + sizeof(b->core)], &length, sizeof(length));
   if (length > 0)
      std::memcpy(&data[at + sizeof(b->core) + sizeof(length)], b->data, length);
}

// Reads the record at \p at of \p data into \p b, whose data is copied into
// \p bytes so that it is aligned as htslib leaves it; returns where the next
// record starts.
std::size_t ReadRecord(const std::vector<char> &data, std::size_t at, bam1_t &b, std::vector<uint8_t> &bytes)
{
   std::memcpy(&b.core, &data[at], sizeof(b.core));
   at += sizeof(b.core);
   uint32_t length = 0;
   std::memcpy(&length, &data[at], sizeof(length));
   at += sizeof(length);
   bytes.resize(length);
   if (length > 0)
      std::memcpy(bytes.data(), &data[at], length);
   b.data = bytes.data();
   b.l_data = static_cast<int>(length);
   b.m_data = length;
   return at + length;
}

BlockOrder ProcessBlock(const Block &block, BlockWriter &out, const sam_hdr_t *hdr, uint32_t quality_policy,
                        SharedRefs &refs, RefCache &rname_cache, RefCache &rnext_cache, RAMTagWriter &tags, bam1_t &b,
                        std::vector<uint8_t> &bytes, size_t &records)
{
   BlockOrder order;
   for (std::size_t at = 0; at < block.data.size();) {
      at = ReadRecord(block.data, at, b, bytes);
      RAMNTupleRecord &rec = out.Record();
      FillRecord(&rec, &b, hdr, quality_policy, refs, rname_cache, rnext_cache, tags);
      order.max_span = std::max(order.max_span, rec.GetRefSpan());
      order.Note(rec.refid, rec.pos);
      order.counts.Note(rec.refid, rec.flag);
      // The mate's reference for flagstat, as htslib's mtid.
      int32_t mate = -1;
      if (b.core.mtid == b.core.tid)
         mate = rec.refid;
      else if (b.core.mtid >= 0)
         mate = refs.rnames.FindRefId(sam_hdr_tid2name(hdr, b.core.mtid));
      order.flags.Note(rec.flag, rec.refid, mate, rec.mapq);
      out.Add();
      records++;
   }
   return order;
}

void WorkerMain(ramcore::blocks::BlockQueue &queue, Progress &progress, ROOT::RNTupleFillContext &ctx,
                const sam_hdr_t *hdr, SharedRefs &refs, uint32_t quality_policy)
{
   RefCache rname_cache;
   RefCache rnext_cache;
   RAMTagWriter tags(refs.tags);
   bam1_t record{};
   std::vector<uint8_t> bytes;
   ramcore::blocks::FillBlocks(queue, progress, ctx,
                               [&](Block &block, BlockWriter &out, size_t &records, std::vector<std::string> &) {
                                  return ProcessBlock(block, out, hdr, quality_policy, refs, rname_cache, rnext_cache,
                                                      tags, record, bytes, records);
                               });
}

} // namespace

void bamtoramntuple(const char *bamfile, const char *treefile, bool /*split*/, bool /*cache*/,
                    int compression_algorithm, uint32_t quality_policy, int threads, size_t block_bytes)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   threads = std::max(threads, 1);

   samFile *bamIn = sam_open(bamfile, "r");
   if (bamIn == nullptr) {
      std::cerr << "Cannot open BAM: " << bamfile << "\n";
      return;
   }

   // htslib inflates the input on as many threads as fill the records.
   hts_set_threads(bamIn, threads);

   sam_hdr_t *hdr = sam_hdr_read(bamIn);
   if (hdr == nullptr) {
      std::cerr << "Cannot read BAM header\n";
      sam_close(bamIn);
      return;
   }

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (rootFile == nullptr || !rootFile->IsOpen()) {
      std::cerr << "Failed to create: " << treefile << "\n";
      sam_hdr_destroy(hdr);
      sam_close(bamIn);
      return;
   }

   ROOT::EnableThreadSafety();

   // The tables of the file being written, filled in as the records are
   // written and stored in METADATA at the end.
   SharedRefs refs;
   const int nRef = sam_hdr_nref(hdr);
   for (int i = 0; i < nRef; ++i)
      refs.rnames.GetRefId(sam_hdr_tid2name(hdr, i));
   // RNEXT ids in header order, so they do not depend on worker timing.
   refs.rnexts.GetRefId("=");
   for (const auto &name : refs.rnames.GetRefs())
      refs.rnexts.GetRefId(name);

   ROOT::RNTupleWriteOptions opts;
   opts.SetCompression(compression_algorithm);
   opts.SetMaxUnzippedPageSize(64000);
   // Required by RNTupleParallelWriter.
   opts.SetUseBufferedWrite(true);

   Progress progress;
   bool failed = false;
   {
      auto writer = ROOT::RNTupleParallelWriter::Append(ramcore::blocks::MakeBareModel(), "RAM", *rootFile, opts);
      const std::unique_ptr<bam1_t, void (*)(bam1_t *)> rec(bam_init1(), bam_destroy1);
      std::size_t count = 0;
      ramcore::blocks::FillInBlocks(
         *writer, threads, progress,
         [hdr, &refs, quality_policy](ramcore::blocks::BlockQueue &queue, Progress &p, ROOT::RNTupleFillContext &ctx) {
            WorkerMain(queue, p, ctx, hdr, refs, quality_policy);
         },
         [&](Block &block) {
            while (block.data.size() < block_bytes) {
               const int status = sam_read1(bamIn, hdr, rec.get());
               if (status < -1)
                  failed = true;
               if (status < 0)
                  break;
               AppendRecord(block.data, rec.get());
               if (++count % 1000000 == 0)
                  std::cerr << "Read " << count << " records...\n";
            }
            return !block.data.empty();
         });
      writer.reset();
   }
   if (failed) {
      std::cerr << "Failed to read BAM file " << bamfile << "\n";
      sam_hdr_destroy(hdr);
      sam_close(bamIn);
      return;
   }

   RAMFileMetadata meta;
   meta.max_ref_span = progress.order.max_span;
   meta.counts = std::move(progress.order.counts);
   meta.flags = progress.order.flags;
   meta.block_ends = std::move(progress.block_ends);
   // Region queries can only seek on a sorted file; the file records which it is.
   meta.coordinate_sorted = progress.order.check.sorted;
   if (!meta.coordinate_sorted)
      std::cerr << bamfile << " is not in coordinate order; region queries will read it in full.\n";

   const char *text = sam_hdr_str(hdr);
   if (text != nullptr) {
      std::istringstream ss(text);
      std::string line;
      while (std::getline(ss, line)) {
         if (line.size() >= 3)
            meta.header_lines.push_back(line);
      }
   }
   meta.rname_refs = refs.rnames.GetRefs();
   meta.rnext_refs = refs.rnexts.GetRefs();
   refs.tags.Store(meta);
   RAMNTupleRecord::WriteMetadata(*rootFile, meta);

   rootFile->Close();
   sam_hdr_destroy(hdr);
   sam_close(bamIn);

   std::cout << "\nRAM file created: " << treefile << "\n"
             << "Number of entries: " << progress.records << "\n";
   refs.rnames.Print();
   refs.rnexts.Print();

   stopwatch.Print();
}