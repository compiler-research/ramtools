#include "ramcore/BamtoNTuple.h"
#include "ramcore/QualityBlocks.h"

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>
#include <TStopwatch.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

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

void FillRecord(RAMNTupleRecord *rec, const bam1_t *b, const sam_hdr_t *hdr, uint32_t quality_policy,
                RAMNTupleRefs &rnames, RAMNTupleRefs &rnexts)
{
   const char *rname = (b->core.tid < 0) ? "*" : sam_hdr_tid2name(hdr, b->core.tid);

   const char *rnext = "*";
   if (b->core.mtid >= 0)
      rnext = (b->core.mtid == b->core.tid) ? "=" : sam_hdr_tid2name(hdr, b->core.mtid);

   rec->SetBit(quality_policy);
   rec->SetQNAME(bam_get_qname(b));
   rec->SetFLAG(b->core.flag);
   rec->refid = rnames.GetRefId(rname);
   rec->SetPOS(static_cast<int32_t>(b->core.pos + 1));
   rec->SetMAPQ(b->core.qual);
   rec->SetCIGAR(GetCigar(b));
   rec->refnext = rnexts.GetRefId(rnext);
   rec->SetPNEXT(static_cast<int32_t>(b->core.mpos + 1));
   rec->SetTLEN(static_cast<int32_t>(b->core.isize));
   rec->SetSEQ(GetSeq(b));
   rec->SetQUAL(GetQual(b));

   rec->ResetNOPT();
   for (const auto &tag : GetTags(b)) {
      rec->SetOPT(tag);
   }
}

} // namespace

void bamtoramntuple(const char *bamfile, const char *treefile, bool /*split*/, bool /*cache*/,
                    int compression_algorithm, uint32_t quality_policy)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   samFile *bamIn = sam_open(bamfile, "r");
   if (bamIn == nullptr) {
      std::cerr << "Cannot open BAM: " << bamfile << "\n";
      return;
   }

   hts_set_threads(bamIn, static_cast<int>(std::thread::hardware_concurrency()));

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

   // The state of the file being written: filled in as the records are written,
   // stored in METADATA at the end.
   RAMNTupleRefs rnames;
   RAMNTupleRefs rnexts;
   RAMFileMetadata meta;
   RAMCoordinateOrder order;

   const int nRef = sam_hdr_nref(hdr);
   for (int i = 0; i < nRef; ++i)
      rnames.GetRefId(sam_hdr_tid2name(hdr, i));

   auto model = RAMNTupleRecord::MakeModel();
   ROOT::RNTupleWriteOptions opts;
   opts.SetCompression(compression_algorithm);
   opts.SetMaxUnzippedPageSize(64000);

   auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, opts);
   QualityBlockWriter out(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                          [&writer](ROOT::REntry &e) { writer->Fill(e); });

   bam1_t *rec = bam_init1();
   std::size_t count = 0;

   while (sam_read1(bamIn, hdr, rec) >= 0) {
      RAMNTupleRecord &record = out.Record();
      FillRecord(&record, rec, hdr, quality_policy, rnames, rnexts);
      meta.max_ref_span = std::max(meta.max_ref_span, record.GetRefSpan());
      order.Note(record.refid, record.pos);
      meta.counts.Note(record.refid, record.flag);
      out.Add();

      ++count;
      if (count % 1000000 == 0)
         std::cerr << "Converted " << count << " records...\n";
   }

   bam_destroy1(rec);
   out.Finish();
   meta.qual_block_ends = out.TakeBlockEnds();
   writer.reset();

   // Region queries can only seek on a sorted file; the file records which it is.
   meta.coordinate_sorted = order.sorted;
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
   meta.rname_refs = rnames.GetRefs();
   meta.rnext_refs = rnexts.GetRefs();
   RAMNTupleRecord::WriteMetadata(*rootFile, meta);

   rootFile->Close();
   sam_hdr_destroy(hdr);
   sam_close(bamIn);

   std::cout << "\nRAM file created: " << treefile << "\n"
             << "Number of entries: " << count << "\n";
   rnames.Print();
   rnexts.Print();

   stopwatch.Print();
}