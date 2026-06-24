#include "ramcore/BamtoNTuple.h"

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>
#include <TStopwatch.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>

namespace {

constexpr uint16_t kUnmapped = 0x4;
constexpr int32_t kPositionInterval = 10000;
constexpr int64_t kMappedInterval = 100;

std::string GetSeq(const bam1_t *b)
{
   const int len = b->core.l_qseq;
   const uint8_t *enc = bam_get_seq(b);
   std::string s(static_cast<std::size_t>(len), 'N');
   for (int i = 0; i < len; ++i)
      s[static_cast<std::size_t>(i)] = seq_nt16_str[bam_seqi(enc, i)];
   return s;
}

std::string GetQual(const bam1_t *b)
{
   const int len = b->core.l_qseq;
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
         tags.emplace_back(FormatTag(tag.data(), 'f', std::to_string(v)));
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

void FillRecord(RAMNTupleRecord *rec, const bam1_t *b, const sam_hdr_t *hdr, uint32_t quality_policy)
{
   const char *rname = (b->core.tid < 0) ? "*" : sam_hdr_tid2name(hdr, b->core.tid);

   const char *rnext = "*";
   if (b->core.mtid >= 0)
      rnext = (b->core.mtid == b->core.tid) ? "=" : sam_hdr_tid2name(hdr, b->core.mtid);

   rec->SetBit(quality_policy);
   rec->SetQNAME(bam_get_qname(b));
   rec->SetFLAG(b->core.flag);
   rec->SetREFID(rname);
   rec->SetPOS(static_cast<int32_t>(b->core.pos + 1));
   rec->SetMAPQ(b->core.qual);
   rec->SetCIGAR(GetCigar(b));
   rec->SetREFNEXT(rnext);
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

void bamtoramntuple(const char *bamfile, const char *treefile, bool index, bool /*split*/, bool /*cache*/,
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

   RAMNTupleRecord::InitializeRefs();

   const int nRef = sam_hdr_nref(hdr);
   for (int i = 0; i < nRef; ++i)
      RAMNTupleRecord::GetRnameRefs()->SetRefId(sam_hdr_tid2name(hdr, i));

   auto model = RAMNTupleRecord::MakeModel();
   ROOT::RNTupleWriteOptions opts;
   opts.SetCompression(compression_algorithm);
   opts.SetMaxUnzippedPageSize(64000);

   auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, opts);
   auto entry = writer->GetModel().CreateEntry();
   auto recordPtr = entry->GetPtr<RAMNTupleRecord>("record");

   bam1_t *rec = bam_init1();
   std::size_t count = 0;
   int64_t mapped_count = 0;
   int32_t last_refid = -1;
   int32_t last_indexed_pos = -kPositionInterval;

   while (sam_read1(bamIn, hdr, rec) >= 0) {
      FillRecord(recordPtr.get(), rec, hdr, quality_policy);
      writer->Fill(*entry);

      if (index && !(rec->core.flag & kUnmapped) && recordPtr->GetREFID() >= 0) {
         int32_t current_refid = recordPtr->GetREFID();
         int32_t current_pos = recordPtr->GetPOS() - 1;

         bool new_chrom = (current_refid != last_refid);
         bool far_enough = (current_pos - last_indexed_pos >= kPositionInterval);
         bool periodic = (mapped_count % kMappedInterval == 0);
         bool duplicate = (!new_chrom && current_pos == last_indexed_pos);

         if ((new_chrom || far_enough || periodic) && !duplicate) {
            RAMNTupleRecord::GetIndex()->AddItem(current_refid, current_pos, static_cast<int64_t>(count));
            last_refid = current_refid;
            last_indexed_pos = current_pos;
         }
         mapped_count++;
      }

      ++count;
      if (count % 1000000 == 0)
         std::cerr << "Converted " << count << " records...\n";
   }

   bam_destroy1(rec);
   writer.reset();

   if (index)
      RAMNTupleRecord::WriteIndex(*rootFile);
   RAMNTupleRecord::WriteAllRefs(*rootFile);

   TList headers;
   headers.SetName("headers");
   const char *text = sam_hdr_str(hdr);
   if (text != nullptr) {
      std::istringstream ss(text);
      std::string line;
      while (std::getline(ss, line)) {
         if (line.size() >= 3) {
            auto named_obj = std::make_unique<TNamed>(line.substr(0, 3).c_str(), line.c_str());
            headers.Add(named_obj.release());
         }
      }
   }
   headers.Write();

   rootFile->Close();
   sam_hdr_destroy(hdr);
   sam_close(bamIn);

   std::cout << "\nRAM file created: " << treefile << "\n"
             << "Number of entries: " << count << "\n";
   RAMNTupleRecord::GetRnameRefs()->Print();
   RAMNTupleRecord::GetRnextRefs()->Print();
   if (index)
      std::cout << "Index entries: " << RAMNTupleRecord::GetIndex()->Size() << "\n";

   stopwatch.Print();
}
