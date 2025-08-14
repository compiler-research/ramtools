//
// RAMNTupleRecord.cxx
// Complete implementation of RAM format using RNTuple

#include "rntuple/RAMNTupleRecord.h"
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TError.h>
#include <TFile.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>

using namespace ROOT::Experimental;

std::unique_ptr<RAMNTupleRefs> RAMNTupleRecord::fgRnameRefs = nullptr;
std::unique_ptr<RAMNTupleRefs> RAMNTupleRecord::fgRnextRefs = nullptr;
std::unique_ptr<RAMNTupleIndex> RAMNTupleRecord::fgIndex = nullptr;

static const char *kCodeToSeq = "=ACMGRSVTWYHKDBN";
static uint8_t kSeqToCode[256] = {0};
static bool kSeqTableInit = false;

// CIGAR encoding/decoding tables
static const char *kCodeToCigar = "MIDNSHP=X";
static uint8_t kCigarToCode[256] = {0};
static bool kCigarTableInit = false;

// Illumina 8-level quality binning: maps Q0-40+ to 8 values (0,1,6,15,22,27,33,37,40)
// Reduces quality data ~80% with minimal accuracy loss
const uint8_t RAMNTupleUtils::kIlluminaBinning[256] = {
   0,  1,  6,  6,  6,  6,  6,  6,  6,  6,  15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 22, 22, 22, 22, 22, 27, 27, 27,
   27, 27, 33, 33, 33, 33, 33, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};

// RAMNTupleRefs Implementation
RAMNTupleRefs::RAMNTupleRefs() : fLastId(-1)
{
   fRefVec.reserve(100);
}

int RAMNTupleRefs::GetRefId(const std::string &rname)
{
   if (rname == "*") {
      return -1;
   }

   if (rname == fLastName) {
      return fLastId;
   }

   auto it = std::find(fRefVec.begin(), fRefVec.end(), rname);
   if (it != fRefVec.end()) {
      fLastId = static_cast<int>(std::distance(fRefVec.begin(), it));
      fLastName = rname;
      return fLastId;
   }

   if (static_cast<int>(fRefVec.size()) >= static_cast<int>(fRefVec.capacity())) {
      fRefVec.reserve(fRefVec.capacity() * 2);
   }

   fRefVec.push_back(rname);
   fLastId = static_cast<int>(fRefVec.size() - 1);
   fLastName = rname;
   return fLastId;
}

const std::string &RAMNTupleRefs::GetRefName(int rid) const
{
   static const std::string star = "*";
   if (rid == -1) {
      return star;
   }
   if (rid < 0 || rid >= static_cast<int>(fRefVec.size())) {
      return star;
   }
   return fRefVec[rid];
}

void RAMNTupleRefs::Print() const
{
   int size = static_cast<int>(fRefVec.size());
   printf("RAMNTupleRefs vector:\n");
   for (int i = 0; i < size; i++) {
      printf("%d: %s\n", i, fRefVec[i].c_str());
   }
}
// RAMNTupleIndex Implementation
void RAMNTupleIndex::AddItem(int32_t refid, int32_t pos, int64_t row)
{
   fIndex.push_back({refid, pos, row});
   auto key = std::make_pair(refid, pos);
   fIndexMap[key] = row;
}

int64_t RAMNTupleIndex::GetRow(int32_t refid, int32_t pos) const
{
   if (fIndexMap.empty() && !fIndex.empty()) {
      const_cast<RAMNTupleIndex *>(this)->RebuildMap();
   }

   auto key = std::make_pair(refid, pos);
   auto low = fIndexMap.lower_bound(key);

   if (low == fIndexMap.end()) {
      return -1;
   } else if (low == fIndexMap.begin()) {
      return low->second;
   } else {
      if (low->first.first == refid && low->first.second == pos) {
         return low->second;
      } else {
         --low;
         return low->second;
      }
   }
}

void RAMNTupleIndex::RebuildMap() const
{
   fIndexMap.clear();
   for (const auto &entry : fIndex) {
      fIndexMap[{entry.refid, entry.pos}] = entry.entry;
   }
}

std::vector<int64_t> RAMNTupleIndex::GetRowsInRange(int32_t refid, int32_t start, int32_t end) const
{
   std::vector<int64_t> rows;

   for (const auto &entry : fIndex) {
      if (entry.refid == refid && entry.pos >= start && entry.pos <= end) {
         rows.push_back(entry.entry);
      }
   }

   return rows;
}

void RAMNTupleIndex::Print() const
{
   printf("RAMNTupleIndex map:\n");
   size_t count = 0;
   for (const auto &entry : fIndex) {
      printf("%lld: refid=%d, pos=%d\n", static_cast<long long>(entry.entry), entry.refid, entry.pos);
      if (++count > 10 && fIndex.size() > 20) {
         printf("... (%zu more entries)\n", fIndex.size() - 10);
         break;
      }
   }
}
// RAMNTupleRecord Implementation

RAMNTupleRecord::RAMNTupleRecord()
   : flag(0), refid(-1), pos(0), mapq(0), refnext(-1), pnext(0), tlen(0), compression_flags(kPhred33)
{
   InitializeRefs();
}

void RAMNTupleRecord::InitializeRefs()
{
   if (!fgRnameRefs)
      fgRnameRefs = std::make_unique<RAMNTupleRefs>();
   if (!fgRnextRefs)
      fgRnextRefs = std::make_unique<RAMNTupleRefs>();
   if (!fgIndex)
      fgIndex = std::make_unique<RAMNTupleIndex>();
}

std::unique_ptr<RNTupleReader> RAMNTupleRecord::OpenRAMFile(const std::string &filename, const std::string &ntupleName)
{

   InitializeRefs();

   try {
      auto reader = RNTupleReader::Open(ntupleName, filename);
      ReadAllRefs(filename);
      ReadIndex(filename);
      return reader;
   } catch (const std::exception &e) {
      ::Error("RAMNTupleRecord::OpenRAMFile", "Failed to open file: %s", e.what());
      return nullptr;
   }
}

void RAMNTupleRecord::WriteAllRefs(TFile &file)
{
   // Write refs as separate RNTuple in the same file
   if (!file.IsOpen())
      return;
   file.cd();

   // Create a simple model for metadata
   auto metaModel = RNTupleModel::Create();
   auto rnameField = metaModel->MakeField<std::vector<std::string>>("rname_refs");
   auto rnextField = metaModel->MakeField<std::vector<std::string>>("rnext_refs");

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(505);

   auto metaWriter = RNTupleWriter::Append(std::move(metaModel), "METADATA", file, writeOptions);
   auto metaEntry = metaWriter->CreateEntry();
   auto rnamePtr = metaEntry->GetPtr<std::vector<std::string>>("rname_refs");
   auto rnextPtr = metaEntry->GetPtr<std::vector<std::string>>("rnext_refs");

   *rnamePtr = fgRnameRefs->GetRefs();
   *rnextPtr = fgRnextRefs->GetRefs();
   metaWriter->Fill(*metaEntry);
}

void RAMNTupleRecord::ReadAllRefs(const std::string &filename)
{
   try {
      auto reader = RNTupleReader::Open("METADATA", filename);
      if (!reader || reader->GetNEntries() == 0)
         return;

      try {
         auto refs_view = reader->GetView<std::vector<std::string>>("rname_refs");
         const auto &refs = refs_view(0);
         fgRnameRefs->Clear();
         for (const auto &ref : refs) {
            fgRnameRefs->AddRef(ref);
         }
      } catch (...) {
         // Field doesn't exist
      }

      // Read next reference names
      try {
         auto rnext_view = reader->GetView<std::vector<std::string>>("rnext_refs");
         const auto &refs = rnext_view(0);
         fgRnextRefs->Clear();
         for (const auto &ref : refs) {
            fgRnextRefs->AddRef(ref);
         }
      } catch (...) {
         // Field doesn't exist
      }
   } catch (...) {
      // Metadata might not exist
   }
}

void RAMNTupleRecord::WriteIndex(TFile &file)
{
   if (!fgIndex || fgIndex->Size() == 0 || !file.IsOpen())
      return;
   file.cd();

   // Create index model
   auto indexModel = RNTupleModel::Create();
   auto indexField = indexModel->MakeField<std::vector<RAMNTupleIndex::IndexEntry>>("index_entries");

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(505);

   auto indexWriter = RNTupleWriter::Append(std::move(indexModel), "INDEX", file, writeOptions);
   auto indexEntry = indexWriter->CreateEntry();
   auto indexPtr = indexEntry->GetPtr<std::vector<RAMNTupleIndex::IndexEntry>>("index_entries");

   *indexPtr = fgIndex->GetEntries();
   indexWriter->Fill(*indexEntry);
}

void RAMNTupleRecord::ReadIndex(const std::string &filename)
{
   try {
      auto reader = RNTupleReader::Open("INDEX", filename);
      if (!reader || reader->GetNEntries() == 0)
         return;

      try {
         auto index_view = reader->GetView<std::vector<RAMNTupleIndex::IndexEntry>>("index_entries");
         const auto &entries = index_view(0);
         fgIndex->SetEntries(entries);
      } catch (...) {
         // Field doesn't exist
      }
   } catch (...) {
      // Index might not exist
   }
}

void RAMNTupleRecord::SetRNAME(const std::string &rname)
{
   refid = fgRnameRefs->GetRefId(rname);
}

void RAMNTupleRecord::SetRNEXT(const std::string &rnext)
{
   refnext = fgRnextRefs->GetRefId(rnext);
}

std::string RAMNTupleRecord::GetRNAME() const
{
   return fgRnameRefs->GetRefName(refid);
}

std::string RAMNTupleRecord::GetRNEXT() const
{
   return fgRnextRefs->GetRefName(refnext);
}

void RAMNTupleRecord::SetCIGAR(const std::string &cigar_str)
{
   cigar = RAMNTupleUtils::ParseCIGAR(cigar_str);
}

std::string RAMNTupleRecord::GetCIGAR() const
{
   return RAMNTupleUtils::FormatCIGAR(cigar);
}

void RAMNTupleRecord::SetSEQ(const std::string &seq_str)
{
   seq = RAMNTupleUtils::EncodeSequence(seq_str);
}

std::string RAMNTupleRecord::GetSEQ() const
{
   if (seq.size() < 4)
      return "";
   uint32_t length = *reinterpret_cast<const uint32_t *>(seq.data());
   return RAMNTupleUtils::DecodeSequence(seq.substr(4), length);
}

void RAMNTupleRecord::SetQUAL(const std::string &qual_str)
{
   qual = RAMNTupleUtils::EncodeQuality(qual_str, compression_flags);
}

std::string RAMNTupleRecord::GetQUAL() const
{
   return RAMNTupleUtils::DecodeQuality(qual, compression_flags);
}

int32_t RAMNTupleRecord::GetCIGAROPLEN(size_t idx) const
{
   if (idx >= cigar.size()) {
      ::Error("GetCIGAROPLEN", "idx=%zu out of range, max=%zu", idx, cigar.size());
      return 0;
   }
   return cigar[idx] >> 4;
}

int32_t RAMNTupleRecord::GetCIGAROP(size_t idx) const
{
   if (idx >= cigar.size()) {
      ::Error("GetCIGAROP", "idx=%zu out of range, max=%zu", idx, cigar.size());
      return 0;
   }
   return cigar[idx] & 0xf;
}

void RAMNTupleRecord::Print(const char *) const
{
   std::cout << GetQNAME() << "\t" << GetFLAG() << "\t" << GetRNAME() << "\t" << GetPOS() << "\t"
             << static_cast<int>(GetMAPQ()) << "\t" << GetCIGAR() << "\t" << GetRNEXT() << "\t" << GetPNEXT() << "\t"
             << GetTLEN() << "\t" << GetSEQ() << "\t" << GetQUAL();

   for (const auto &tag : tags) {
      std::cout << "\t" << tag;
   }
   std::cout << std::endl;
}

bool RAMNTupleRecord::IsValid() const
{
   return !qname.empty() && refid >= -1;
}

std::unique_ptr<RNTupleModel> RAMNTupleRecord::MakeModel()
{
   auto model = RNTupleModel::Create();

   model->MakeField<RAMNTupleRecord>("record");

   return model;
}
// Utility Functions Implementation

namespace RAMNTupleUtils {

void InitializeTables()
{
   if (!kSeqTableInit) {
      std::memset(kSeqToCode, 0, 256);
      for (int i = 1; i < 16; i++) {
         kSeqToCode[static_cast<uint8_t>(kCodeToSeq[i])] = i;
      }
      kSeqTableInit = true;
   }

   if (!kCigarTableInit) {
      std::memset(kCigarToCode, 0, 256);
      for (int i = 0; i < 9; i++) {
         kCigarToCode[static_cast<uint8_t>(kCodeToCigar[i])] = i;
      }
      kCigarTableInit = true;
   }
}

std::string EncodeSequence(const std::string &seq)
{
   InitializeTables();

   uint32_t length = seq.length();
   size_t encoded_size = 4 + (length + 1) / 2;
   std::string encoded;
   encoded.resize(encoded_size);

   std::memcpy(&encoded[0], &length, 4);

   size_t j = 4;
   for (size_t i = 0; i + 1 < length; i += 2) {
      encoded[j++] = (kSeqToCode[static_cast<uint8_t>(seq[i])] << 4) | kSeqToCode[static_cast<uint8_t>(seq[i + 1])];
   }
   if (length % 2) {
      encoded[j] = kSeqToCode[static_cast<uint8_t>(seq[length - 1])] << 4;
   }

   return encoded;
}

std::string DecodeSequence(const std::string &encoded_seq, size_t length)
{
   std::string seq;
   seq.resize(length);

   size_t pairs = length / 2;
   for (size_t i = 0; i < pairs; i++) {
      uint8_t byte = encoded_seq[i];
      seq[i * 2] = kCodeToSeq[byte >> 4];
      seq[i * 2 + 1] = kCodeToSeq[byte & 0xf];
   }
   if (length % 2) {
      seq[length - 1] = kCodeToSeq[encoded_seq[length / 2] >> 4];
   }

   return seq;
}

std::string EncodeQuality(const std::string &qual, uint32_t compression_flags)
{
   std::string encoded;

   if (compression_flags & RAMNTupleRecord::kDrop) {
      encoded = "*";
   } else if (compression_flags & RAMNTupleRecord::kIlluminaBinning) {
      encoded.resize(qual.size());
      for (size_t i = 0; i < qual.size(); i++) {
         encoded[i] = kIlluminaBinning[static_cast<uint8_t>(qual[i])];
      }
   } else {
      encoded = qual;
   }

   return encoded;
}

std::string DecodeQuality(const std::string &encoded_qual, uint32_t compression_flags)
{
   if (compression_flags & RAMNTupleRecord::kDrop) {
      return "*";
   }

   std::string qual = encoded_qual;

   if (compression_flags & RAMNTupleRecord::kIlluminaBinning) {
      for (auto &q : qual) {
         q += 33;
      }
   }

   return qual;
}

std::vector<uint32_t> ParseCIGAR(const std::string &cigar_str)
{
   InitializeTables();

   std::vector<uint32_t> cigar_ops;
   std::string num_str;

   for (char c : cigar_str) {
      if (std::isdigit(c)) {
         num_str += c;
      } else if (kCigarToCode[static_cast<uint8_t>(c)] < 9) {
         if (!num_str.empty()) {
            uint32_t len = std::stoul(num_str);
            uint8_t op = kCigarToCode[static_cast<uint8_t>(c)];
            cigar_ops.push_back((len << 4) | op);
            num_str.clear();
         }
      }
   }

   return cigar_ops;
}

std::string FormatCIGAR(const std::vector<uint32_t> &cigar_ops)
{
   std::ostringstream oss;

   for (uint32_t op : cigar_ops) {
      oss << (op >> 4) << kCodeToCigar[op & 0xf];
   }

   return oss.str();
}

} // namespace RAMNTupleUtils
// RAMNTupleConverter Implementation
void RAMNTupleConverter::ConvertSAMToRAMNTuple(const std::string &sam_file, const std::string &ram_file,
                                               uint32_t compression_flags)
{

   RAMNTupleRecord::InitializeRefs();

   auto file = std::unique_ptr<TFile>(TFile::Open(ram_file.c_str(), "RECREATE"));
   if (!file || !file->IsOpen()) {
      ::Error("ConvertSAMToRAMNTuple", "Cannot create RAM file: %s", ram_file.c_str());
      return;
   }

   auto model = RAMNTupleRecord::MakeModel();

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(505); // ZSTD level 5

   auto writer = RNTupleWriter::Append(std::move(model), "RAM", *file, writeOptions);
   auto defaultEntry = writer->CreateEntry();
   auto recordPtr = defaultEntry->GetPtr<RAMNTupleRecord>("record");

   std::ifstream sam(sam_file);
   if (!sam) {
      ::Error("ConvertSAMToRAMNTuple", "Cannot open SAM file: %s", sam_file.c_str());
      return;
   }

   std::string line;
   int64_t entry_number = 0;

   while (std::getline(sam, line)) {
      if (line.empty() || line[0] == '@') {

         if (line.substr(0, 3) == "@SQ") {
            size_t sn_pos = line.find("SN:");
            if (sn_pos != std::string::npos) {
               size_t tab_pos = line.find('\t', sn_pos);
               std::string seq_name =
                  line.substr(sn_pos + 3, tab_pos != std::string::npos ? tab_pos - sn_pos - 3 : std::string::npos);
               RAMNTupleRecord::GetRnameRefs()->GetRefId(seq_name);
            }
         }
         continue;
      }

      RAMNTupleRecord rec;
      rec.SetCompressionMode(compression_flags);

      std::istringstream iss(line);
      std::string qname, rname, cigar_str, rnext, seq_str, qual_str;
      int flag_int, pos_int, mapq_int, pnext_int, tlen_int;

      if (!(iss >> qname >> flag_int >> rname >> pos_int >> mapq_int >> cigar_str >> rnext >> pnext_int >> tlen_int >>
            seq_str >> qual_str)) {
         ::Warning("ConvertSAMToRAMNTuple", "Failed to parse line %ld", (long)(entry_number + 1));
         continue;
      }

      rec.SetQNAME(qname);
      rec.SetFLAG(static_cast<uint16_t>(flag_int));
      rec.SetRNAME(rname);
      rec.SetPOS(pos_int);
      rec.SetMAPQ(static_cast<uint8_t>(mapq_int));
      rec.SetCIGAR(cigar_str);
      rec.SetRNEXT(rnext);
      rec.SetPNEXT(pnext_int);
      rec.SetTLEN(tlen_int);
      rec.SetSEQ(seq_str);
      rec.SetQUAL(qual_str);

      std::string tag;
      while (iss >> tag) {
         rec.AddTag(tag);
      }

      if (rec.refid >= 0 && rec.pos >= 0) {
         RAMNTupleRecord::GetIndex()->AddItem(rec.refid, rec.pos, entry_number);
      }

      *recordPtr = std::move(rec);
      writer->Fill(*defaultEntry);

      entry_number++;

      if (entry_number % 100000 == 0) {
         std::cout << "Processed " << entry_number << " records..." << std::endl;
      }
   }

   writer.reset();

   RAMNTupleRecord::WriteAllRefs(*file);
   RAMNTupleRecord::WriteIndex(*file);

   std::cout << "Conversion complete. Processed " << entry_number << " records." << std::endl;
}

void RAMNTupleConverter::ConvertRAMNTupleToSAM(const std::string &ram_file, const std::string &sam_file)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(ram_file);
   if (!reader) {
      ::Error("ConvertRAMNTupleToSAM", "Cannot open RAM file: %s", ram_file.c_str());
      return;
   }

   std::ofstream sam(sam_file);
   if (!sam) {
      ::Error("ConvertRAMNTupleToSAM", "Cannot create SAM file: %s", sam_file.c_str());
      return;
   }

   sam << "@HD\tVN:1.6\tSO:unsorted\n";

   auto refs = RAMNTupleRecord::GetRnameRefs();
   if (refs) {
      for (size_t i = 0; i < refs->Size(); ++i) {
         sam << "@SQ\tSN:" << refs->GetRefName(i) << "\tLN:1\n";
      }
   }

   auto view = reader->GetView<RAMNTupleRecord>("record");

   for (auto i : reader->GetEntryRange()) {
      const auto &rec = view(i);
      rec.Print();
   }

   std::cout << "Conversion complete. Wrote " << reader->GetNEntries() << " records to SAM file." << std::endl;
}

void RAMNTupleConverter::BuildIndex(const std::string &ram_file)
{
   auto reader = RNTupleReader::Open("RAM", ram_file);
   if (!reader) {
      ::Error("BuildIndex", "Cannot open RAM file: %s", ram_file.c_str());
      return;
   }

   auto index = std::make_unique<RAMNTupleIndex>();
   auto view = reader->GetView<RAMNTupleRecord>("record");

   std::cout << "Building index..." << std::endl;

   for (auto i : reader->GetEntryRange()) {
      const auto &rec = view(i);

      if (rec.refid >= 0 && rec.pos >= 0) {
         index->AddItem(rec.refid, rec.pos, i);
      }

      if (i > 0 && i % 1000000 == 0) {
         std::cout << "Indexed " << i << " records..." << std::endl;
      }
   }

   std::string index_file = ram_file + ".idx";
   auto index_model = RNTupleModel::Create();
   auto index_field = index_model->MakeField<std::vector<RAMNTupleIndex::IndexEntry>>("index");

   auto writer = RNTupleWriter::Recreate(std::move(index_model), "INDEX", index_file);
   *index_field = index->GetEntries();
   writer->Fill();

   std::cout << "Index built with " << index->Size() << " entries." << std::endl;
}

void RAMNTupleConverter::ViewRegion(const std::string &ram_file, const std::string &region)
{
   std::string ref_name;
   int32_t start = 0, end = INT32_MAX;

   size_t colon_pos = region.find(':');
   if (colon_pos != std::string::npos) {
      ref_name = region.substr(0, colon_pos);

      size_t dash_pos = region.find('-', colon_pos);
      if (dash_pos != std::string::npos) {
         start = std::stoi(region.substr(colon_pos + 1, dash_pos - colon_pos - 1)) - 1;
         end = std::stoi(region.substr(dash_pos + 1)) - 1;
      } else {
         start = std::stoi(region.substr(colon_pos + 1)) - 1;
         end = start;
      }
   } else {
      ref_name = region;
   }

   auto reader = RAMNTupleRecord::OpenRAMFile(ram_file);
   if (!reader) {
      ::Error("ViewRegion", "Cannot open RAM file: %s", ram_file.c_str());
      return;
   }

   auto refs = RAMNTupleRecord::GetRnameRefs();
   int32_t ref_id = refs ? refs->GetRefId(ref_name) : -1;

   if (ref_id < 0) {
      ::Error("ViewRegion", "Unknown reference sequence: %s", ref_name.c_str());
      return;
   }

   auto index = RAMNTupleRecord::GetIndex();
   if (index && index->Size() > 0) {
      auto rows = index->GetRowsInRange(ref_id, start, end);
      auto view = reader->GetView<RAMNTupleRecord>("record");

      std::cout << "Found " << rows.size() << " reads in region " << region << std::endl;

      for (int64_t row : rows) {
         const auto &rec = view(row);
         rec.Print();
      }
   } else {
      std::cout << "No index found, performing linear scan..." << std::endl;

      auto view = reader->GetView<RAMNTupleRecord>("record");
      int count = 0;

      for (auto i : reader->GetEntryRange()) {
         const auto &rec = view(i);

         if (rec.refid == ref_id && rec.pos >= start && rec.pos <= end) {
            rec.Print();
            count++;
         }
      }

      std::cout << "Found " << count << " reads in region " << region << std::endl;
   }
}

