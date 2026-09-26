//
// RAMNTupleRecord.cxx
// Complete implementation of RAM format using RNTuple

#include "rntuple/RAMNTupleRecord.h"
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TError.h>
#include <TFile.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

using namespace ROOT;

std::unique_ptr<RAMNTupleRefs> RAMNTupleRecord::fgRnameRefs = nullptr;
std::unique_ptr<RAMNTupleRefs> RAMNTupleRecord::fgRnextRefs = nullptr;
uint32_t RAMNTupleRecord::fgMaxRefSpan = 0;
RAMCoordinateOrder RAMNTupleRecord::fgOrder{};
std::vector<uint64_t> RAMNTupleRecord::fgQualBlockEnds{};

static constexpr std::array<char, 16> kCodeToSeq{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                                 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
static uint8_t kSeqToCode[256] = {0};

// CIGAR encoding/decoding tables
static constexpr std::array<char, 9> kCodeToCigar{'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'};
static uint8_t kCigarToCode[256] = {0};

// Any byte that is not a base code maps here. 15 is 'N', as in htslib's
// seq_nt16_table, so an unrepresentable base degrades to "unknown" rather than
// to '=' ("identical to the reference").
static constexpr uint8_t kSeqCodeUnknown = 15;

// Marks a byte that is not one of MIDNSHP=X; 0 would look like a valid 'M'.
static constexpr uint8_t kCigarCodeInvalid = 0xFF;

// (length << 4) | opcode leaves 28 bits for the length.
static constexpr uint32_t kMaxCigarOpLen = 0x0FFFFFFF;

// The 4-byte length prefix of a packed sequence, stored little-endian so the
// byte order is part of the format and reads need no aligned load.
static void StoreLE32(char *dst, uint32_t value)
{
   const std::array<unsigned char, 4> bytes{
      static_cast<unsigned char>(value & 0xFF),
      static_cast<unsigned char>((value >> 8) & 0xFF),
      static_cast<unsigned char>((value >> 16) & 0xFF),
      static_cast<unsigned char>((value >> 24) & 0xFF),
   };
   std::memcpy(dst, bytes.data(), bytes.size());
}

static uint32_t LoadLE32(const char *src)
{
   std::array<unsigned char, 4> bytes{};
   std::memcpy(bytes.data(), src, bytes.size());
   return static_cast<uint32_t>(bytes[0]) | (static_cast<uint32_t>(bytes[1]) << 8) |
          (static_cast<uint32_t>(bytes[2]) << 16) | (static_cast<uint32_t>(bytes[3]) << 24);
}

// Illumina 8-level quality binning: maps Q0-40+ to 8 values (0,1,6,15,22,27,33,37,40)
// Reduces quality data ~80% with minimal accuracy loss
const uint8_t RAMNTupleUtils::kIlluminaBinning[256] = {
   0,  1,  6,  6,  6,  6,  6,  6,  6,  6,  15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 22, 22, 22, 22, 22, 27, 27, 27,
   27, 27, 33, 33, 33, 33, 33, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
   40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40};

// RAMNTupleRefs Implementation
RAMNTupleRefs::RAMNTupleRefs()
{
   m_refVec.reserve(100);
}

void RAMNTupleRefs::Rebuild()
{
   m_index.clear();
   for (size_t i = 0; i < m_refVec.size(); i++)
      m_index.emplace(m_refVec[i], static_cast<int>(i));
}

int RAMNTupleRefs::GetRefId(const std::string &rname)
{
   if (rname == "*") {
      return -1;
   }

   const std::lock_guard<std::mutex> lock(m_mutex);
   auto it = m_index.find(rname);
   if (it != m_index.end())
      return it->second;

   m_refVec.push_back(rname);
   const int id = static_cast<int>(m_refVec.size() - 1);
   m_index.emplace(rname, id);
   return id;
}

int RAMNTupleRefs::FindRefId(const std::string &rname) const
{
   if (rname == "*")
      return -1;
   const std::lock_guard<std::mutex> lock(m_mutex);
   auto it = m_index.find(rname);
   return it != m_index.end() ? it->second : -1;
}

const std::string &RAMNTupleRefs::GetRefName(int rid) const
{
   static const std::string star = "*";
   const std::lock_guard<std::mutex> lock(m_mutex);
   if (rid < 0 || rid >= static_cast<int>(m_refVec.size())) {
      return star;
   }
   return m_refVec[rid];
}

size_t RAMNTupleRefs::Size() const
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   return m_refVec.size();
}

void RAMNTupleRefs::Clear()
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   m_refVec.clear();
   m_index.clear();
}

void RAMNTupleRefs::AddRef(const std::string &ref)
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   m_refVec.push_back(ref);
   m_index.emplace(ref, static_cast<int>(m_refVec.size() - 1));
}

std::vector<std::string> RAMNTupleRefs::GetRefs() const
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   return m_refVec;
}

void RAMNTupleRefs::SetRefs(const std::vector<std::string> &refs)
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   m_refVec = refs;
   Rebuild();
}

void RAMNTupleRefs::Print() const
{
   const std::lock_guard<std::mutex> lock(m_mutex);
   int size = static_cast<int>(m_refVec.size());
   printf("RAMNTupleRefs vector:\n");
   for (int i = 0; i < size; i++) {
      printf("%d: %s\n", i, m_refVec[i].c_str());
   }
}
// RAMNTupleRecord Implementation

RAMNTupleRecord::RAMNTupleRecord()
   : flag(0), refid(-1), pos(0), mapq(0), refnext(-1), pnext(0), tlen(0), compression_flags(kPhred33)
{
   EnsureTables();
}

void RAMNTupleRecord::EnsureTables()
{
   if (!fgRnameRefs)
      fgRnameRefs = std::make_unique<RAMNTupleRefs>();
   if (!fgRnextRefs)
      fgRnextRefs = std::make_unique<RAMNTupleRefs>();
}

// Resets the per-file state. Only the writers and OpenRAMFile() may call this:
// RNTuple constructs a record whenever a view or a writer model is created,
// so the constructor must not.
void RAMNTupleRecord::InitializeRefs()
{
   EnsureTables();
   fgMaxRefSpan = 0;
   fgOrder = RAMCoordinateOrder{};
   fgQualBlockEnds.clear();
}

// Whether the file holds an RNTuple of that name; false for a missing file.
static bool HasNTuple(const std::string &filename, const std::string &ntupleName)
{
   std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
   return file && !file->IsZombie() && file->Get<ROOT::RNTuple>(ntupleName.c_str()) != nullptr;
}

std::unique_ptr<RNTupleReader> RAMNTupleRecord::OpenRAMFile(const std::string &filename, const std::string &ntupleName)
{
   InitializeRefs();

   if (!HasNTuple(filename, ntupleName)) {
      ::Error("RAMNTupleRecord::OpenRAMFile", "%s has no RNTuple %s", filename.c_str(), ntupleName.c_str());
      return nullptr;
   }
   auto reader = RNTupleReader::Open(ntupleName, filename);
   ReadAllRefs(filename);
   return reader;
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
   auto spanField = metaModel->MakeField<uint32_t>("max_ref_span");
   auto sortedField = metaModel->MakeField<bool>("coordinate_sorted");
   auto blocksField = metaModel->MakeField<std::vector<uint64_t>>("qual_block_ends");

   RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(505);

   auto metaWriter = RNTupleWriter::Append(std::move(metaModel), "METADATA", file, writeOptions);
   auto metaEntry = metaWriter->GetModel().CreateEntry();
   auto rnamePtr = metaEntry->GetPtr<std::vector<std::string>>("rname_refs");
   auto rnextPtr = metaEntry->GetPtr<std::vector<std::string>>("rnext_refs");

   auto spanPtr = metaEntry->GetPtr<uint32_t>("max_ref_span");
   auto sortedPtr = metaEntry->GetPtr<bool>("coordinate_sorted");
   auto blocksPtr = metaEntry->GetPtr<std::vector<uint64_t>>("qual_block_ends");

   *rnamePtr = fgRnameRefs->GetRefs();
   *rnextPtr = fgRnextRefs->GetRefs();
   *spanPtr = fgMaxRefSpan;
   *sortedPtr = fgOrder.sorted;
   *blocksPtr = fgQualBlockEnds;
   metaWriter->Fill(*metaEntry);
}

void RAMNTupleRecord::ReadAllRefs(const std::string &filename)
{
   if (!HasNTuple(filename, "METADATA"))
      return;
   auto reader = RNTupleReader::Open("METADATA", filename);
   if (reader->GetNEntries() == 0)
      return;
   const auto &desc = reader->GetDescriptor();
   auto has = [&](const char *field) { return desc.FindFieldId(field) != ROOT::kInvalidDescriptorId; };

   if (has("rname_refs"))
      fgRnameRefs->SetRefs(reader->GetView<std::vector<std::string>>("rname_refs")(0));
   if (has("rnext_refs"))
      fgRnextRefs->SetRefs(reader->GetView<std::vector<std::string>>("rnext_refs")(0));

   // Absent in files written before the fields existed: 0 means "unknown", and
   // such files are read as sorted.
   fgMaxRefSpan = has("max_ref_span") ? reader->GetView<uint32_t>("max_ref_span")(0) : 0;
   fgOrder.sorted = has("coordinate_sorted") ? reader->GetView<bool>("coordinate_sorted")(0) : true;
   if (has("qual_block_ends"))
      fgQualBlockEnds = reader->GetView<std::vector<uint64_t>>("qual_block_ends")(0);
}

void RAMNTupleRecord::SetRNAME(const std::string &rname)
{
   refid = fgRnameRefs->GetRefId(rname);
}

void RAMNTupleRecord::SetRNEXT(const std::string &rnext)
{
   refnext = fgRnextRefs->GetRefId(rnext);
}

const std::string &RAMNTupleRecord::GetRNAME() const
{
   return fgRnameRefs->GetRefName(refid);
}

const std::string &RAMNTupleRecord::GetRNEXT() const
{
   return fgRnextRefs->GetRefName(refnext);
}

uint32_t RAMNTupleRecord::GetRefSpan() const
{
   uint32_t span = 0;
   for (uint32_t op : cigar) {
      switch (op & 0xF) {
      case RAM_CIGAR_M:
      case RAM_CIGAR_D:
      case RAM_CIGAR_N:
      case RAM_CIGAR_EQUAL:
      case RAM_CIGAR_X: span += (op >> 4); break;
      default: break;
      }
   }
   return span;
}

int RAMNTupleRecord::GetSEQLEN() const
{
   if (seq.size() < 4)
      return 0;
   return static_cast<int>(LoadLE32(seq.data()));
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
   // Restores the "*" that EncodeSequence folded into an empty payload.
   if (seq.size() < 4)
      return "*";
   const uint32_t length = LoadLE32(seq.data());
   return RAMNTupleUtils::DecodeSequence(seq.data() + 4, seq.size() - 4, length);
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
   model->MakeField<std::vector<std::uint8_t>>(kQualBlockField);

   return model;
}
// Utility Functions Implementation

namespace RAMNTupleUtils {

void InitializeTables()
{
   static const bool initialised = [] {
      std::memset(kSeqToCode, kSeqCodeUnknown, 256);
      for (size_t i = 0; i < kCodeToSeq.size(); i++) {
         const char base = kCodeToSeq[i];
         kSeqToCode[static_cast<uint8_t>(base)] = static_cast<uint8_t>(i);
         // SAM permits lowercase bases (SEQ is [A-Za-z=.]+).
         kSeqToCode[static_cast<uint8_t>(std::tolower(static_cast<unsigned char>(base)))] = static_cast<uint8_t>(i);
      }
      std::memset(kCigarToCode, kCigarCodeInvalid, 256);
      for (size_t i = 0; i < kCodeToCigar.size(); i++) {
         kCigarToCode[static_cast<uint8_t>(kCodeToCigar[i])] = static_cast<uint8_t>(i);
      }
      return true;
   }();
   (void)initialised;
}

std::string EncodeSequence(const std::string &seq)
{
   InitializeTables();

   // "*" is SAM's "sequence not stored" sentinel, not a one-base read. An empty
   // packed string means "no sequence"; a genuinely empty SEQ still gets its
   // 4-byte length prefix.
   if (seq == "*") {
      return {};
   }

   const uint32_t length = static_cast<uint32_t>(seq.length());
   const size_t encoded_size = 4 + ((static_cast<size_t>(length) + 1) / 2);
   std::string encoded;
   encoded.resize(encoded_size);

   StoreLE32(&encoded[0], length);

   size_t j = 4;
   for (size_t i = 0; i + 1 < length; i += 2) {
      encoded[j++] = static_cast<char>((kSeqToCode[static_cast<uint8_t>(seq[i])] << 4) |
                                       kSeqToCode[static_cast<uint8_t>(seq[i + 1])]);
   }
   if (length % 2) {
      encoded[j] = static_cast<char>(kSeqToCode[static_cast<uint8_t>(seq[length - 1])] << 4);
   }

   return encoded;
}

// NOLINTNEXTLINE(misc-use-internal-linkage) -- declared in RAMNTupleRecord.h
std::string DecodeSequence(const char *packed, size_t packed_size, size_t length)
{
   InitializeTables();

   // A truncated payload would otherwise be read past its end and yield bases
   // that were never stored.
   const size_t needed = (length + 1) / 2;
   if (packed == nullptr || packed_size < needed) {
      ::Error("DecodeSequence", "packed sequence holds %zu bytes, %zu needed for %zu bases", packed_size, needed,
              length);
      return {};
   }

   std::string seq;
   seq.resize(length);

   const std::string_view packed_bytes(packed, packed_size);
   const size_t pairs = length / 2;
   for (size_t i = 0; i < pairs; i++) {
      const uint8_t byte = static_cast<uint8_t>(packed_bytes[i]);
      seq[i * 2] = kCodeToSeq[byte >> 4];
      seq[i * 2 + 1] = kCodeToSeq[byte & 0xf];
   }
   if (length % 2) {
      seq[length - 1] = kCodeToSeq[static_cast<uint8_t>(packed_bytes[length / 2]) >> 4];
   }

   return seq;
}

std::string EncodeQuality(const std::string &qual, uint32_t compression_flags)
{
   if (compression_flags & RAMNTupleRecord::kDrop) {
      return "*";
   }

   if (compression_flags & RAMNTupleRecord::kIlluminaBinning) {
      // "*" means "quality not available". It is a sentinel, not a Phred
      // string, so it must not be fed through the binning table
      if (qual == "*")
         return {};

      std::string encoded(qual.size(), '\0');
      for (size_t i = 0; i < qual.size(); i++) {
         // SAM stores quality as Phred+33 ASCII, but kIlluminaBinning is
         // indexed by the Phred VALUE. Without the -33 every lookup lands 33
         // slots too far right
         // Clamp to 0..93: 93 is SAM's maximum Phred, and it also keeps the
         // index inside the initialised part of the table (entries 110..255
         // are zero-filled, so an out-of-range value would silently decode as
         // Q0 -- the opposite error, but still an error).
         const int phred = std::clamp(static_cast<int>(static_cast<unsigned char>(qual[i])) - 33, 0, 93);
         encoded[i] = static_cast<char>(kIlluminaBinning[phred]);
      }
      return encoded;
   }

   return qual;
}

std::string DecodeQuality(const std::string &encoded_qual, uint32_t compression_flags)
{
   if (compression_flags & RAMNTupleRecord::kDrop) {
      return "*";
   }

   if (compression_flags & RAMNTupleRecord::kIlluminaBinning) {
      // Empty is the "quality not available" sentinel written by
      // EncodeQuality; restore the "*" it stood for.
      if (encoded_qual.empty())
         return "*";

      std::string qual = encoded_qual;
      for (auto &q : qual) {
         q = static_cast<char>(static_cast<unsigned char>(q) + 33);
      }
      return qual;
   }

   return encoded_qual;
}

std::vector<uint32_t> ParseCIGAR(const std::string &cigar_str)
{
   InitializeTables();

   std::vector<uint32_t> cigar_ops;

   // "*" (and an empty field) mean the alignment has no CIGAR.
   if (cigar_str.empty() || cigar_str == "*")
      return cigar_ops;

   uint64_t length = 0;
   bool have_length = false;

   auto fail = [&]() {
      ::Error("ParseCIGAR", "malformed CIGAR '%s'", cigar_str.c_str());
      cigar_ops.clear();
      return cigar_ops;
   };

   for (char c : cigar_str) {
      if (c >= '0' && c <= '9') {
         length = length * 10 + static_cast<uint64_t>(c - '0');
         // Bounded here rather than by std::stoul, which threw std::out_of_range
         // on a long run of digits and took the whole conversion down with it.
         if (length > kMaxCigarOpLen)
            return fail();
         have_length = true;
         continue;
      }

      const uint8_t op = kCigarToCode[static_cast<uint8_t>(c)];
      // An unrecognised operator used to index a zero-filled table and come back
      // as 0 -- 'M' -- so "10Q" was silently stored as ten matches.
      if (op == kCigarCodeInvalid || !have_length)
         return fail();

      cigar_ops.push_back((static_cast<uint32_t>(length) << 4) | op);
      length = 0;
      have_length = false;
   }

   if (have_length)
      return fail();
   return cigar_ops;
}

std::string FormatCIGAR(const std::vector<uint32_t> &cigar_ops)
{
   // Without this an unmapped read comes back with an empty CIGAR column
   // producing a malformed SAM record
   if (cigar_ops.empty())
      return "*";

   std::string out{};
   for (uint32_t op : cigar_ops) {
      out += std::to_string(op >> 4);
      out += kCodeToCigar.at(op & 0xf);
   }
   return out;
}

} // namespace RAMNTupleUtils
