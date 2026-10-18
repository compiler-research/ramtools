#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static4x16.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

// After <cstdint>: the header uses uint8_t without including it.
#include <htscodecs/tokenise_name3.h>

namespace {

// CRAM 3.1's variant: reverse-strand qualities are modelled in sequencing order.
constexpr int kFqzVersion = 3;
// fqzcomp's parameter set for MiSeq-like data; the smallest on HG00154.
constexpr int kFqzStrategy = 2;
// The name tokeniser's level, of 1 to 9. On NA12878 chr20, level 9 with
// arithmetic coding saves 2 MB more of 61 MB, but encodes 4.7 times slower.
constexpr int kTok3Level = 3;
// rANS Nx16 order 1 with RLE and PACK, as CRAM 3.1 compresses its quality
// series. On NA12878 chr20: 72.4 MB, against 71.6 MB with fqzcomp, which
// decodes 25 times slower, and 114.6 MB for the qualities as a zstd column.
constexpr int kRansOrder = 1 | RANS_ORDER_RLE | RANS_ORDER_PACK;

// The first byte of a quality block says how the rest is compressed. A rANS
// block holds the size of its lengths part in four bytes, little endian, the
// lengths as varints, then the qualities, each part compressed with kRansOrder.
// An fqzcomp block is fqzcomp's own stream, which holds the lengths.
constexpr char kQualRans = 1;
constexpr char kQualFqzcomp = 2;

bool IsSamQuality(const std::string &qual)
{
   return std::all_of(qual.begin(), qual.end(), [](char c) { return c >= '!' && c <= '~'; });
}

// fqzcomp and the name tokeniser take char, rANS unsigned char.
unsigned char *AsBytes(char *p)
{
   return static_cast<unsigned char *>(static_cast<void *>(p));
}

std::string Rans(std::string &in)
{
   unsigned int size = 0;
   unsigned char *out = rans_compress_4x16(AsBytes(in.data()), static_cast<unsigned int>(in.size()), &size, kRansOrder);
   if (!out)
      throw std::runtime_error("rANS could not compress a quality block");
   std::string bytes(size, '\0');
   std::memcpy(bytes.data(), out, size);
   free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   return bytes;
}

std::string Unrans(char *in, std::size_t size)
{
   unsigned int n = 0;
   unsigned char *out = rans_uncompress_4x16(AsBytes(in), static_cast<unsigned int>(size), &n);
   if (!out)
      throw std::runtime_error("rANS could not decode a quality block");
   std::string bytes(n, '\0');
   std::memcpy(bytes.data(), out, n);
   free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   return bytes;
}

// The block of \p lengths and \p quals, as kQualRans describes it.
std::string RansBlock(const std::vector<uint32_t> &lengths, std::string &quals)
{
   std::string varints;
   for (uint32_t n : lengths) {
      for (; n >= 0x80; n >>= 7)
         varints.push_back(static_cast<char>((n & 0x7f) | 0x80));
      varints.push_back(static_cast<char>(n));
   }
   const std::string packedLengths = Rans(varints);
   std::string block(1, kQualRans);
   for (int shift = 0; shift < 32; shift += 8)
      block.push_back(static_cast<char>((packedLengths.size() >> shift) & 0xff));
   block += packedLengths;
   block += Rans(quals);
   return block;
}

// Fills \p lengths and \p quals from the block after its codec byte, for \p records qualities.
void DecodeRansBlock(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths, std::string &quals)
{
   constexpr std::size_t kHeader = 5;
   if (bytes.size() < kHeader)
      throw std::runtime_error("a quality block is damaged");
   std::size_t lengthsSize = 0;
   for (std::size_t i = 0; i < 4; i++)
      lengthsSize |= static_cast<std::size_t>(static_cast<unsigned char>(bytes.at(1 + i))) << (8 * i);
   if (lengthsSize >= bytes.size() - kHeader)
      throw std::runtime_error("a quality block is damaged");

   const std::string varints = Unrans(&bytes.at(kHeader), lengthsSize);
   lengths.clear();
   std::size_t total = 0;
   uint32_t n = 0;
   int shift = 0;
   for (const char c : varints) {
      const auto byte = static_cast<unsigned char>(c);
      if (shift > 28)
         throw std::runtime_error("a quality block is damaged");
      n |= static_cast<uint32_t>(byte & 0x7f) << shift;
      shift += 7;
      if (byte < 0x80) {
         if (n > static_cast<uint32_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("a quality block is damaged");
         lengths.push_back(static_cast<int>(n));
         total += n;
         n = 0;
         shift = 0;
      }
   }
   quals = Unrans(&bytes.at(kHeader + lengthsSize), bytes.size() - kHeader - lengthsSize);
   if (shift != 0 || lengths.size() != records || quals.size() != total)
      throw std::runtime_error("a quality block is damaged");
}

// Ends that do not rise, or that run past the last row, would send a row to a
// block that does not hold it.
void CheckBlockEnds(const std::vector<uint64_t> &ends, ROOT::NTupleSize_t rows)
{
   for (std::size_t i = 0; i < ends.size(); i++) {
      if ((i > 0 && ends[i] <= ends[i - 1]) || ends[i] >= rows)
         throw std::runtime_error("the block ends in the file are damaged");
   }
}

} // namespace

BlockWriter::BlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second, FillFn fill,
                         std::size_t block_records)
   : m_entries{std::move(first), std::move(second)},
     m_fill(std::move(fill)),
     m_blockRecords(std::max<std::size_t>(block_records, 1))
{
   for (std::size_t i = 0; i < 2; i++) {
      m_records.at(i) = m_entries.at(i)->GetPtr<RAMNTupleRecord>("record").get();
      m_qualBlobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField).get();
      m_nameBlobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField).get();
   }
}

void BlockWriter::Add()
{
   RAMNTupleRecord &rec = *m_records.at(m_current);
   m_names += rec.qname;
   m_names += '\0';

   // kFqzcomp comes in from the quality policy and picks the block's codec;
   // kQualInBlock goes out meaning the quality is in the block. fqzcomp takes
   // no empty strings, so "*" stays in the record. So do the lossy policies:
   // their bytes only read back through DecodeQuality, and binned qualities
   // from Q30 up look like SAM quality text.
   if (rec.TestBit(RAMNTupleRecord::kFqzcomp))
      m_fqzcomp = true;
   rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kFqzcomp);
   const bool lossless = !rec.TestBit(RAMNTupleRecord::kIlluminaBinning | RAMNTupleRecord::kDrop);
   if (lossless && rec.qual != "*" && !rec.qual.empty() && IsSamQuality(rec.qual)) {
      for (const char c : rec.qual)
         m_quals.push_back(static_cast<char>(c - '!'));
      m_lengths.push_back(static_cast<uint32_t>(rec.qual.size()));
      m_flags.push_back(((rec.flag & 0x10) ? FQZ_FREVERSE : 0) | ((rec.flag & 0x80) ? FQZ_FREAD2 : 0));
      rec.qual.clear();
      rec.SetBit(RAMNTupleRecord::kQualInBlock);
   } else {
      rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kQualInBlock);
   }
   m_rows++;

   if (m_pending)
      m_fill(*m_entries.at(1 - m_current));
   m_pending = true;
   m_current = 1 - m_current;
   if (m_rows == m_blockRecords)
      Finish();
}

void BlockWriter::Finish()
{
   if (!m_pending)
      return;
   const std::size_t last = 1 - m_current;
   if (!m_lengths.empty() && m_fqzcomp) {
      fqz_slice slice{};
      slice.num_records = static_cast<int>(m_lengths.size());
      slice.len = m_lengths.data();
      slice.flags = m_flags.data();
      std::size_t size = 0;
      char *packed =
         fqz_compress(kFqzVersion, &slice, m_quals.data(), m_quals.size(), &size, kFqzStrategy, /*gp=*/nullptr);
      if (!packed)
         throw std::runtime_error("fqzcomp could not compress a quality block");
      const std::string_view bytes(packed, size);
      m_qualBlobs.at(last)->assign(1, kQualFqzcomp);
      m_qualBlobs.at(last)->insert(m_qualBlobs.at(last)->end(), bytes.begin(), bytes.end());
      free(packed); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   } else if (!m_lengths.empty()) {
      const std::string bytes = RansBlock(m_lengths, m_quals);
      m_qualBlobs.at(last)->assign(bytes.begin(), bytes.end());
   }
   int size = 0;
   int lastStart = 0;
   uint8_t *names = tok3_encode_names(m_names.data(), static_cast<int>(m_names.size()), kTok3Level,
                                      /*use_arith=*/0, &size, &lastStart);
   if (!names)
      throw std::runtime_error("the name tokeniser could not compress a block of names");
   m_nameBlobs.at(last)->resize(static_cast<std::size_t>(size));
   std::memcpy(m_nameBlobs.at(last)->data(), names, static_cast<std::size_t>(size));
   free(names); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs

   m_fill(*m_entries.at(last));
   m_qualBlobs.at(last)->clear();
   m_nameBlobs.at(last)->clear();

   m_rowsSinceTake += m_rows;
   m_blockEnds.push_back(m_rowsSinceTake - 1);
   m_pending = false;
   m_rows = 0;
   m_names.clear();
   m_quals.clear();
   m_lengths.clear();
   m_flags.clear();
   m_fqzcomp = false;
}

std::vector<uint64_t> BlockWriter::TakeBlockEnds()
{
   m_rowsSinceTake = 0;
   return std::exchange(m_blockEnds, {});
}

QualityBlockReader::QualityBlockReader(RAMFile &file)
   : m_flagsView(file.GetReader().GetView<uint32_t>("record.compression_flags")),
     m_ends(file.GetMetadata().block_ends),
     m_readAhead(std::min<std::size_t>(kQualityReadAhead, std::max(1U, std::thread::hardware_concurrency())))
{
   ROOT::RNTupleReader &reader = file.GetReader();
   if (reader.GetDescriptor().FindFieldId(RAMNTupleRecord::kQualBlockField) != ROOT::kInvalidDescriptorId)
      m_view.emplace(reader.GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField));
   CheckBlockEnds(m_ends, reader.GetNEntries());
}

std::string QualityBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kQualInBlock))
      return rec.GetQUAL();

   const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
   if (it == m_ends.end() || !m_view)
      throw std::runtime_error("record " + std::to_string(row) + " has no quality block");
   const auto block = static_cast<std::size_t>(it - m_ends.begin());
   if (block != m_block)
      Load(block);

   const std::size_t slot = m_current.slots[static_cast<std::size_t>(row - m_firstRow)];
   return m_current.quals.substr(m_current.offsets[slot], static_cast<std::size_t>(m_current.lengths[slot]));
}

QualityBlockReader::Packed QualityBlockReader::Read(std::size_t block)
{
   const uint64_t first = block == 0 ? 0 : m_ends[block - 1] + 1;
   Packed p;
   p.slots.assign(static_cast<std::size_t>(m_ends[block] - first + 1), 0);
   for (std::size_t i = 0; i < p.slots.size(); i++) {
      if (m_flagsView(first + i) & RAMNTupleRecord::kQualInBlock)
         p.slots[i] = p.records++;
   }
   if (p.records > 0) {
      const std::vector<std::uint8_t> &bytes = (*m_view)(m_ends[block]);
      p.bytes.assign(bytes.begin(), bytes.end());
   }
   return p;
}

QualityBlockReader::Block QualityBlockReader::Decode(Packed packed)
{
   Block b;
   b.slots = std::move(packed.slots);
   b.lengths.assign(packed.records, 0);
   if (packed.records > 0) {
      std::vector<char> &bytes = packed.bytes;
      if (bytes.size() < 2)
         throw std::runtime_error("a quality block is damaged");
      if (bytes.front() == kQualRans) {
         DecodeRansBlock(bytes, packed.records, b.lengths, b.quals);
      } else if (bytes.front() == kQualFqzcomp) {
         std::size_t size = 0;
         char *out =
            fqz_decompress(&bytes.at(1), bytes.size() - 1, &size, b.lengths.data(), static_cast<int>(packed.records));
         if (!out)
            throw std::runtime_error("fqzcomp could not decode a quality block");
         b.quals.assign(out, size);
         free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
      } else {
         throw std::runtime_error("a quality block has an unknown codec");
      }
      for (auto &c : b.quals)
         c = static_cast<char>(c + '!');
   }
   b.offsets.resize(packed.records);
   std::size_t offset = 0;
   for (std::size_t i = 0; i < packed.records; i++) {
      b.offsets[i] = offset;
      offset += static_cast<std::size_t>(b.lengths[i]);
   }
   return b;
}

void QualityBlockReader::Load(std::size_t block)
{
   // Read ahead one block more for each consecutive block, so random lookups start nothing.
   m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;

   auto ahead = m_ahead.find(block);
   if (ahead != m_ahead.end()) {
      m_current = ahead->second.get();
      m_ahead.erase(ahead);
   } else {
      m_current = Decode(Read(block));
   }
   m_block = block;
   m_firstRow = block == 0 ? 0 : m_ends[block - 1] + 1;

   const std::size_t depth = std::min(m_run, m_readAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth; next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next)));
   }
}

NameBlockReader::NameBlockReader(RAMFile &file)
   : m_view(file.GetReader().GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField)),
     m_ends(file.GetMetadata().block_ends)
{
   CheckBlockEnds(m_ends, file.GetReader().GetNEntries());
}

std::string_view NameBlockReader::Get(ROOT::NTupleSize_t row)
{
   if (m_block == kNoBlock || row < m_firstRow || row > m_ends[m_block]) {
      const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
      if (it == m_ends.end())
         throw std::runtime_error("record " + std::to_string(row) + " has no block of names");
      Load(static_cast<std::size_t>(it - m_ends.begin()));
   }
   const auto i = static_cast<std::size_t>(row - m_firstRow);
   const std::vector<std::size_t> &starts = m_current.starts;
   return std::string_view(m_current.names).substr(starts[i], starts[i + 1] - starts[i] - 1);
}

std::vector<std::uint8_t> NameBlockReader::Read(std::size_t block)
{
   const std::vector<std::uint8_t> &bytes = m_view(m_ends[block]);
   return {bytes.begin(), bytes.end()};
}

NameBlockReader::Block NameBlockReader::Decode(std::vector<std::uint8_t> packed, std::size_t block, std::size_t rows)
{
   uint32_t size = 0;
   uint8_t *names = tok3_decode_names(packed.data(), static_cast<uint32_t>(packed.size()), &size);
   if (!names)
      throw std::runtime_error("the name tokeniser could not decode block " + std::to_string(block));
   Block b;
   b.names.resize(size);
   std::memcpy(b.names.data(), names, size);
   free(names); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs

   b.starts.reserve(rows + 1);
   b.starts.push_back(0);
   for (std::size_t end = b.names.find('\0'); end != std::string::npos; end = b.names.find('\0', end + 1))
      b.starts.push_back(end + 1);
   if (b.starts.size() != rows + 1)
      throw std::runtime_error("block " + std::to_string(block) + " has " + std::to_string(b.starts.size() - 1) +
                               " names for " + std::to_string(rows) + " records");
   return b;
}

void NameBlockReader::Load(std::size_t block)
{
   // Read ahead only while the blocks come in order, so random lookups start nothing.
   m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;
   const auto rows = [this](std::size_t b) {
      return static_cast<std::size_t>(m_ends[b] - (b == 0 ? 0 : m_ends[b - 1] + 1) + 1);
   };

   auto ahead = m_ahead.find(block);
   if (ahead != m_ahead.end()) {
      m_current = ahead->second.get();
      m_ahead.erase(ahead);
   } else {
      m_current = Decode(Read(block), block, rows(block));
   }
   m_block = block;
   m_firstRow = block == 0 ? 0 : m_ends[block - 1] + 1;

   const std::size_t depth = std::min(m_run, kNameReadAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth; next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next), next, rows(next)));
   }
}
