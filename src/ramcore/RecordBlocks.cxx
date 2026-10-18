#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <htscodecs/fqzcomp_qual.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <future>
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

bool IsSamQuality(const std::string &qual)
{
   return std::all_of(qual.begin(), qual.end(), [](char c) { return c >= '!' && c <= '~'; });
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

   // kQualInBlock comes in from the quality policy as a request (-fqzcomp) and
   // goes out meaning the quality is in the block. fqzcomp takes no empty
   // strings, so "*" stays in the record. So do the lossy policies: their
   // bytes only read back through DecodeQuality, and binned qualities from
   // Q30 up look like SAM quality text.
   const bool wanted = rec.TestBit(RAMNTupleRecord::kQualInBlock);
   const bool lossless = !rec.TestBit(RAMNTupleRecord::kIlluminaBinning | RAMNTupleRecord::kDrop);
   if (wanted && lossless && rec.qual != "*" && !rec.qual.empty() && IsSamQuality(rec.qual)) {
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
   if (!m_lengths.empty()) {
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
      m_qualBlobs.at(last)->assign(bytes.begin(), bytes.end());
      free(packed); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
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
      std::size_t size = 0;
      char *out = fqz_decompress(packed.bytes.data(), packed.bytes.size(), &size, b.lengths.data(),
                                 static_cast<int>(packed.records));
      if (!out)
         throw std::runtime_error("fqzcomp could not decode a quality block");
      b.quals.assign(out, size);
      free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
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
