#include "ramcore/QualityBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <htscodecs/fqzcomp_qual.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace {

// CRAM 3.1's variant: reverse-strand qualities are modelled in sequencing order.
constexpr int kFqzVersion = 3;
// fqzcomp's parameter set for MiSeq-like data; the smallest on HG00154.
constexpr int kFqzStrategy = 2;

bool IsSamQuality(const std::string &qual)
{
   return std::all_of(qual.begin(), qual.end(), [](char c) { return c >= '!' && c <= '~'; });
}

} // namespace

QualityBlockWriter::QualityBlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second,
                                       FillFn fill, std::size_t block_records)
   : m_entries{std::move(first), std::move(second)},
     m_fill(std::move(fill)),
     m_blockRecords(std::max<std::size_t>(block_records, 1))
{
   for (std::size_t i = 0; i < 2; i++) {
      m_records.at(i) = m_entries.at(i)->GetPtr<RAMNTupleRecord>("record").get();
      m_blobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField).get();
   }
}

void QualityBlockWriter::Add()
{
   RAMNTupleRecord &rec = *m_records.at(m_current);
   // fqzcomp takes no empty strings, so "*" stays in the record.
   if (rec.TestBit(RAMNTupleRecord::kPhred33) && rec.qual != "*" && !rec.qual.empty() && IsSamQuality(rec.qual)) {
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

void QualityBlockWriter::Finish()
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
      m_blobs.at(last)->assign(bytes.begin(), bytes.end());
      free(packed); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   }
   m_fill(*m_entries.at(last));
   m_blobs.at(last)->clear();

   m_rowsSinceTake += m_rows;
   m_blockEnds.push_back(m_rowsSinceTake - 1);
   m_pending = false;
   m_rows = 0;
   m_quals.clear();
   m_lengths.clear();
   m_flags.clear();
}

std::vector<uint64_t> QualityBlockWriter::TakeBlockEnds()
{
   m_rowsSinceTake = 0;
   return std::exchange(m_blockEnds, {});
}

QualityBlockReader::QualityBlockReader(ROOT::RNTupleReader &reader)
   : m_flagsView(reader.GetView<uint32_t>("record.compression_flags")),
     m_readAhead(std::max(1U, std::thread::hardware_concurrency()))
{
   if (reader.GetDescriptor().FindFieldId(RAMNTupleRecord::kQualBlockField) != ROOT::kInvalidDescriptorId)
      m_view.emplace(reader.GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField));
}

std::string QualityBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kQualInBlock))
      return rec.GetQUAL();

   const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
   const auto it = std::lower_bound(ends.begin(), ends.end(), static_cast<uint64_t>(row));
   if (it == ends.end() || !m_view)
      throw std::runtime_error("record " + std::to_string(row) + " has no quality block");
   const auto block = static_cast<std::size_t>(it - ends.begin());
   if (block != m_block)
      Load(block);

   const std::size_t slot = m_current.slots[static_cast<std::size_t>(row - m_firstRow)];
   return m_current.quals.substr(m_current.offsets[slot], static_cast<std::size_t>(m_current.lengths[slot]));
}

QualityBlockReader::Packed QualityBlockReader::Read(std::size_t block)
{
   const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
   const uint64_t first = block == 0 ? 0 : ends[block - 1] + 1;
   Packed p;
   p.slots.assign(static_cast<std::size_t>(ends[block] - first + 1), 0);
   for (std::size_t i = 0; i < p.slots.size(); i++) {
      if (m_flagsView(first + i) & RAMNTupleRecord::kQualInBlock)
         p.slots[i] = p.records++;
   }
   if (p.records > 0) {
      const std::vector<std::uint8_t> &bytes = (*m_view)(ends[block]);
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
   const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
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
   m_firstRow = block == 0 ? 0 : ends[block - 1] + 1;

   const std::size_t depth = std::min(m_run, m_readAhead);
   for (std::size_t next = block + 1; next < ends.size() && next <= block + depth; next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next)));
   }
}
