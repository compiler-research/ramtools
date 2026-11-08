#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <htscodecs/fqzcomp_qual.h>
#include <htscodecs/rANS_static4x16.h>
#include <libdeflate.h>

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

// DEFLATE, which CRAM 3.1 uses for the bases of a file without a reference, at
// level 6. On NA12878 chr20 the bases take 84.6 MB, against 96.1 MB as a column
// with ROOT's 505 (zstd at level 10, since ROOT gives zstd twice its level), in
// half the CPU, and decode as fast. The gain is the compressor's, not the
// blocks': zstd at level 10 in blocks gives 95.2 MB. Level 9 gives 81.9 MB in
// six times the CPU.
constexpr int kDeflateLevel = 6;

// The first byte of a block of qualities or bases says how the rest is
// compressed. A rANS or DEFLATE block holds the size of its lengths part in
// four bytes, little endian, the lengths as varints compressed with kRansOrder,
// then the values, compressed with rANS (kRansOrder) or with raw DEFLATE. A
// quality block has the lengths of the records whose quality it holds, a block
// of bases one for every record of the block. An fqzcomp block is fqzcomp's own
// stream, which holds the lengths.
constexpr char kQualRans = 1;
constexpr char kQualFqzcomp = 2;
constexpr char kSeqDeflate = 1;

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
      throw std::runtime_error("rANS could not compress a block");
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
      throw std::runtime_error("rANS could not decode a block");
   std::string bytes(n, '\0');
   std::memcpy(bytes.data(), out, n);
   free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   return bytes;
}

// A block with the codec byte \p codec, \p lengths and the compressed values
// \p packed, as kQualRans and kSeqDeflate describe it.
std::string BlockWithLengths(char codec, const std::vector<uint32_t> &lengths, const std::string &packed)
{
   std::string varints;
   for (uint32_t n : lengths) {
      for (; n >= 0x80; n >>= 7)
         varints.push_back(static_cast<char>((n & 0x7f) | 0x80));
      varints.push_back(static_cast<char>(n));
   }
   const std::string packedLengths = Rans(varints);
   std::string block(1, codec);
   for (int shift = 0; shift < 32; shift += 8)
      block.push_back(static_cast<char>((packedLengths.size() >> shift) & 0xff));
   block += packedLengths;
   block += packed;
   return block;
}

// Fills \p lengths from the lengths part of a block that holds \p records
// values, and returns where its values start and, in \p total, their size.
// \p what names the block in errors.
std::size_t ReadLengths(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths, std::size_t &total,
                        const char *what)
{
   const std::string damaged = std::string(what) + " is damaged";
   constexpr std::size_t kHeader = 5;
   if (bytes.size() < kHeader)
      throw std::runtime_error(damaged);
   std::size_t lengthsSize = 0;
   for (std::size_t i = 0; i < 4; i++)
      lengthsSize |= static_cast<std::size_t>(static_cast<unsigned char>(bytes.at(1 + i))) << (8 * i);
   if (lengthsSize >= bytes.size() - kHeader)
      throw std::runtime_error(damaged);

   const std::string varints = Unrans(&bytes.at(kHeader), lengthsSize);
   lengths.clear();
   total = 0;
   uint32_t n = 0;
   int shift = 0;
   for (const char c : varints) {
      const auto byte = static_cast<unsigned char>(c);
      if (shift > 28)
         throw std::runtime_error(damaged);
      n |= static_cast<uint32_t>(byte & 0x7f) << shift;
      shift += 7;
      if (byte < 0x80) {
         if (n > static_cast<uint32_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error(damaged);
         lengths.push_back(static_cast<int>(n));
         total += n;
         n = 0;
         shift = 0;
      }
   }
   if (shift != 0 || lengths.size() != records)
      throw std::runtime_error(damaged);
   return kHeader + lengthsSize;
}

// Fills \p lengths and \p quals, as SAM text, from a quality block that holds
// \p records qualities.
void DecodeQualityBlock(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths, std::string &quals)
{
   if (bytes.size() < 2)
      throw std::runtime_error("a quality block is damaged");
   if (bytes.front() == kQualRans) {
      std::size_t total = 0;
      const std::size_t start = ReadLengths(bytes, records, lengths, total, /*what=*/"a quality block");
      // Straight into quals, which may hold an earlier block: no copy, and only
      // a larger block's tail is zero-filled.
      quals.resize(total);
      auto size = static_cast<unsigned int>(total);
      if (total > std::numeric_limits<unsigned int>::max() ||
          !rans_uncompress_to_4x16(AsBytes(&bytes.at(start)), static_cast<unsigned int>(bytes.size() - start),
                                   AsBytes(quals.data()), &size) ||
          size != total)
         throw std::runtime_error("a quality block is damaged");
   } else if (bytes.front() == kQualFqzcomp) {
      std::size_t size = 0;
      lengths.assign(records, 0);
      char *out = fqz_decompress(&bytes.at(1), bytes.size() - 1, &size, lengths.data(), static_cast<int>(records));
      if (!out)
         throw std::runtime_error("fqzcomp could not decode a quality block");
      quals.assign(out, size);
      free(out); // NOLINT(cppcoreguidelines-owning-memory,cppcoreguidelines-no-malloc): allocated by htscodecs
   } else {
      throw std::runtime_error("a quality block has an unknown codec");
   }
   for (auto &c : quals)
      c = static_cast<char>(c + '!');
}

// A quality block of fqzcomp, which decodes 25 times slower than rANS.
bool IsFqzcompBlock(const std::vector<char> &bytes)
{
   return !bytes.empty() && bytes.front() == kQualFqzcomp;
}

// Fills \p lengths and \p seqs from a block that holds \p records sequences.
void DecodeSeqBlock(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths, std::string &seqs)
{
   if (bytes.size() < 2)
      throw std::runtime_error("a block of bases is damaged");
   if (bytes.front() != kSeqDeflate)
      throw std::runtime_error("a block of bases has an unknown codec");
   std::size_t total = 0;
   const std::size_t start = ReadLengths(bytes, records, lengths, total, /*what=*/"a block of bases");
   // seqs may hold an earlier block, and inflating writes every byte, so only a
   // larger block's tail is zero-filled.
   seqs.resize(total);
   const std::unique_ptr<libdeflate_decompressor, void (*)(libdeflate_decompressor *)> inflate(
      libdeflate_alloc_decompressor(), libdeflate_free_decompressor);
   if (!inflate)
      throw std::bad_alloc();
   // Without a place for the actual size, libdeflate fails unless the stream
   // gives exactly `total` bytes.
   if (libdeflate_deflate_decompress(inflate.get(), &bytes.at(start), bytes.size() - start, seqs.data(), total,
                                     nullptr) != LIBDEFLATE_SUCCESS)
      throw std::runtime_error("a block of bases is damaged");
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

// For each block, a number for the cluster that holds its last row, and so its
// bytes: equal numbers, the same cluster.
std::vector<std::size_t> BlockClusters(const ROOT::RNTupleDescriptor &desc, const std::vector<uint64_t> &ends)
{
   std::vector<uint64_t> firsts;
   for (const auto &cluster : desc.GetClusterIterable())
      firsts.push_back(cluster.GetFirstEntryIndex());
   std::sort(firsts.begin(), firsts.end());
   std::vector<std::size_t> clusters;
   clusters.reserve(ends.size());
   for (const uint64_t end : ends)
      clusters.push_back(static_cast<std::size_t>(std::upper_bound(firsts.begin(), firsts.end(), end) - firsts.begin()));
   return clusters;
}

} // namespace

void BlockWriter::FreeDeflate::operator()(libdeflate_compressor *c) const
{
   libdeflate_free_compressor(c);
}

BlockWriter::BlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second, FillFn fill,
                         std::size_t block_records)
   : m_entries{std::move(first), std::move(second)},
     m_fill(std::move(fill)),
     m_blockRecords(std::max<std::size_t>(block_records, 1)),
     m_deflate(libdeflate_alloc_compressor(kDeflateLevel))
{
   if (!m_deflate)
      throw std::bad_alloc();
   for (std::size_t i = 0; i < 2; i++) {
      m_records.at(i) = m_entries.at(i)->GetPtr<RAMNTupleRecord>("record").get();
      m_qualBlobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField).get();
      m_nameBlobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField).get();
      m_seqBlobs.at(i) = m_entries.at(i)->GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kSeqBlockField).get();
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

   // kSeqInBlock goes out meaning the bases are in the block; "*" stays in the
   // record. The block has a length for every record, 0 for those, so that a
   // reader finds a record's bases without reading the flags.
   if (!rec.seq.empty() && rec.seq != "*") {
      m_seqs += rec.seq;
      m_seqLengths.push_back(static_cast<uint32_t>(rec.seq.size()));
      rec.seq.clear();
      rec.SetBit(RAMNTupleRecord::kSeqInBlock);
   } else {
      m_seqLengths.push_back(0);
      rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kSeqInBlock);
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
      const std::string bytes = BlockWithLengths(kQualRans, m_lengths, Rans(m_quals));
      m_qualBlobs.at(last)->assign(bytes.begin(), bytes.end());
   }
   if (!m_seqs.empty()) {
      std::string packed(libdeflate_deflate_compress_bound(m_deflate.get(), m_seqs.size()), '\0');
      const std::size_t size =
         libdeflate_deflate_compress(m_deflate.get(), m_seqs.data(), m_seqs.size(), packed.data(), packed.size());
      if (size == 0)
         throw std::runtime_error("libdeflate could not compress a block of bases");
      packed.resize(size);
      const std::string bytes = BlockWithLengths(kSeqDeflate, m_seqLengths, packed);
      m_seqBlobs.at(last)->assign(bytes.begin(), bytes.end());
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
   m_seqBlobs.at(last)->clear();

   m_rowsSinceTake += m_rows;
   m_blockEnds.push_back(m_rowsSinceTake - 1);
   m_pending = false;
   m_rows = 0;
   m_names.clear();
   m_quals.clear();
   m_lengths.clear();
   m_flags.clear();
   m_fqzcomp = false;
   m_seqs.clear();
   m_seqLengths.clear();
}

std::vector<uint64_t> BlockWriter::TakeBlockEnds()
{
   m_rowsSinceTake = 0;
   return std::exchange(m_blockEnds, {});
}

FlaggedBlockReader::FlaggedBlockReader(RAMFile &file, const char *field, uint32_t bit, DecodeFn decode, SlowFn slow,
                                       std::size_t readAhead, bool everyRow)
   : m_field(field),
     m_bit(bit),
     m_decode(decode),
     m_slow(slow),
     m_everyRow(everyRow),
     m_flagsView(file.GetReader().GetView<uint32_t>("record.compression_flags")),
     m_ends(file.GetMetadata().block_ends),
     m_readAhead(std::min<std::size_t>(readAhead, std::max(1U, std::thread::hardware_concurrency())))
{
   ROOT::RNTupleReader &reader = file.GetReader();
   if (reader.GetDescriptor().FindFieldId(field) != ROOT::kInvalidDescriptorId)
      m_view.emplace(reader.GetView<std::vector<std::uint8_t>>(field));
   CheckBlockEnds(m_ends, reader.GetNEntries());
   if (file.UsesClusterCache())
      m_clusters = BlockClusters(reader.GetDescriptor(), m_ends);
}

std::string_view FlaggedBlockReader::FromBlock(ROOT::NTupleSize_t row)
{
   if (m_block == kNoBlock || row < m_firstRow || row > m_ends[m_block]) {
      const auto it = std::lower_bound(m_ends.begin(), m_ends.end(), static_cast<uint64_t>(row));
      if (it == m_ends.end() || !m_view)
         throw std::runtime_error("record " + std::to_string(row) + " has no block in " + m_field);
      Load(static_cast<std::size_t>(it - m_ends.begin()));
   }

   const auto i = static_cast<std::size_t>(row - m_firstRow);
   const std::size_t slot = m_everyRow ? i : m_current.slots[i];
   const std::string_view values = m_current.values;
   return values.substr(m_current.offsets[slot], static_cast<std::size_t>(m_current.lengths[slot]));
}

FlaggedBlockReader::Block FlaggedBlockReader::Read(std::size_t block)
{
   Block b;
   if (!m_spare.empty()) {
      b = std::move(m_spare.back());
      m_spare.pop_back();
   }
   const uint64_t first = block == 0 ? 0 : m_ends[block - 1] + 1;
   const auto rows = static_cast<std::size_t>(m_ends[block] - first + 1);
   b.records = 0;
   b.bytes.clear();
   b.slots.clear();
   if (m_everyRow) {
      // A block without values is empty, and no row of it has the bit.
      const std::vector<std::uint8_t> &bytes = (*m_view)(m_ends[block]);
      if (!bytes.empty()) {
         b.records = rows;
         b.bytes.assign(bytes.begin(), bytes.end());
      }
      return b;
   }
   b.slots.assign(rows, 0);
   for (std::size_t i = 0; i < rows; i++) {
      if (m_flagsView(first + i) & m_bit)
         b.slots[i] = b.records++;
   }
   if (b.records > 0) {
      const std::vector<std::uint8_t> &bytes = (*m_view)(m_ends[block]);
      b.bytes.assign(bytes.begin(), bytes.end());
   }
   return b;
}

FlaggedBlockReader::Block FlaggedBlockReader::Decode(Block b, DecodeFn decode)
{
   if (b.records > 0) {
      decode(b.bytes, b.records, b.lengths, b.values);
   } else {
      b.lengths.clear();
      b.values.clear();
   }
   b.offsets.resize(b.records);
   std::size_t offset = 0;
   for (std::size_t i = 0; i < b.records; i++) {
      b.offsets[i] = offset;
      offset += static_cast<std::size_t>(b.lengths[i]);
   }
   return b;
}

void FlaggedBlockReader::Load(std::size_t block)
{
   // Read ahead one block more for each consecutive block, so random lookups start nothing.
   m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;

   Block next;
   auto ahead = m_ahead.find(block);
   if (ahead != m_ahead.end()) {
      // Out of the map first, so that a block that failed fails again if asked for again.
      std::future<Block> done = std::move(ahead->second);
      m_ahead.erase(ahead);
      next = done.get();
   } else {
      next = Decode(Read(block), m_decode);
   }
   // One spare is enough during a scan: each block loaded frees one, and each
   // block read ahead takes one.
   if (m_spare.size() <= m_readAhead)
      m_spare.push_back(std::move(m_current));
   m_current = std::move(next);
   m_block = block;
   m_firstRow = block == 0 ? 0 : m_ends[block - 1] + 1;

   // With ROOT's cluster cache, only blocks of this block's cluster: the cache
   // keeps the cluster it is asked for and the next one, so a block read ahead
   // from the next cluster dropped the one the records were still in, and they
   // read it again. A full read of NA12878's HiSeq chromosome 20 read 6.6 GB of
   // its 2.2 GB file. Without the cache pages are read one by one, and waiting
   // for the first block of each cluster only made a region slower. So did it
   // for fqzcomp, whose blocks take about 60 ms each: 5.5 to 5.9 s instead of
   // 4.4 to 4.5 s for 2 million records, so a slow block still reads ahead.
   const bool anyCluster = m_clusters.empty() || (m_slow != nullptr && m_slow(m_current.bytes));
   const std::size_t depth = std::min(m_run, m_readAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth &&
                                      (anyCluster || m_clusters[next] == m_clusters[block]);
        next++) {
      // A thread for each block. Threads kept for the reader's life took more
      // CPU: 28.4 s instead of 25.7 s for a full read of chr20.
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next), m_decode));
   }
}

QualityBlockReader::QualityBlockReader(RAMFile &file)
   : FlaggedBlockReader(file, RAMNTupleRecord::kQualBlockField, RAMNTupleRecord::kQualInBlock, DecodeQualityBlock,
                        IsFqzcompBlock, kQualityReadAhead, /*everyRow=*/false)
{
}

std::string_view QualityBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kQualInBlock)) {
      m_text = rec.GetQUAL();
      return m_text;
   }
   return FromBlock(row);
}

SeqBlockReader::SeqBlockReader(RAMFile &file)
   : FlaggedBlockReader(file, RAMNTupleRecord::kSeqBlockField, RAMNTupleRecord::kSeqInBlock, DecodeSeqBlock,
                        /*slow=*/nullptr, kSeqReadAhead, /*everyRow=*/true)
{
}

std::string_view SeqBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kSeqInBlock))
      return rec.GetSEQ();
   return FromBlock(row);
}

NameBlockReader::NameBlockReader(RAMFile &file)
   : m_view(file.GetReader().GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kNameBlockField)),
     m_ends(file.GetMetadata().block_ends)
{
   CheckBlockEnds(m_ends, file.GetReader().GetNEntries());
   if (file.UsesClusterCache())
      m_clusters = BlockClusters(file.GetReader().GetDescriptor(), m_ends);
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

   // With ROOT's cluster cache, only blocks of this block's cluster, as
   // FlaggedBlockReader::Load() says.
   const std::size_t depth = std::min(m_run, kNameReadAhead);
   for (std::size_t next = block + 1; next < m_ends.size() && next <= block + depth &&
                                      (m_clusters.empty() || m_clusters[next] == m_clusters[block]);
        next++) {
      if (m_ahead.count(next) == 0)
         m_ahead.emplace(next, std::async(std::launch::async, Decode, Read(next), next, rows(next)));
   }
}
