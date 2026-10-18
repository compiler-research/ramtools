#ifndef RAMCORE_RECORDBLOCKS_H
#define RAMCORE_RECORDBLOCKS_H

// Records in blocks. The last record of a block carries the block's names,
// compressed with htscodecs' name tokeniser (tok3, CRAM 3.1's codec), in
// "nameblock", and its lossless qualities in "qualblock", compressed with
// htscodecs' rANS Nx16 as CRAM 3.1 compresses its quality series, or with
// fqzcomp when asked for. METADATA lists the row each block ends on.

#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

/// Records per block, as many as a CRAM slice. A region query decodes every
/// block it touches: on NA12878 chr20 about 3 ms for the names, 2 ms for rANS
/// qualities and 60 ms for fqzcomp ones. For names 10,000 gives 60.9 MB,
/// against 62.9 MB with 5,000 and 58.9 MB with 100,000; for fqzcomp it keeps
/// all but a percent of 25,000's gain at less than half the decode time.
constexpr std::size_t kBlockRecords = 10000;
/// Blocks decoded ahead during a scan. Four keep up with a full dump of fqzcomp
/// blocks, and each fqzcomp decoder holds a 68 MB model, so more only costs
/// memory.
constexpr std::size_t kQualityReadAhead = 4;
/// Blocks of names decoded ahead during a scan. A block decodes in about 3 ms,
/// and a dump writes its records in about 15, so two are plenty.
constexpr std::size_t kNameReadAhead = 2;

/// Holds one record back so the record that ends a block can carry it. Every
/// record's name goes into its block. A quality goes into a block when it is
/// lossless SAM text in '!'..'~'; the rest stay in the record. The block is
/// compressed with fqzcomp when its records' policy has kFqzcomp, and with rANS
/// otherwise.
class BlockWriter {
public:
   using FillFn = std::function<void(ROOT::REntry &)>;

   BlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second, FillFn fill,
               std::size_t block_records = kBlockRecords);

   RAMNTupleRecord &Record() { return *m_records.at(m_current); }
   void Add();
   /// Ends the current block; needed before the written rows are flushed.
   void Finish();
   /// Last row of each block finished since the previous call, counted from the
   /// first row after it.
   std::vector<uint64_t> TakeBlockEnds();

private:
   std::array<std::unique_ptr<ROOT::REntry>, 2> m_entries;
   std::array<RAMNTupleRecord *, 2> m_records{};
   std::array<std::vector<std::uint8_t> *, 2> m_qualBlobs{};
   std::array<std::vector<std::uint8_t> *, 2> m_nameBlobs{};
   FillFn m_fill;
   std::size_t m_blockRecords;
   std::size_t m_current = 0;
   bool m_pending = false;
   std::size_t m_rows = 0;
   uint64_t m_rowsSinceTake = 0;
   std::string m_names; ///< The block's names, each ended by '\0'.
   std::string m_quals;
   std::vector<uint32_t> m_lengths;
   std::vector<uint32_t> m_flags;
   bool m_fqzcomp = false; ///< A record of the block asked for fqzcomp.
   std::vector<uint64_t> m_blockEnds;
};

/// Returns QUAL as SAM text, decoding one block at a time and decoding ahead on
/// other threads during a scan.
class QualityBlockReader {
public:
   explicit QualityBlockReader(RAMFile &file);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Packed {
      std::vector<std::size_t> slots; ///< Per row: index of its quality in the block.
      std::size_t records = 0;
      std::vector<char> bytes;
   };
   struct Block {
      std::vector<std::size_t> slots;
      std::string quals;
      std::vector<int> lengths;
      std::vector<std::size_t> offsets;
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Views are not thread-safe, so only Decode() runs on other threads.
   Packed Read(std::size_t block);
   static Block Decode(Packed packed);
   void Load(std::size_t block);

   ROOT::RNTupleView<uint32_t> m_flagsView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> m_view;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_readAhead;
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
};

/// Returns QNAME, decoding one block at a time and decoding ahead on another
/// thread during a scan. Decoding takes about 0.3 us of CPU a name, most of it
/// in the name tokeniser's allocations.
class NameBlockReader {
public:
   explicit NameBlockReader(RAMFile &file);

   /// Valid until the next call.
   std::string_view Get(ROOT::NTupleSize_t row);

private:
   struct Block {
      std::string names;               ///< The block's names, each ended by '\0'.
      std::vector<std::size_t> starts; ///< Where each name starts, and one past the last.
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Views are not thread-safe, so only Decode() runs on other threads.
   std::vector<std::uint8_t> Read(std::size_t block);
   static Block Decode(std::vector<std::uint8_t> packed, std::size_t block, std::size_t rows);
   void Load(std::size_t block);

   ROOT::RNTupleView<std::vector<std::uint8_t>> m_view;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
};

#endif // RAMCORE_RECORDBLOCKS_H
