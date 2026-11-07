#ifndef RAMCORE_RECORDBLOCKS_H
#define RAMCORE_RECORDBLOCKS_H

// Records in blocks. The last record of a block carries the block's names,
// compressed with htscodecs' name tokeniser (tok3, CRAM 3.1's codec), in
// "nameblock"; its lossless qualities in "qualblock", compressed with
// htscodecs' rANS Nx16 as CRAM 3.1 compresses its quality series, or with
// fqzcomp when asked for; and its bases in "seqblock", compressed with
// libdeflate. METADATA lists the row each block ends on.

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

struct libdeflate_compressor;

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
/// Blocks of bases decoded ahead during a scan. A block of 10,000 reads of 150
/// bases decodes in about 0.6 ms, and a dump writes its records in about 8, so
/// one is enough: two took as long on chr20, with 1 to 2 MB more.
constexpr std::size_t kSeqReadAhead = 1;

/// Holds one record back so the record that ends a block can carry it. Every
/// record's name goes into its block. A quality goes into a block when it is
/// lossless SAM text in '!'..'~'; the rest stay in the record. The block is
/// compressed with fqzcomp when its records' policy has kFqzcomp, and with rANS
/// otherwise. Bases go into a block unless they are "*" or empty.
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
   std::array<std::vector<std::uint8_t> *, 2> m_seqBlobs{};
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
   std::string m_seqs;
   std::vector<uint32_t> m_seqLengths;
   struct FreeDeflate {
      void operator()(libdeflate_compressor *c) const;
   };
   std::unique_ptr<libdeflate_compressor, FreeDeflate> m_deflate; ///< Kept from block to block.
   std::vector<uint64_t> m_blockEnds;
};

/// Values that the records with a bit in compression_flags keep in the blocks
/// of a field: the qualities with kQualInBlock in qualblock, the bases with
/// kSeqInBlock in seqblock. Decodes one block at a time, and ahead on other
/// threads during a scan.
class FlaggedBlockReader {
protected:
   /// Fills \p lengths and \p values from the bytes of a block that holds
   /// \p records values. Runs on other threads.
   using DecodeFn = void (*)(std::vector<char> &bytes, std::size_t records, std::vector<int> &lengths,
                             std::string &values);

   /// With \p everyRow, a block holds a value for every row, empty for those
   /// without the bit, so finding a row's value needs no flags.
   FlaggedBlockReader(RAMFile &file, const char *field, uint32_t bit, DecodeFn decode, std::size_t readAhead,
                      bool everyRow);

   /// The value of \p row, whose record has the bit. Valid until the next call.
   /// A view into the decoded block, not a copy: copying every quality out cost
   /// 7% of a full dump's CPU.
   std::string_view FromBlock(ROOT::NTupleSize_t row);

private:
   struct Block {
      std::vector<char> bytes;        ///< The block as read.
      std::size_t records = 0;        ///< Values it holds.
      std::vector<std::size_t> slots; ///< Per row: index of its value; empty with everyRow.
      std::string values;
      std::vector<int> lengths;
      std::vector<std::size_t> offsets;
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Views are not thread-safe, so only Decode() runs on other threads. Read()
   /// fills a block that is done with, if there is one, so that decoding writes
   /// into memory the process already has: with new memory for every block,
   /// glibc gave each freed block back to the kernel, and reading only the bases
   /// of chr20 took 306,454 page faults and 1.1 s of system time.
   Block Read(std::size_t block);
   static Block Decode(Block b, DecodeFn decode);
   void Load(std::size_t block);

   std::string m_field;
   uint32_t m_bit;
   DecodeFn m_decode;
   bool m_everyRow;
   ROOT::RNTupleView<uint32_t> m_flagsView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> m_view;
   std::vector<uint64_t> m_ends; ///< Last row of each block.
   std::size_t m_readAhead;
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
   std::vector<Block> m_spare; ///< Blocks done with, for Read() to fill again.
};

/// Returns QUAL as SAM text, from its block or from the record.
class QualityBlockReader : public FlaggedBlockReader {
public:
   explicit QualityBlockReader(RAMFile &file);

   /// Valid until the next call.
   std::string_view Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   std::string m_text; ///< The last quality Get() took from its record, not from a block.
};

/// Returns SEQ, from its block or from the record.
class SeqBlockReader : public FlaggedBlockReader {
public:
   explicit SeqBlockReader(RAMFile &file);

   /// Valid until the next call, and while \p rec is.
   std::string_view Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);
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
