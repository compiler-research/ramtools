#include "ramcore/RecordBlocks.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct Read {
   uint16_t flag;
   std::string qual;
   std::string name = "r";
};

const char *const kFile = "record_blocks_test.root";
// The policies samtoramntuple passes for lossless qualities, by default and
// with -fqzcomp.
const uint32_t kRansPolicy = RAMNTupleRecord::kPhred33;
const uint32_t kFqzcompPolicy = RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kFqzcomp;
// The first byte of a quality block, by codec.
const uint8_t kRansCodec = 1;
const uint8_t kFqzcompCodec = 2;

using Damage = std::function<void(std::vector<std::uint8_t> &)>;

// Writes the reads; \p ends, when given, replaces the block ends, and \p damage
// changes every quality block, as in a damaged file.
void Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy, const std::vector<uint64_t> &ends = {},
           const Damage &damage = nullptr)
{
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   RAMFileMetadata meta;
   meta.rname_refs = {"chr1"};
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      BlockWriter out(
         writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
         [&writer, &damage](ROOT::REntry &e) {
            auto &blob = *e.GetPtr<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);
            if (damage && !blob.empty())
               damage(blob);
            writer->Fill(e);
         },
         block);
      for (const auto &read : reads) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetBit(policy);
         rec.SetQNAME(read.name);
         rec.SetFLAG(read.flag);
         rec.refid = 0;
         rec.SetQUAL(read.qual);
         out.Add();
      }
      out.Finish();
      meta.block_ends = ends.empty() ? out.TakeBlockEnds() : ends;
   }
   RAMNTupleRecord::WriteMetadata(*file, meta);
   file->Close();
}

void ExpectQualities(const std::vector<Read> &reads, const std::vector<std::size_t> &order)
{
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   for (const std::size_t row : order)
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
}

// The first byte of the quality block that row \p end carries, or 0 for none.
uint8_t CodecAt(uint64_t end)
{
   auto file = RAMFile::Open(kFile);
   EXPECT_NE(file, nullptr);
   if (!file)
      return 0;
   auto blocks = file->GetReader().GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kQualBlockField);
   const std::vector<std::uint8_t> &blob = blocks(end);
   return blob.empty() ? 0 : blob.front();
}

class RecordBlocksTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(RecordBlocksTest, EveryQualityComesBackInAnyOrder)
{
   const std::vector<Read> reads = {
      {0, "IIIIHHHGGF"},
      {0x10, "#####ABCDE"},
      {0x80, "!~!~!~"},
      {0, "*"}, // reuses row 1's record object, so its in-block bit must be cleared
      {0x90, "ABCDEFGHIJKLMNO"},
      {0, "AB C"}, // not SAM quality
      {0, "FFFFFFFFFF"},
      {0x10, "5"},
      {0, "*"},
      {0, "JJJJJJJJJJJJ"},
   };
   for (const uint32_t policy : {kRansPolicy, kFqzcompPolicy}) {
      Write(reads, /*block=*/3, policy);
      const auto written = RAMFile::Open(kFile);
      ASSERT_NE(written, nullptr);
      EXPECT_EQ(written->GetMetadata().block_ends, (std::vector<uint64_t>{2, 5, 8, 9}));

      std::vector<std::size_t> forward;
      std::vector<std::size_t> backward;
      for (std::size_t i = 0; i < reads.size(); i++) {
         forward.push_back(i);
         backward.push_back(reads.size() - 1 - i);
      }
      ExpectQualities(reads, forward);
      ExpectQualities(reads, backward);

      auto file = RAMFile::Open(kFile);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      for (std::size_t row = 0; row < reads.size(); row++) {
         const bool inBlock = reads[row].qual != "*" && reads[row].qual != "AB C";
         EXPECT_EQ(view(row).TestBit(RAMNTupleRecord::kQualInBlock), inBlock) << "row " << row;
         EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kFqzcomp)) << "a request, not kept in the record";
         if (inBlock)
            EXPECT_TRUE(view(row).qual.empty()) << "row " << row;
      }
   }
}

// Lossless qualities go into blocks by default, compressed with rANS; with
// -fqzcomp the same blocks are fqzcomp's.
TEST_F(RecordBlocksTest, EachBlockSaysItsCodec)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0x10, "#+5?"}, {0, "FFFF"}};
   Write(reads, /*block=*/2, kRansPolicy);
   EXPECT_EQ(CodecAt(1), kRansCodec);
   EXPECT_EQ(CodecAt(2), kRansCodec);
   Write(reads, /*block=*/2, kFqzcompPolicy);
   EXPECT_EQ(CodecAt(1), kFqzcompCodec);
   EXPECT_EQ(CodecAt(2), kFqzcompCodec);
}

TEST_F(RecordBlocksTest, ABlockWithoutQualitiesNeedsNoData)
{
   const std::vector<Read> reads = {{0, "*"}, {0, "*"}, {0, "*"}, {0, "HHHH"}};
   for (const uint32_t policy : {kRansPolicy, kFqzcompPolicy}) {
      Write(reads, /*block=*/3, policy);
      EXPECT_EQ(CodecAt(2), 0) << "no qualities, no block";
      ExpectQualities(reads, {3, 0, 1, 2});
   }
}

TEST_F(RecordBlocksTest, LossyPoliciesStayInTheRecord)
{
   for (const uint32_t policy :
        std::vector<uint32_t>{RAMNTupleRecord::kDrop, RAMNTupleRecord::kDrop | RAMNTupleRecord::kFqzcomp}) {
      Write({{0, "IIII"}}, /*block=*/3, policy);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      QualityBlockReader quals(*file);
      EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
      EXPECT_EQ(quals.Get(view(0), 0), "*");
   }
}

// Binned qualities from Q30 up encode to printable bytes ('F', Q37, becomes
// byte 37, '%'), which must not be taken for SAM quality and put in a block.
TEST_F(RecordBlocksTest, IlluminaBinnedQualitiesStayInTheRecord)
{
   Write({{0, "FFFFFFFF"}, {0, "#+5?FFII"}}, /*block=*/3, RAMNTupleRecord::kIlluminaBinning);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   for (std::size_t row = 0; row < 2; row++) {
      EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kQualInBlock)) << "row " << row;
      EXPECT_EQ(quals.Get(view(row), row), view(row).GetQUAL()) << "row " << row;
   }
   EXPECT_EQ(quals.Get(view(0), 0), "FFFFFFFF");
}

// A damaged file is reported, not read past: block ends that fall back, and
// an end beyond the last row.
TEST_F(RecordBlocksTest, DamagedBlockEndsAreReported)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "HHHH"}, {0, "GGGG"}};
   for (const auto &ends : {std::vector<uint64_t>{1, 0}, std::vector<uint64_t>{1, 3}}) {
      Write(reads, /*block=*/2, kRansPolicy, ends);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      EXPECT_THROW(QualityBlockReader quals(*file), std::runtime_error) << ends[0] << "," << ends[1];
      EXPECT_THROW(NameBlockReader names(*file), std::runtime_error) << ends[0] << "," << ends[1];
   }
}

// A damaged quality block is reported when it is read: an unknown codec, a
// block cut short, a lengths part longer than the block, and streams that do
// not decode.
TEST_F(RecordBlocksTest, DamagedQualityBlocksAreReported)
{
   const std::vector<Read> reads = {{0, "IIIIHHHGGF"}, {0x10, "#####ABCDE"}};
   const std::vector<std::pair<uint32_t, Damage>> damages = {
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(0) = 9; }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(1); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(4); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(4) = 0x7f; }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.resize(b.size() - 1); }},
      {kRansPolicy, [](std::vector<std::uint8_t> &b) { b.at(1) = 1; }},
      {kFqzcompPolicy, [](std::vector<std::uint8_t> &b) { b.resize(3); }},
   };
   for (std::size_t i = 0; i < damages.size(); i++) {
      Write(reads, /*block=*/2, damages[i].first, /*ends=*/{}, damages[i].second);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
      QualityBlockReader quals(*file);
      EXPECT_THROW(quals.Get(view(0), 0), std::runtime_error) << "damage " << i;
   }
}

// Names go into blocks whatever the quality policy, and come back in any order:
// a repeated name, "*", an empty one, the longest SAM allows.
TEST_F(RecordBlocksTest, EveryNameComesBackInAnyOrder)
{
   const std::vector<Read> reads = {
      {0, "IIII", "A00217:76:HFLT3DSXX:4:1428:30599:8797"},
      {0, "IIII", "A00217:76:HFLT3DSXX:4:1428:30599:8797"},
      {0, "IIII", "*"},
      {0, "IIII", ""},
      {0, "IIII", std::string(254, 'q')},
      {0, "IIII", "SRR123.4"},
      {0, "IIII", "r!#$%&'()*+,-./0;<=>?@[]^_`{|}~"},
   };
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   NameBlockReader names(*file);
   for (const std::size_t row : std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 6, 3, 0, 5, 1})
      EXPECT_EQ(names.Get(row), reads[row].name) << "row " << row;
}

TEST_F(RecordBlocksTest, NamesAreNotAColumn)
{
   Write({{0, "IIII", "a"}}, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   const auto &desc = file->GetReader().GetDescriptor();
   EXPECT_EQ(desc.FindFieldId("qname", desc.FindFieldId("record")), ROOT::kInvalidDescriptorId);
}

// Block ends that put more rows in a block than it has names for.
TEST_F(RecordBlocksTest, ABlockWithTooFewNamesIsReported)
{
   const std::vector<Read> reads = {{0, "IIII", "a"}, {0, "IIII", "b"}, {0, "IIII", "c"}};
   Write(reads, /*block=*/2, RAMNTupleRecord::kPhred33, /*ends=*/{2});
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   NameBlockReader names(*file);
   EXPECT_THROW(names.Get(0), std::runtime_error);
}

} // namespace
