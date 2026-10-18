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
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Read {
   uint16_t flag;
   std::string qual;
   std::string name = "r";
};

const char *const kFile = "record_blocks_test.root";
// The policy samtoramntuple -fqzcomp passes: lossless, in blocks.
const uint32_t kFqzcomp = RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kQualInBlock;

// Writes the reads; \p ends, when given, replaces the block ends, as in a damaged file.
void Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy, const std::vector<uint64_t> &ends = {})
{
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   RAMFileMetadata meta;
   meta.rname_refs = {"chr1"};
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      BlockWriter out(
         writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
         [&writer](ROOT::REntry &e) { writer->Fill(e); }, block);
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
   Write(reads, /*block=*/3, kFqzcomp);
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
   }
}

TEST_F(RecordBlocksTest, ABlockWithoutQualitiesNeedsNoData)
{
   const std::vector<Read> reads = {{0, "*"}, {0, "*"}, {0, "*"}, {0, "HHHH"}};
   Write(reads, /*block=*/3, kFqzcomp);
   ExpectQualities(reads, {3, 0, 1, 2});
}

TEST_F(RecordBlocksTest, LossyPoliciesStayInTheRecord)
{
   Write({{0, "IIII"}}, /*block=*/3, RAMNTupleRecord::kDrop | RAMNTupleRecord::kQualInBlock);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
   EXPECT_EQ(quals.Get(view(0), 0), "*");
}

// Binned qualities from Q30 up encode to printable bytes ('F', Q37, becomes
// byte 37, '%'), which must not be taken for SAM quality and put in a block.
TEST_F(RecordBlocksTest, IlluminaBinnedQualitiesStayInTheRecord)
{
   Write({{0, "FFFFFFFF"}, {0, "#+5?FFII"}}, /*block=*/3,
         RAMNTupleRecord::kIlluminaBinning | RAMNTupleRecord::kQualInBlock);
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
      Write(reads, /*block=*/2, kFqzcomp, ends);
      auto file = RAMFile::Open(kFile);
      ASSERT_NE(file, nullptr);
      EXPECT_THROW(QualityBlockReader quals(*file), std::runtime_error) << ends[0] << "," << ends[1];
      EXPECT_THROW(NameBlockReader names(*file), std::runtime_error) << ends[0] << "," << ends[1];
   }
}

// Without -fqzcomp the policy does not ask for blocks, and lossless qualities
// stay in the record as before.
TEST_F(RecordBlocksTest, QualitiesStayInTheRecordUnlessAsked)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0x10, "#+5?"}};
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto file = RAMFile::Open(kFile);
   ASSERT_NE(file, nullptr);
   auto view = file->GetReader().GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*file);
   for (std::size_t row = 0; row < reads.size(); row++) {
      EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kQualInBlock)) << "row " << row;
      EXPECT_EQ(view(row).qual, reads[row].qual) << "row " << row;
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
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
