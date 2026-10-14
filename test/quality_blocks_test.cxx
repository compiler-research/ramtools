#include "ramcore/QualityBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
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
};

const char *const kFile = "quality_blocks_test.root";
// The policy samtoramntuple -fqzcomp passes: lossless, in blocks.
const uint32_t kFqzcomp = RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kQualInBlock;

// Writes the reads; \p ends, when given, replaces the block ends, as in a damaged file.
void Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy, const std::vector<uint64_t> &ends = {})
{
   RAMNTupleRecord::InitializeRefs();
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      QualityBlockWriter out(
         writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
         [&writer](ROOT::REntry &e) { writer->Fill(e); }, block);
      for (const auto &read : reads) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetBit(policy);
         rec.SetQNAME("r");
         rec.SetFLAG(read.flag);
         rec.SetRNAME("chr1");
         rec.SetQUAL(read.qual);
         out.Add();
      }
      out.Finish();
      RAMNTupleRecord::SetQualBlockEnds(ends.empty() ? out.TakeBlockEnds() : ends);
   }
   RAMNTupleRecord::WriteAllRefs(*file);
   file->Close();
}

void ExpectQualities(const std::vector<Read> &reads, const std::vector<std::size_t> &order)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   for (const std::size_t row : order)
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
}

class QualityBlocksTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(QualityBlocksTest, EveryQualityComesBackInAnyOrder)
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
   EXPECT_EQ(RAMNTupleRecord::GetQualBlockEnds(), (std::vector<uint64_t>{2, 5, 8, 9}));

   std::vector<std::size_t> forward;
   std::vector<std::size_t> backward;
   for (std::size_t i = 0; i < reads.size(); i++) {
      forward.push_back(i);
      backward.push_back(reads.size() - 1 - i);
   }
   ExpectQualities(reads, forward);
   ExpectQualities(reads, backward);

   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   for (std::size_t row = 0; row < reads.size(); row++) {
      const bool inBlock = reads[row].qual != "*" && reads[row].qual != "AB C";
      EXPECT_EQ(view(row).TestBit(RAMNTupleRecord::kQualInBlock), inBlock) << "row " << row;
   }
}

TEST_F(QualityBlocksTest, ABlockWithoutQualitiesNeedsNoData)
{
   const std::vector<Read> reads = {{0, "*"}, {0, "*"}, {0, "*"}, {0, "HHHH"}};
   Write(reads, /*block=*/3, kFqzcomp);
   ExpectQualities(reads, {3, 0, 1, 2});
}

TEST_F(QualityBlocksTest, LossyPoliciesStayInTheRecord)
{
   Write({{0, "IIII"}}, /*block=*/3, RAMNTupleRecord::kDrop | RAMNTupleRecord::kQualInBlock);
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
   EXPECT_EQ(quals.Get(view(0), 0), "*");
}

// Binned qualities from Q30 up encode to printable bytes ('F', Q37, becomes
// byte 37, '%'), which must not be taken for SAM quality and put in a block.
TEST_F(QualityBlocksTest, IlluminaBinnedQualitiesStayInTheRecord)
{
   Write({{0, "FFFFFFFF"}, {0, "#+5?FFII"}}, /*block=*/3,
         RAMNTupleRecord::kIlluminaBinning | RAMNTupleRecord::kQualInBlock);
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   for (std::size_t row = 0; row < 2; row++) {
      EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kQualInBlock)) << "row " << row;
      EXPECT_EQ(quals.Get(view(row), row), view(row).GetQUAL()) << "row " << row;
   }
   EXPECT_EQ(quals.Get(view(0), 0), "FFFFFFFF");
}

// A damaged file is reported, not read past: block ends that fall back, and
// an end beyond the last row.
TEST_F(QualityBlocksTest, DamagedBlockEndsAreReported)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "HHHH"}, {0, "GGGG"}};
   for (const auto &ends : {std::vector<uint64_t>{1, 0}, std::vector<uint64_t>{1, 3}}) {
      Write(reads, /*block=*/2, kFqzcomp, ends);
      auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
      ASSERT_NE(reader, nullptr);
      EXPECT_THROW(QualityBlockReader quals(*reader), std::runtime_error) << ends[0] << "," << ends[1];
   }
}

// Without -fqzcomp the policy does not ask for blocks, and lossless qualities
// stay in the record as before.
TEST_F(QualityBlocksTest, QualitiesStayInTheRecordUnlessAsked)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0x10, "#+5?"}};
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   for (std::size_t row = 0; row < reads.size(); row++) {
      EXPECT_FALSE(view(row).TestBit(RAMNTupleRecord::kQualInBlock)) << "row " << row;
      EXPECT_EQ(view(row).qual, reads[row].qual) << "row " << row;
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
   }
}

} // namespace
