#include "ramcore/QualityBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Read {
   uint16_t flag;
   std::string qual;
};

const char *const kFile = "quality_blocks_test.root";

void Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy)
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
      RAMNTupleRecord::SetQualBlockEnds(out.TakeBlockEnds());
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
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
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
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   ExpectQualities(reads, {3, 0, 1, 2});
}

TEST_F(QualityBlocksTest, LossyPoliciesStayInTheRecord)
{
   Write({{0, "IIII"}}, /*block=*/3, RAMNTupleRecord::kDrop);
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
   EXPECT_EQ(quals.Get(view(0), 0), "*");
}

} // namespace
