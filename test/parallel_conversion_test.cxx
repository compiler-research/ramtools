// Many small blocks on several threads must equal one block on one thread.
#include "ramcore/QualityBlocks.h"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr size_t kTinyBlock = 700;

std::string Header(const char *order = "coordinate")
{
   std::string h{};
   h += "@HD\tVN:1.6\tSO:";
   h += order;
   h += "\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n@SQ\tSN:chr3\tLN:100000\n";
   h += "@PG\tID:test\tPN:test\n";
   return h;
}

std::string
Record(int n, const char *rname, int flag, int pos, const char *rnext = "*", int pnext = 0, const char *cigar = "50M")
{
   const bool mapped = !(flag & 4);
   std::string r = "r" + std::to_string(n) + "\t" + std::to_string(flag) + "\t" + rname + "\t" + std::to_string(pos) +
                   "\t60\t" + (mapped ? cigar : "*") + "\t" + rnext + "\t" + std::to_string(pnext) + "\t0\t" +
                   (mapped ? std::string(50, 'A') : "*") + "\t" + (mapped ? std::string(50, 'F') : "*") +
                   "\tNM:i:" + std::to_string(n % 3) + "\tRG:Z:g" + std::to_string(n % 2) + "\n";
   return r;
}

std::string SortedRecords(int per_ref)
{
   std::string s{};
   int n = 0;
   for (int i = 0; i < per_ref; ++i)
      s += Record(n++, /*rname=*/"chr1", (i % 2) ? 99 : 0, 1 + (i * 37), (i % 5) ? "=" : "chr2", 1 + (i * 41));
   s += Record(n++, /*rname=*/"chr1", /*flag=*/4,
               1 + (per_ref * 37) + 100); // placed unmapped read after the last mapped one
   for (int i = 0; i < per_ref; ++i)
      s += Record(n++, /*rname=*/"chr2", /*flag=*/0, 5 + (i * 41), (i % 7) ? "=" : "chr1", 1 + (i * 37));
   s += Record(n++, /*rname=*/"*", /*flag=*/4, /*pos=*/0);
   s += Record(n++, /*rname=*/"*", /*flag=*/4, /*pos=*/0);
   return s;
}

void WriteSam(const char *path, const std::string &content)
{
   std::ofstream sam(path);
   sam << content;
}

struct Dump {
   std::vector<std::string> header;
   std::vector<std::string> lines;
   bool sorted = false;
   uint32_t max_span = 0;
   std::vector<std::string> rname_refs;
};

Dump ReadBack(const char *file)
{
   Dump d;
   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   EXPECT_NE(reader, nullptr) << file;
   if (!reader)
      return d;
   d.header = RAMNTupleRecord::GetHeaderLines();
   d.sorted = RAMNTupleRecord::IsCoordinateSorted();
   d.max_span = RAMNTupleRecord::GetMaxRefSpan();
   d.rname_refs = RAMNTupleRecord::GetRnameRefs()->GetRefs();
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   for (auto i : reader->GetEntryRange()) {
      const RAMNTupleRecord &r = view(i);
      std::string line = r.GetQNAME() + "\t" + std::to_string(r.GetFLAG()) + "\t" + r.GetRNAME() + "\t" +
                         std::to_string(r.GetPOS()) + "\t" + std::to_string(r.GetMAPQ()) + "\t" + r.GetCIGAR() + "\t" +
                         r.GetRNEXT() + "\t" + std::to_string(r.GetPNEXT()) + "\t" + std::to_string(r.GetTLEN()) +
                         "\t" + r.GetSEQ() + "\t" + quals.Get(r, i);
      for (const auto &t : r.GetTags())
         line += "\t" + t;
      d.lines.push_back(line);
   }
   return d;
}

class ParallelConversionTest : public ::testing::Test {
protected:
   static constexpr const char *kSam = "parallel_test.sam";
   static constexpr const char *kSeq = "parallel_test_seq.root";
   static constexpr const char *kPar = "parallel_test_par.root";

   void TearDown() override
   {
      std::remove(kSam);
      std::remove(kSeq);
      std::remove(kPar);
   }
};

// With -fqzcomp every input block ends its quality blocks, and the block ends
// gathered from all threads have to line up with the rows. Each read has its
// own qualities, so a quality read back from the wrong row shows.
TEST_F(ParallelConversionTest, FqzcompOnSeveralThreadsKeepsEveryQuality)
{
   std::string sam = Header();
   for (int n = 0; n < 200; ++n) {
      std::string qual(50, '!');
      for (size_t i = 0; i < qual.size(); ++i)
         qual[i] = static_cast<char>('!' + (((n * 7) + static_cast<int>(i)) % 60));
      sam += "r" + std::to_string(n) + "\t0\tchr1\t" + std::to_string(1 + (n * 37)) + "\t60\t50M\t*\t0\t0\t" +
             std::string(50, 'A') + "\t" + ((n % 9 == 0) ? std::string("*") : qual) + "\n";
   }
   WriteSam(kSam, sam);
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   const Dump plain = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505,
                  /*quality_policy=*/RAMNTupleRecord::kPhred33 | RAMNTupleRecord::kQualInBlock, /*threads=*/3,
                  kTinyBlock);
   const Dump fqz = ReadBack(kPar);
   EXPECT_GT(RAMNTupleRecord::GetQualBlockEnds().size(), 10U) << "many input blocks, so many quality blocks";

   ASSERT_EQ(fqz.lines.size(), 200U);
   EXPECT_EQ(fqz.lines, plain.lines) << "every quality comes back from its block, on its row";
}

TEST_F(ParallelConversionTest, ManyBlocksOnSeveralThreadsMatchOneBlockRecordForRecord)
{
   WriteSam(kSam, Header() + SortedRecords(/*per_ref=*/300));
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/3, kTinyBlock);
   const Dump par = ReadBack(kPar);

   ASSERT_EQ(par.lines.size(), 603U);
   EXPECT_EQ(par.lines, seq.lines) << "same records in the same order as with one block";
   EXPECT_TRUE(par.sorted);
   EXPECT_EQ(par.sorted, seq.sorted);
   EXPECT_EQ(par.max_span, seq.max_span);
   EXPECT_EQ(par.max_span, 50U);
   EXPECT_EQ(par.rname_refs, seq.rname_refs) << "reference ids follow the header in both";

   auto reader = ROOT::RNTupleReader::Open("RAM", kPar);
   EXPECT_GT(reader->GetDescriptor().GetNClusters(), 20U);
   const RAMNTupleViewOpts opts = {true, false, ""};
   EXPECT_EQ(ramntupleview(kPar, "chr1:11150-11250", opts), ramntupleview(kSeq, "chr1:11150-11250", opts));
   EXPECT_EQ(ramntupleview(kPar, "chr1:11150-11250", opts), 1) << "the placed unmapped read";
   EXPECT_EQ(ramntupleview(kPar, "chr2:4000-8200", opts), ramntupleview(kSeq, "chr2:4000-8200", opts));
   EXPECT_EQ(ramntupleview(kPar, "chr2", opts), 300);
   EXPECT_EQ(ramntupleview(kPar, "chr3", opts), 0);
}

TEST_F(ParallelConversionTest, OneThreadAndOneBlockGiveOneCluster)
{
   WriteSam(kSam, Header() + SortedRecords(/*per_ref=*/40));
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/1);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);

   auto reader = ROOT::RNTupleReader::Open("RAM", kPar);
   EXPECT_EQ(reader->GetDescriptor().GetNClusters(), 1U);
}

// Each block is sorted on its own; only the step between blocks goes back.
TEST_F(ParallelConversionTest, DisorderAcrossABlockBoundaryMarksTheFileUnsorted)
{
   std::string s = Header(/*order=*/"unsorted");
   int n = 0;
   for (int i = 0; i < 12; ++i)
      s += Record(n++, /*rname=*/"chr1", /*flag=*/0, 5000 + (i * 10)); // 12 records at ~150 bytes each: a few blocks
   for (int i = 0; i < 12; ++i)
      s +=
         Record(n++, /*rname=*/"chr1", /*flag=*/0, 1000 + (i * 10)); // sorted within itself, but before the block above
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/2, kTinyBlock);
   EXPECT_NE(testing::internal::GetCapturedStderr().find("not in coordinate order"), std::string::npos);

   const Dump par = ReadBack(kPar);
   EXPECT_FALSE(par.sorted);
   ASSERT_EQ(par.lines.size(), 24U);
   EXPECT_EQ(par.lines.front().substr(0, 3), "r0\t");
   EXPECT_EQ(par.lines.back().substr(0, 4), "r23\t");
   const RAMNTupleViewOpts opts = {true, false, ""};
   EXPECT_EQ(ramntupleview(kPar, "chr1:1000-1200", opts), 12);
}

TEST_F(ParallelConversionTest, UnplacedRecordsBeforePlacedOnesInALaterBlockMeanUnsorted)
{
   std::string s = Header(/*order=*/"unsorted");
   int n = 0;
   for (int i = 0; i < 10; ++i)
      s += Record(n++, /*rname=*/"*", /*flag=*/4, /*pos=*/0);
   for (int i = 0; i < 10; ++i)
      s += Record(n++, /*rname=*/"chr1", /*flag=*/0, 1000 + (i * 10));
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/2, kTinyBlock);
   testing::internal::GetCapturedStderr();

   const Dump par = ReadBack(kPar);
   EXPECT_FALSE(par.sorted);
   EXPECT_EQ(par.lines.size(), 20U);
}

TEST_F(ParallelConversionTest, HeaderLongerThanABlockIsReadBeforeTheRecords)
{
   std::string s = Header();
   for (int i = 0; i < 40; ++i)
      s += "@CO\tcomment line number " + std::to_string(i) + " padding padding padding padding\n";
   s += SortedRecords(/*per_ref=*/20);
   WriteSam(kSam, s);
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/2, kTinyBlock);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_EQ(par.rname_refs, seq.rname_refs);

   EXPECT_EQ(par.header, seq.header);
   ASSERT_EQ(par.header.size(), 45U) << "@HD, three @SQ, @PG and forty @CO lines";
   EXPECT_EQ(par.header.front(), "@HD\tVN:1.6\tSO:coordinate");
   EXPECT_EQ(par.header.back(), "@CO\tcomment line number 39 padding padding padding padding");
}

TEST_F(ParallelConversionTest, MalformedAndEmptyLinesAreSkippedInEveryBlock)
{
   std::string s = Header();
   s += Record(/*n=*/0, /*rname=*/"chr1", /*flag=*/0, /*pos=*/100);
   s += "\n";
   s += "bad\tx\tchr1\t200\t60\t50M\t*\t0\t0\t*\t*\n"; // flag is not a number
   s += Record(/*n=*/1, /*rname=*/"chr1", /*flag=*/0, /*pos=*/300);
   s += "short\t0\tchr1\n";
   s += Record(/*n=*/2, /*rname=*/"chr1", /*flag=*/0, /*pos=*/400, /*rnext=*/"*", /*pnext=*/0,
               /*cigar=*/"10M5Q35M"); // bad CIGAR operator
   s += Record(/*n=*/3, /*rname=*/"chr1", /*flag=*/0, /*pos=*/500);
   s += "@CO\ta comment among the records\n";
   s += Record(/*n=*/4, /*rname=*/"chr1", /*flag=*/0, /*pos=*/600);
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/2, kTinyBlock);
   testing::internal::GetCapturedStderr();

   const Dump seq = ReadBack(kSeq);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(seq.lines.size(), 4U);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);

   EXPECT_EQ(par.header, seq.header);
   ASSERT_EQ(par.header.size(), 6U) << "the late @CO line is kept";
   EXPECT_EQ(par.header.back(), "@CO\ta comment among the records");
}

TEST_F(ParallelConversionTest, LastLineWithoutNewlineIsKept)
{
   std::string s = Header() + SortedRecords(/*per_ref=*/30);
   s.pop_back();
   WriteSam(kSam, s);
   samtoramntuple(kSam, kSeq, /*compression_algorithm=*/505, /*quality_policy=*/0);
   const Dump seq = ReadBack(kSeq);
   ASSERT_EQ(seq.lines.size(), 63U);

   samtoramntuple(kSam, kPar, /*compression_algorithm=*/505, /*quality_policy=*/0, /*threads=*/2, kTinyBlock);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);
}

// A directory opens but fails on read.
TEST_F(ParallelConversionTest, ReadErrorIsReportedAfterTheWorkersAreStopped)
{
   testing::internal::CaptureStdout();
   EXPECT_FALSE(samtoramntuple(".", kPar, 505, 0, 3, kTinyBlock));
   EXPECT_NE(testing::internal::GetCapturedStdout().find("Failed to read SAM file"), std::string::npos);
}

TEST_F(ParallelConversionTest, MissingInputIsReported)
{
   testing::internal::CaptureStdout();
   EXPECT_FALSE(samtoramntuple("does_not_exist.sam", kPar, 505, 0, 2));
   EXPECT_NE(testing::internal::GetCapturedStdout().find("Failed to parse SAM file"), std::string::npos);
}

} // namespace
