#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include "../benchmark/generate_sam_benchmark.h"
#include "../tools/ramview.cxx"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/SamToTTree.h"
namespace {

class ramcoreTest : public ::testing::Test {
protected:
    void SetUp() override {
       GenerateSAMFile("samexample.sam", 100);
       std::remove("test_ttree.root");
       std::remove("test_rntuple.root");
    }

    void TearDown() override {
        std::remove("test_ttree.root");
        std::remove("test_rntuple.root");
        std::remove("samexample.sam");
    }
};

TEST_F(ramcoreTest, ConversionProducesEqualEntries) {
   const char *samFile = "samexample.sam";
   const char *ttreeFile = "test_ttree.root";
   const char *rntupleFile = "test_rntuple.root";

   samtoram(samFile, ttreeFile, true, true, true, 1, 0);
   samtoramntuple(samFile, rntupleFile, true, true, true, 505, 0);

   auto ft = std::unique_ptr<TFile>(TFile::Open(ttreeFile));
   ASSERT_TRUE(ft && !ft->IsZombie());

   auto ttree = dynamic_cast<TTree *>(ft->Get("RAM"));
   Long64_t ttreeEntries = ttree->GetEntries();

   auto reader = ROOT::RNTupleReader::Open("RAM", rntupleFile);
   Long64_t rntupleEntries = reader->GetNEntries();

   EXPECT_EQ(ttreeEntries, rntupleEntries);
   EXPECT_EQ(ttreeEntries, 100);
}

TEST_F(ramcoreTest, RNTupleViewRegionQueries)
{
   const char *rntupleFile = "test_rntuple.root";
   samtoramntuple("samexample.sam", rntupleFile, true, true, true, 505, 0);

   Long64_t hit = ramntupleview(rntupleFile, "chr1:1-1000000", true, false, nullptr);
   EXPECT_GE(hit, 0);

   Long64_t miss = ramntupleview(rntupleFile, "chrNonExistent:1-100", true, false, nullptr);
   EXPECT_EQ(miss, 0);

   Long64_t wildcard = ramntupleview(rntupleFile, "*", true, false, nullptr);
   EXPECT_EQ(wildcard, 100);

   Long64_t empty = ramntupleview(rntupleFile, "", true, false, nullptr);
   EXPECT_EQ(empty, 100);

   Long64_t null = ramntupleview(rntupleFile, nullptr, true, false, nullptr);
   EXPECT_EQ(null, 100);

   Long64_t whole = ramntupleview(rntupleFile, "chr1", true, false, nullptr);
   EXPECT_GE(whole, 0);

   Long64_t single = ramntupleview(rntupleFile, "chr1:500", true, false, nullptr);
   EXPECT_GE(single, 0);

   Long64_t invalid = ramntupleview(rntupleFile, "chr1:abc-def", true, false, nullptr);
   EXPECT_EQ(invalid, 0);

   Long64_t lateChr = ramntupleview(rntupleFile, "chrX:1-100", true, false, nullptr);
   EXPECT_GE(lateChr, 0);

   Long64_t zeroStart = ramntupleview(rntupleFile, "chr1:0-100", true, false, nullptr);
   EXPECT_GE(zeroStart, 0);
}

TEST_F(ramcoreTest, RNTupleViewOpenFailure)
{
   Long64_t count = ramntupleview("nonexistent_file.root", "chr1:1-100", true, false, nullptr);
   EXPECT_EQ(count, 0);
}

TEST_F(ramcoreTest, RNTupleDataIntegrity)
{
   const char *rntupleFile = "test_rntuple.root";
   samtoramntuple("samexample.sam", rntupleFile, true, true, true, 505, 0);

   auto reader = ROOT::RNTupleReader::Open("RAM", rntupleFile);
   ASSERT_NE(reader, nullptr);

   auto viewPos = reader->GetView<int>("record.pos");
   auto viewSeq = reader->GetView<std::string>("record.seq");

   ASSERT_GT(reader->GetNEntries(), 0);

   int firstPos = viewPos(0);
   std::string firstSeq = viewSeq(0);

   EXPECT_GT(firstPos, 0);

   size_t expectedSize = 18 + 4;
   EXPECT_EQ(firstSeq.size(), expectedSize);

   uint32_t storedLen = 0;
   std::memcpy(&storedLen, firstSeq.data(), sizeof(storedLen));
   EXPECT_EQ(storedLen, 36);

   std::cout << "[   INFO   ] Data Integrity - Pos: " << firstPos << ", Encoded Seq Size: " << firstSeq.size()
             << " (Header: " << storedLen << ")" << std::endl;
}

TEST_F(ramcoreTest, RNTupleViewCigarOverlap)
{
   const char *customSam = "test_cigar.sam";
   const char *rntupleFile = "test_cigar.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      sam << "@SQ\tSN:chr1\tLN:10000\n";
      sam << "cigar_read\t0\tchr1\t100\t60\t10M80D10M\t*\t0\t0\t"
          << "AAAAAAAAAAAAAAAAAAAA\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, false, false, false, 505, 0);

   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:140-160", true, false, nullptr), 1);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:201-210", true, false, nullptr), 0);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:50-110", true, false, nullptr), 1);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:1-50", true, false, nullptr), 0);

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, RNTupleViewFlagFiltering)
{
   const char *customSam = "test_flags.sam";
   const char *rntupleFile = "test_flags.root";
   const std::string seq(50, 'A');

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      sam << "@SQ\tSN:chr1\tLN:1000\n";
      sam << "primary\t0\tchr1\t100\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      sam << "secondary\t256\tchr1\t100\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      sam << "supplementary\t2048\tchr1\t100\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      sam << "unmapped\t4\tchr1\t100\t0\t50M\t*\t0\t0\t" << seq << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, false, false, false, 505, 0);

   Long64_t count = ramntupleview(rntupleFile, "chr1:1-500", true, false, nullptr);
   EXPECT_EQ(count, 1) << "Only primary read should pass FLAG_FILTER";

   std::remove(customSam);
   std::remove(rntupleFile);
}

} // namespace