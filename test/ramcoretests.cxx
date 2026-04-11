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
#include "rntuple/RAMNTupleRecord.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamParser.h"
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
        RAMNTupleRecord::GetIndex()->Clear();
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

TEST_F(ramcoreTest, IndexGetRowsInRange)
{
   RAMNTupleRecord::InitializeRefs();
   auto *index = RAMNTupleRecord::GetIndex();

   //  test with hard coded entries
   index->AddItem(/*refid=*/0, /*pos=*/100, /*row=*/0);
   index->AddItem(/*refid=*/0, /*pos=*/200, /*row=*/1);
   index->AddItem(/*refid=*/0, /*pos=*/300, /*row=*/2);
   index->AddItem(/*refid=*/1, /*pos=*/150, /*row=*/3);

   auto rows = index->GetRowsInRange(/*refid=*/0, /*start=*/150, /*end=*/250);
   ASSERT_EQ(rows.size(), 1U);
   EXPECT_EQ(rows[0], 1);

   auto all = index->GetRowsInRange(/*refid=*/0, /*start=*/0, /*end=*/400);
   EXPECT_EQ(all.size(), 3U);

   auto none = index->GetRowsInRange(/*refid=*/0, /*start=*/400, /*end=*/500);
   EXPECT_TRUE(none.empty());

   auto otherChrom = index->GetRowsInRange(/*refid=*/1, /*start=*/100, /*end=*/200);
   ASSERT_EQ(otherChrom.size(), 1U);
   EXPECT_EQ(otherChrom[0], 3);

   // test with generated entries
   const char *mockFile = "test_mock_index.root";
   samtoramntuple(/*datafile=*/"samexample.sam", mockFile, /*index=*/true, /*split=*/true, /*cache=*/true,
                  /*compression_algorithm=*/505, /*quality_policy=*/0);

   auto reader = RAMNTupleRecord::OpenRAMFile(mockFile);
   ASSERT_NE(reader, nullptr);
   EXPECT_GT(index->Size(), 0U);

   int chr1_refid = RAMNTupleRecord::GetRnameRefs()->GetRefId("chr1");
   EXPECT_GE(chr1_refid, 0);

   auto wideRows = index->GetRowsInRange(/*refid=*/chr1_refid, /*start=*/0, /*end=*/1000000000);
   for (int64_t row : wideRows) {
      EXPECT_GE(row, 0);
      EXPECT_LT(row, 100); // record length of samexample.sam is 100 in SetUp()
   }

   auto invalidRows = index->GetRowsInRange(/*refid=*/-1, /*start=*/0, /*end=*/1000000000);
   EXPECT_TRUE(invalidRows.empty());

   std::remove(mockFile);
}

} // namespace

TEST_F(ramcoreTest, SmartIndexSkipsUnmappedReads)
{
   const char *customSam = "test_unmapped_index.sam";
   const char *rntupleFile = "test_unmapped_index.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      sam << "mapped1\t0\tchr1\t1000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "mapped2\t0\tchr1\t2000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);

   Long64_t count = ramntupleview(rntupleFile, "chr1:900-2100", true, false, nullptr);
   EXPECT_EQ(count, 2) << "Both mapped reads should be queryable";

   Long64_t unmapped = ramntupleview(rntupleFile, "*:0-100", true, false, nullptr);
   EXPECT_EQ(unmapped, 0) << "Unmapped reads should not appear in index queries";

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, SmartIndexCreatesEntryAtChromosomeBoundary)
{
   const char *customSam = "test_chrom_boundary.sam";
   const char *rntupleFile = "test_chrom_boundary.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      sam << "@SQ\tSN:chr2\tLN:100000\n";
      for (int i = 0; i < 50; ++i)
         sam << "chr1_r" << i << "\t0\tchr1\t" << (1000 + i * 100) << "\t60\t50M\t*\t0\t0\t" << std::string(50, 'A')
             << "\t*\n";
      for (int i = 0; i < 50; ++i)
         sam << "chr2_r" << i << "\t0\tchr2\t" << (500 + i * 100) << "\t60\t50M\t*\t0\t0\t" << std::string(50, 'A')
             << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);

   Long64_t chr1_hits = ramntupleview(rntupleFile, "chr1:1000-6000", true, false, nullptr);
   EXPECT_GT(chr1_hits, 0) << "chr1 reads should be queryable";

   Long64_t chr2_hits = ramntupleview(rntupleFile, "chr2:500-5500", true, false, nullptr);
   EXPECT_GT(chr2_hits, 0) << "chr2 reads should be immediately queryable at boundary";

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, SmartIndexRespectsPositionInterval)
{
   const char *customSam = "test_pos_interval.sam";
   const char *rntupleFile = "test_pos_interval.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000000\n";
      // Cluster of 200 reads at position 1000 (same position, should not generate 200 index entries)
      for (int i = 0; i < 200; ++i)
         sam << "cluster_r" << i << "\t0\tchr1\t1000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      // Read far away at 50000 (>10kb gap, should trigger position interval index)
      sam << "far_read\t0\tchr1\t50000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);

   Long64_t cluster = ramntupleview(rntupleFile, "chr1:900-1100", true, false, nullptr);
   EXPECT_EQ(cluster, 200);

   Long64_t far = ramntupleview(rntupleFile, "chr1:49900-50100", true, false, nullptr);
   EXPECT_EQ(far, 1) << "Distant read should be indexed via position interval";

   std::remove(customSam);
   std::remove(rntupleFile);
}
TEST_F(ramcoreTest, MalformedFlagSkipsRecord) {
    const char* path = "/tmp/ram_test_flag.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "bad\tabc\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fprintf(f, "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fclose(f);

    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord&, int) { count++; });
    std::remove(path);

    EXPECT_EQ(count, 1);
}

TEST_F(ramcoreTest, MalformedPosSkipsRecord) {
    const char* path = "/tmp/ram_test_pos.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "bad\t0\tchr1\tXYZ\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fprintf(f, "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fclose(f);

    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord&, int) { count++; });
    std::remove(path);

    EXPECT_EQ(count, 1);
}

TEST_F(ramcoreTest, MalformedMapqSkipsRecord) {
    const char* path = "/tmp/ram_test_mapq.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "bad\t0\tchr1\t100\tBAD\t10M\t*\t0\t0\tACGT\tIIII\n");
    fprintf(f, "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fclose(f);

    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord&, int) { count++; });
    std::remove(path);

    EXPECT_EQ(count, 1);
}

TEST_F(ramcoreTest, MalformedPnextSkipsRecord) {
    const char* path = "/tmp/ram_test_pnext.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "bad\t0\tchr1\t100\t60\t10M\t*\tNOP\t0\tACGT\tIIII\n");
    fprintf(f, "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fclose(f);

    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord&, int) { count++; });
    std::remove(path);

    EXPECT_EQ(count, 1);
}

TEST_F(ramcoreTest, MalformedTlenSkipsRecord) {
    const char* path = "/tmp/ram_test_tlen.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "bad\t0\tchr1\t100\t60\t10M\t*\t0\tBAD\tACGT\tIIII\n");
    fprintf(f, "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII\n");
    fclose(f);

    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord&, int) { count++; });
    std::remove(path);

    EXPECT_EQ(count, 1);
}

TEST_F(ramcoreTest, ValidRecordParsesCorrectly) {
    const char* path = "/tmp/ram_test_valid.sam";
    FILE* f = fopen(path, "w");
    fprintf(f, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n");
    fprintf(f, "r1\t16\tchr1\t500\t30\t5M\t*\t200\t10\tACGTA\tIIIII\n");
    fclose(f);

    ramcore::SamRecord captured;
    int count = 0;
    ramcore::SamParser parser;
    parser.ParseFile(path,
        [](const std::string&, const std::string&) {},
        [&](const ramcore::SamRecord& rec, int) { captured = rec; count++; });
    std::remove(path);

    ASSERT_EQ(count, 1);
    EXPECT_EQ(captured.flag,  16);
    EXPECT_EQ(captured.pos,   500);
    EXPECT_EQ(captured.mapq,  30);
    EXPECT_EQ(captured.pnext, 200);
    EXPECT_EQ(captured.tlen,  10);
}
