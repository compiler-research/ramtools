#include <gtest/gtest.h>
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMFile.h"
#include "rntuple/RAMNTupleRecord.h"
#include "generate_sam_benchmark.h"
#include <ROOT/RNTupleReader.hxx>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>

class ChromosomeSplitTest : public ::testing::Test {
protected:
   void SetUp() override
   {
      if (!std::filesystem::exists("test.sam")) {
         GenerateSAMFile("test.sam", 100);
      }
      CleanupTestFiles();
   }

   void TearDown() override { CleanupTestFiles(); }

   void CleanupTestFiles()
   {
      std::remove("test_regular.root");
      for (const auto &entry : std::filesystem::directory_iterator(".")) {
         std::string filename = entry.path().filename().string();
         if (filename.find("test_split_") == 0 && filename.find(".root") != std::string::npos) {
            std::remove(filename.c_str());
         }
      }
   }
};

TEST_F(ChromosomeSplitTest, NoDataLoss)
{
   samtoramntuple(/*datafile=*/"test.sam", /*treefile=*/"test_regular.root", /*compression_algorithm=*/505,
                  /*quality_policy=*/1);
   auto regularReader = ROOT::RNTupleReader::Open("RAM", "test_regular.root");
   Long64_t totalEntries = regularReader->GetNEntries();

   samtoramntuple_split_by_chromosome("test.sam", "test_split", 505, 1);

   Long64_t splitEntriesSum = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      std::string filename = entry.path().filename().string();
      if (filename.find("test_split_") == 0 && filename.find(".root") != std::string::npos) {
         auto reader = ROOT::RNTupleReader::Open("RAM", filename);
         if (reader) {
            splitEntriesSum += reader->GetNEntries();
         }
      }
   }

   EXPECT_EQ(totalEntries, splitEntriesSum);
}

TEST_F(ChromosomeSplitTest, CorrectChromosomeAssignment)
{
   samtoramntuple_split_by_chromosome("test.sam", "test_split", 505, 1);

   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      std::string filename = entry.path().filename().string();
      if (filename.find("test_split_") == 0 && filename.find(".root") != std::string::npos) {
         size_t pos = filename.find("test_split_");
         size_t end = filename.find(".root");
         std::string expectedChr = filename.substr(pos + 11, end - pos - 11);

         auto file = RAMFile::Open(filename);
         ASSERT_NE(file, nullptr);

         auto viewRecord = file->GetReader().GetView<RAMNTupleRecord>("record");

         for (auto i : file->GetReader().GetEntryRange()) {
            const auto &record = viewRecord(i);
            std::string actualChr = file->GetRname(record.GetREFID());
            EXPECT_EQ(expectedChr, actualChr) << "Wrong chromosome in " << filename << " at entry " << i;
         }
      }
   }
}

TEST_F(ChromosomeSplitTest, MetadataPresent)
{
   samtoramntuple_split_by_chromosome("test.sam", "test_split", 505, 1);

   int filesChecked = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      std::string filename = entry.path().filename().string();
      if (filename.find("test_split_") == 0 && filename.find(".root") != std::string::npos) {
         auto metaReader = ROOT::RNTupleReader::Open("METADATA", filename);
         EXPECT_NE(metaReader, nullptr) << "Missing METADATA in " << filename;

         if (metaReader) {
            EXPECT_GT(metaReader->GetNEntries(), 0);
         }

         filesChecked++;
      }
   }

   EXPECT_GT(filesChecked, 0);
}

TEST_F(ChromosomeSplitTest, RegionCountsMatchTheUnsplitFile)
{
   samtoramntuple(/*datafile=*/"test.sam", /*treefile=*/"test_regular.root", /*compression_algorithm=*/505,
                  /*quality_policy=*/1);
   samtoramntuple_split_by_chromosome("test.sam", "test_split", 505, 1);

   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      std::string filename = entry.path().filename().string();
      if (filename.find("test_split_") != 0 || filename.find(".root") == std::string::npos)
         continue;
      std::string chr = filename.substr(11, filename.size() - 11 - 5);
      EXPECT_EQ(ramntupleview(filename.c_str(), chr.c_str()), ramntupleview("test_regular.root", chr.c_str()))
         << filename;
   }
}

// Opening each chromosome's writer constructs a record. That used to reset the
// longest span seen so far, so earlier chromosomes' files were written with a
// span too small for the region seek to back off enough.
TEST_F(ChromosomeSplitTest, SplitFilesKeepTheLongestSpan)
{
   const char *spanSam = "test_split_span.sam";
   {
      std::ofstream sam(spanSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:1000000\n";
      sam << "spanning\t0\tchr1\t1000\t60\t10M199990N10M\t*\t0\t0\t" << std::string(20, 'A') << "\t*\n";
      sam << "short\t0\tchr1\t300000\t60\t50M\t*\t0\t0\t" << std::string(50, 'C') << "\t*\n";
      sam << "other\t0\tchr2\t100\t60\t4M\t*\t0\t0\tACGT\t*\n";
   }
   samtoramntuple_split_by_chromosome(spanSam, "test_split", 505, 1);

   const auto file = RAMFile::Open("test_split_chr1.root");
   ASSERT_NE(file, nullptr);
   EXPECT_GE(file->GetMetadata().max_ref_span, 200010U);
   EXPECT_EQ(ramntupleview("test_split_chr1.root", "chr1:200000-201000"), 1) << "the spanning read reaches the region";

   std::remove(spanSam);
}

// Each file counts only its own records, as samtools idxstats would on that
// file; the other references are listed with none.
TEST_F(ChromosomeSplitTest, SplitFilesCountTheirOwnRecords)
{
   const char *countSam = "test_split_counts.sam";
   {
      std::ofstream sam(countSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";
      sam << "a\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t*\n";
      sam << "b\t256\tchr1\t200\t60\t4M\t*\t0\t0\tACGT\t*\n";
      sam << "c\t4\tchr1\t200\t0\t*\t*\t0\t0\tACGT\t*\n";
      sam << "d\t0\tchr2\t100\t60\t4M\t*\t0\t0\tACGT\t*\n";
      sam << "e\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n";
   }
   samtoramntuple_split_by_chromosome(countSam, /*output_prefix=*/"test_split", /*compression_algorithm=*/505,
                                      /*quality_policy=*/1);

   const auto chr1 = RAMFile::Open("test_split_chr1.root");
   ASSERT_NE(chr1, nullptr);
   EXPECT_EQ(chr1->GetMetadata().counts.mapped, (std::vector<uint64_t>{2, 0}));
   EXPECT_EQ(chr1->GetMetadata().counts.unmapped, (std::vector<uint64_t>{1, 0}));
   EXPECT_EQ(chr1->GetMetadata().counts.unplaced, 0U) << "-split writes no file for unplaced records";

   const auto chr2 = RAMFile::Open("test_split_chr2.root");
   ASSERT_NE(chr2, nullptr);
   EXPECT_EQ(chr2->GetMetadata().counts.mapped, (std::vector<uint64_t>{0, 1}));
   EXPECT_EQ(chr2->GetMetadata().counts.unmapped, (std::vector<uint64_t>{0, 0}));

   std::remove(countSam);
}
