#include <gtest/gtest.h>
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include "generate_sam_benchmark.h"
#include <ROOT/RNTupleReader.hxx>
#include <cstdio>
#include <filesystem>

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
   samtoramntuple("test.sam", "test_regular.root", false, true, true, 505, 1);
   auto regularReader = ROOT::Experimental::RNTupleReader::Open("RAM", "test_regular.root");
   Long64_t totalEntries = regularReader->GetNEntries();

   samtoramntuple_split_by_chromosome("test.sam", "test_split", 505, 1);

   Long64_t splitEntriesSum = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      std::string filename = entry.path().filename().string();
      if (filename.find("test_split_") == 0 && filename.find(".root") != std::string::npos) {
         auto reader = ROOT::Experimental::RNTupleReader::Open("RAM", filename);
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

         auto reader = ROOT::Experimental::RNTupleReader::Open("RAM", filename);
         ASSERT_NE(reader, nullptr);

         auto viewRecord = reader->GetView<RAMNTupleRecord>("record");

         for (auto i : reader->GetEntryRange()) {
            const auto &record = viewRecord(i);
            std::string actualChr = record.GetRNAME();
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
         auto metaReader = ROOT::Experimental::RNTupleReader::Open("METADATA", filename);
         EXPECT_NE(metaReader, nullptr) << "Missing METADATA in " << filename;

         if (metaReader) {
            EXPECT_GT(metaReader->GetNEntries(), 0);
         }

         filesChecked++;
      }
   }

   EXPECT_GT(filesChecked, 0);
}
