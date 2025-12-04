#include <gtest/gtest.h>
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/RAMNTupleView.h"
#include "generate_sam_benchmark.h"
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RNTupleReader.hxx>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "../tools/ramview.cxx"

namespace {

constexpr int kNumReadsForTest = 100;
constexpr const char kSamFile[] = "samexample.sam";
constexpr const char kTTreeFile[] = "test_ttree.root";
constexpr const char kRNTupleFile[] = "test_rntuple.root";

class ramcoreTest : public ::testing::Test {
protected:

   void SetUp() override
   {
      GenerateSAMFile(kSamFile, kNumReadsForTest);

      std::remove(kTTreeFile);
      std::remove(kRNTupleFile);
   }

   void TearDown() override
   {
      std::remove(kTTreeFile);
      std::remove(kRNTupleFile);
   }
};

TEST_F(ramcoreTest, ConversionProducesEqualEntries)
{
   samtoram(kSamFile, kTTreeFile, true, true, true, 1, 0);
   samtoramntuple(kSamFile, kRNTupleFile, true, true, true, 505, 0);

   auto ft = std::unique_ptr<TFile>(TFile::Open(kTTreeFile));
   ASSERT_TRUE(ft && !ft->IsZombie()) << "Failed to open TTree file";

   auto ttree = dynamic_cast<TTree *>(ft->Get("RAM"));
   ASSERT_NE(ttree, nullptr) << "Failed to get TTree";
   Long64_t ttreeEntries = ttree->GetEntries();

   auto reader = ROOT::RNTupleReader::Open("RAM", kRNTupleFile);
   ASSERT_NE(reader, nullptr) << "Failed to open RNTuple";
   Long64_t rntupleEntries = reader->GetNEntries();

   EXPECT_EQ(ttreeEntries, rntupleEntries)
      << "Entry count mismatch - TTree: " << ttreeEntries << ", RNTuple: " << rntupleEntries;
   EXPECT_GT(ttreeEntries, 0) << "No entries found";

   const char *region = "chrM:1-100000000";

   testing::internal::CaptureStdout();
   ramview(kTTreeFile, region, /*cache=*/true, /*perfstats=*/false, /*perfstatsfilename=*/nullptr);
   std::string ramview_output{};
   ramview_output = testing::internal::GetCapturedStdout();

   testing::internal::CaptureStdout();
   ramntupleview(kRNTupleFile, region, /*cache=*/true, /*perfstats=*/false, /*perfstatsfilename=*/nullptr);
   std::string ramntupleview_output{};
   ramntupleview_output = testing::internal::GetCapturedStdout();

   EXPECT_TRUE(ramview_output.find("Found") != std::string::npos);
   EXPECT_TRUE(ramntupleview_output.find("Found") != std::string::npos);
   EXPECT_TRUE(ramntupleview_output.find("records in region") != std::string::npos);
}

} // namespace

