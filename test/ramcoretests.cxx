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

class ramcoreTest : public ::testing::Test {
protected:
   const int m_num_reads_for_test = 100;
   const std::string m_samFile = "samexample.sam";
   const std::string m_ttreeFile = "test_ttree.root";
   const std::string m_rntupleFile = "test_rntuple.root";

   void SetUp() override
   {
      GenerateSAMFile(m_samFile.c_str(), m_num_reads_for_test);

      std::remove(m_ttreeFile.c_str());
      std::remove(m_rntupleFile.c_str());
   }

   void TearDown() override
   {
      std::remove(m_ttreeFile.c_str());
      std::remove(m_rntupleFile.c_str());
   }
};

TEST_F(ramcoreTest, ConversionProducesEqualEntries)
{
   samtoram(m_samFile.c_str(), m_ttreeFile.c_str(), true, true, true, 1, 0);
   samtoramntuple(m_samFile.c_str(), m_rntupleFile.c_str(), true, true, true, 505, 0);

   auto ft = std::unique_ptr<TFile>(TFile::Open(m_ttreeFile.c_str()));
   ASSERT_TRUE(ft && !ft->IsZombie()) << "Failed to open TTree file";

   auto ttree = dynamic_cast<TTree *>(ft->Get("RAM"));
   ASSERT_NE(ttree, nullptr) << "Failed to get TTree";
   Long64_t ttreeEntries = ttree->GetEntries();

   auto reader = ROOT::RNTupleReader::Open("RAM", m_rntupleFile.c_str());
   ASSERT_NE(reader, nullptr) << "Failed to open RNTuple";
   Long64_t rntupleEntries = reader->GetNEntries();

   EXPECT_EQ(ttreeEntries, rntupleEntries)
      << "Entry count mismatch - TTree: " << ttreeEntries << ", RNTuple: " << rntupleEntries;
   EXPECT_GT(ttreeEntries, 0) << "No entries found";

   const char *region = "chrM:1-100000000";

   testing::internal::CaptureStdout();
   ramview(m_ttreeFile.c_str(), region, /*cache=*/true, /*perfstats=*/false, /*perfstatsfilename=*/nullptr);
   std::string ramview_output{};
   ramview_output = testing::internal::GetCapturedStdout();

   testing::internal::CaptureStdout();
   ramntupleview(m_rntupleFile.c_str(), region, /*cache=*/true, /*perfstats=*/false, /*perfstatsfilename=*/nullptr);
   std::string ramntupleview_output{};
   ramntupleview_output = testing::internal::GetCapturedStdout();

   EXPECT_TRUE(ramview_output.find("Found") != std::string::npos);
   EXPECT_TRUE(ramntupleview_output.find("Found") != std::string::npos);
   EXPECT_TRUE(ramntupleview_output.find("records in region") != std::string::npos);
}

} // namespace
