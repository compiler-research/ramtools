#include <gtest/gtest.h>
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include "generate_sam_benchmark.h"
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RNTupleReader.hxx>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

class ramcoreTest : public ::testing::Test {
protected:
    void SetUp() override {
        const int num_reads_for_test = 100;
        if (!std::filesystem::exists("samexample.sam")) {
            GenerateSAMFile("samexample.sam", num_reads_for_test);
        }

        std::remove("test_ttree.root");
        std::remove("test_rntuple.root");
    }

    void TearDown() override {
        std::remove("test_ttree.root");
        std::remove("test_rntuple.root");
    }
};

TEST_F(ramcoreTest, ConversionProducesEqualEntries) {
    const char* samFile = "samexample.sam";
    const char* ttreeFile = "test_ttree.root";
    const char* rntupleFile = "test_rntuple.root";

    samtoram(samFile, ttreeFile, true, true, true, 1, 0);

    samtoramntuple(samFile, rntupleFile, true, true, true, 505, 0);

    auto ft = std::unique_ptr<TFile>(TFile::Open(ttreeFile));
    ASSERT_TRUE(ft && !ft->IsZombie()) << "Failed to open TTree file";

    auto ttree = dynamic_cast<TTree*>(ft->Get("RAM"));
    ASSERT_NE(ttree, nullptr) << "Failed to get TTree";
    Long64_t ttreeEntries = ttree->GetEntries();

    auto reader = ROOT::RNTupleReader::Open("RAM", rntupleFile);
    ASSERT_NE(reader, nullptr) << "Failed to open RNTuple";
    Long64_t rntupleEntries = reader->GetNEntries();

    EXPECT_EQ(ttreeEntries, rntupleEntries)
        << "Entry count mismatch - TTree: " << ttreeEntries
        << ", RNTuple: " << rntupleEntries;
    EXPECT_GT(ttreeEntries, 0) << "No entries found";
}

