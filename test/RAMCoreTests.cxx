#include <gtest/gtest.h>
#include "RAMCore/SamToTTree.h"
#include "RAMCore/SamToNTuple.h"
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RNTupleReader.hxx>
#include <cstdio>

class RAMCoreTest : public ::testing::Test {
protected:
    void SetUp() override {
        
        std::remove("test_ttree.root");
        std::remove("test_rntuple.root");
    }
    
    void TearDown() override {
        
        std::remove("test_ttree.root");
        std::remove("test_rntuple.root");
    }
};

TEST_F(RAMCoreTest, ConversionProducesEqualEntries) {
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
    
    auto reader = ROOT::Experimental::RNTupleReader::Open("RAM", rntupleFile);
    ASSERT_NE(reader, nullptr) << "Failed to open RNTuple";
    Long64_t rntupleEntries = reader->GetNEntries();
    
    EXPECT_EQ(ttreeEntries, rntupleEntries) 
        << "Entry count mismatch - TTree: " << ttreeEntries 
        << ", RNTuple: " << rntupleEntries;
    EXPECT_GT(ttreeEntries, 0) << "No entries found";
}

