#include <gtest/gtest.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RNTupleReader.hxx>
#include <iostream>
#include <string>

class CountEntriesTest : public ::testing::Test {
protected:
    std::string samFile;
    std::string ttreeFile;
    std::string rntupleFile;

    void SetUp() override {
        // Setup test files
        samFile = "samexample.sam";
        ttreeFile = "ttree_test.root";
        rntupleFile = "rntuple_test.root";
        
        // Clean up any existing files
        gSystem->Unlink(ttreeFile.c_str());
        gSystem->Unlink(rntupleFile.c_str());
    }

    void TearDown() override {
        // Clean up test files
        gSystem->Unlink(ttreeFile.c_str());
        gSystem->Unlink(rntupleFile.c_str());
    }
};

TEST_F(CountEntriesTest, CompareEntryCountsFromSamFile) {
    // Use gSystem->Exec to run ROOT macros - we're already in tools/ directory
    int ret1 = gSystem->Exec(Form("root -l -b -q 'samtoram.cxx+(\"%s\", \"%s\")'", 
                                  samFile.c_str(), ttreeFile.c_str()));
    ASSERT_EQ(0, ret1) << "Failed to create TTree file";
    
    int ret2 = gSystem->Exec(Form("root -l -b -q 'samtoramntuple.cxx+(\"%s\", \"%s\")'", 
                                  samFile.c_str(), rntupleFile.c_str()));
    ASSERT_EQ(0, ret2) << "Failed to create RNTuple file";
    
    TFile* ft = TFile::Open(ttreeFile.c_str());
    ASSERT_TRUE(ft != nullptr) << "Failed to open TTree file";
    ASSERT_FALSE(ft->IsZombie()) << "TTree file is zombie";
    
    auto ttree = static_cast<TTree*>(ft->Get("RAM"));
    ASSERT_TRUE(ttree != nullptr) << "Failed to get TTree 'RAM'";
    Long64_t nTT = ttree->GetEntries();
    ft->Close();
    
    auto reader = ROOT::Experimental::RNTupleReader::Open("RAM", rntupleFile);
    ASSERT_TRUE(reader != nullptr) << "Failed to open RNTuple";
    Long64_t nNT = reader->GetNEntries();
    
    EXPECT_EQ(nTT, nNT) << "Entry count mismatch - TTree: " << nTT 
                        << ", RNTuple: " << nNT;
    
    EXPECT_GT(nTT, 0) << "TTree has no entries";
    
    std::cout << "[INFO] Successfully verified entry counts match: " << nTT << std::endl;
}

