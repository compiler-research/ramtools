#include <gtest/gtest.h>
#include <TSystem.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <algorithm>

class RamviewCountsTest : public ::testing::Test {
protected:
    std::string samFile;
    std::string ttreeFile;
    std::string rntupleFile;
    std::string ttreeOutput;
    std::string rntupleOutput;

    void SetUp() override {
        samFile = "samexample.sam";
        ttreeFile = "ttree_view.root";
        rntupleFile = "rntuple_view.root";
        ttreeOutput = "ttree_out.txt";
        rntupleOutput = "rntuple_out.txt";
        
        // Clean up any existing files
        CleanupFiles();
    }

    void TearDown() override {
        CleanupFiles();
    }

    void CleanupFiles() {
        gSystem->Unlink(ttreeFile.c_str());
        gSystem->Unlink(rntupleFile.c_str());
        gSystem->Unlink(ttreeOutput.c_str());
        gSystem->Unlink(rntupleOutput.c_str());
    }

    long long GetCountFromOutput(const char* path) {
        FILE* fp = fopen(path, "r");
        if (!fp) return -1;
        
        char line[256];
        long long cnt = -1;
        while (fgets(line, sizeof(line), fp)) {
            if (strncmp(line, "Found", 5) == 0) {
                sscanf(line, "Found %lld", &cnt);
            }
        }
        fclose(fp);
        return cnt;
    }
};

TEST_F(RamviewCountsTest, CompareViewCountsForRegion) {
    std::string region = "chr1:1-10000000";
    
    // We're already in tools/ directory, so no path needed
    int r1 = gSystem->Exec(Form("root -l -b -q 'samtoram.cxx+(\"%s\", \"%s\")'", 
                                samFile.c_str(), ttreeFile.c_str()));
    ASSERT_EQ(0, r1) << "Failed to create TTree file";
    
    int r2 = gSystem->Exec(Form("root -l -b -q 'samtoramntuple.cxx+(\"%s\", \"%s\")'", 
                                samFile.c_str(), rntupleFile.c_str()));
    ASSERT_EQ(0, r2) << "Failed to create RNTuple file";
    
    region.erase(std::remove_if(region.begin(), region.end(), ::isspace), region.end());
    
    gSystem->Exec(Form("root -l -b -q 'ramview.cxx+(\"%s\", \"%s\")' > %s", 
                      ttreeFile.c_str(), region.c_str(), ttreeOutput.c_str()));
    
    gSystem->Exec(Form("root -l -b -q 'ramntupleview.cxx+(\"%s\", \"%s\")' > %s", 
                      rntupleFile.c_str(), region.c_str(), rntupleOutput.c_str()));
    
    long long ttreeCount = GetCountFromOutput(ttreeOutput.c_str());
    long long rntupleCount = GetCountFromOutput(rntupleOutput.c_str());
    
    EXPECT_GT(ttreeCount, 0) << "No records found by ramview";
    EXPECT_GT(rntupleCount, 0) << "No records found by ramntupleview";
    
    EXPECT_EQ(ttreeCount, rntupleCount) 
        << "View count mismatch - ramview: " << ttreeCount 
        << ", ramntupleview: " << rntupleCount;
    
    std::cout << "[INFO] Successfully verified view counts match: " << ttreeCount 
              << " records in region " << region << std::endl;
}

