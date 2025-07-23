#include <TSystem.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <algorithm>


void test_ramview_counts(const char* samfile = "samexample.sam",
                         const char* region  = "chr1:1-10000000") {
    const char* ttreeFile   = "ttree_view.root";
    const char* rntupleFile = "rntuple_view.root";

    gSystem->Unlink(ttreeFile);
    gSystem->Unlink(rntupleFile);
    gSystem->Unlink("ttree_out.txt");
    gSystem->Unlink("rntuple_out.txt");

    int r1 = gSystem->Exec(Form("root -l -b -q 'samtoram.C+(\"%s\", \"%s\")'", samfile, ttreeFile));
    
    int r2 = gSystem->Exec(Form("root -l -b -q 'rntuple/samtoramntuple.C+(\"%s\", \"%s\")'", samfile, rntupleFile));
    assert(r1 == 0 && r2 == 0);

    std::string regStr(region);
    regStr.erase(std::remove_if(regStr.begin(), regStr.end(), ::isspace), regStr.end());

    gSystem->Exec(Form("root -l -b -q 'ramview.C+(\"%s\", \"%s\")' > ttree_out.txt", ttreeFile, regStr.c_str()));
    gSystem->Exec(Form("root -l -b -q 'rntuple/ramntupleview.C+(\"%s\", \"%s\")' > rntuple_out.txt", rntupleFile, regStr.c_str()));

    auto getCount = [](const char* path) -> long long {
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
    };

    long long cTree   = getCount("ttree_out.txt");
    long long cNTuple = getCount("rntuple_out.txt");

    std::cout << "[INFO] ramview records      : " << cTree   << std::endl;
    std::cout << "[INFO] ramntupleview records: " << cNTuple << std::endl;

    if (cTree == cNTuple) {
        std::cout << "[TEST] PASS – viewers report same count (" << cTree << ")" << std::endl;
    } else {
        std::cerr << "[TEST] FAIL – count mismatch" << std::endl;
        assert(false && "Viewer count mismatch");
    }

    gSystem->Unlink(ttreeFile);
    gSystem->Unlink(rntupleFile);
    gSystem->Unlink("ttree_out.txt");
    gSystem->Unlink("rntuple_out.txt");
} 

