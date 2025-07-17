#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <ROOT/RNTupleReader.hxx>
#include <cassert>
#include <iostream>


void test_count_entries(const char* samfile = "samexample.sam") {
    
    const char* ttreeFile   = "ttree_test.root";
    const char* rntupleFile = "rntuple_test.root";

    
    gSystem->Unlink(ttreeFile);
    gSystem->Unlink(rntupleFile);

    
    int ret1 = gSystem->Exec(Form("root -l -b -q 'samtoram.C+(\"%s\", \"%s\")'", samfile, ttreeFile));
    
    int ret2 = gSystem->Exec(Form("root -l -b -q 'rntuple/samtoramntuple.C+(\"%s\", \"%s\")'", samfile, rntupleFile));
    assert(ret1 == 0 && ret2 == 0);

    
    TFile* ft = TFile::Open(ttreeFile);
    assert(ft && !ft->IsZombie());
    auto ttree = static_cast<TTree*>(ft->Get("RAM"));
    assert(ttree);
    Long64_t nTT = ttree->GetEntries();
    ft->Close();

    
    auto reader = ROOT::Experimental::RNTupleReader::Open("RAM", rntupleFile);
    assert(reader);
    Long64_t nNT = reader->GetNEntries();

    std::cout << "[INFO] TTree entries   : " << nTT << std::endl;
    std::cout << "[INFO] RNTuple entries : " << nNT << std::endl;

    if (nTT == nNT) {
        std::cout << "[TEST] PASS – entry counts match (" << nTT << ")" << std::endl;
    } else {
        std::cerr << "[TEST] FAIL – counts differ (TTree=" << nTT << ", RNTuple=" << nNT << ")" << std::endl;
        assert(false && "Entry count mismatch");
    }

    
    gSystem->Unlink(ttreeFile);
    gSystem->Unlink(rntupleFile);
} 

