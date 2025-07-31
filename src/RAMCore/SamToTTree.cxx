#include "RAMCore/SamToTTree.h"
#include "RAMCore/SamParser.h"
#include "ttree/RAMRecord.h"
#include "ttree/Utils.h"

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TNamed.h>
#include <TStopwatch.h>
#include <Compression.h>

void samtoram(const char *datafile,
              const char *treefile,
              bool index, bool split, bool cache,
              Int_t compression_algorithm,
              UInt_t quality_policy)
{
    TStopwatch stopwatch;
    stopwatch.Start();
    
    auto f = TFile::Open(treefile, "RECREATE");
    if (!f) {
        printf("Failed to create output file %s\n", treefile);
        return;
    }
    
    f->SetCompressionLevel(1);
    f->SetCompressionAlgorithm(compression_algorithm);
    
    auto tree = new TTree("RAM", datafile);
    
    auto* r = new RAMRecord;
    r->SetBit(quality_policy);
    
    int splitlevel = split ? 99 : 0;
    tree->Branch("RAMRecord.", &r, 64000, splitlevel);
    tree->SetMaxTreeSize(500000000000LL);
    
    if (!cache) {
        tree->SetCacheSize(0);
    }
    
    TList* headers = new TList;
    headers->SetName("headers");
    tree->GetUserInfo()->Add(headers);
    
    RAMCore::SamParser parser;
    
    auto header_callback = [&headers](const std::string& tag, const std::string& content) {
        headers->Add(new TNamed(tag.c_str(), content.c_str()));
    };
    
    auto record_callback = [&](const RAMCore::SamRecord& sam_record, size_t record_num) {
        
        r->SetQNAME(sam_record.qname.c_str());
        r->SetFLAG(sam_record.flag);
        r->SetREFID(sam_record.rname.c_str());
        r->SetPOS(sam_record.pos);
        r->SetMAPQ(sam_record.mapq);
        r->SetCIGAR(sam_record.cigar.c_str());
        r->SetREFNEXT(sam_record.rnext.c_str());
        r->SetPNEXT(sam_record.pnext);
        r->SetTLEN(sam_record.tlen);
        r->SetSEQ(sam_record.seq.c_str());
        r->SetQUAL(sam_record.qual.c_str());
        
        r->ResetNOPT();
        for (const auto& opt : sam_record.optional_fields) {
            r->SetOPT(opt.c_str());
        }
        
        tree->Fill();
        
        if (index && record_num % 1000 == 0) {
            RAMRecord::GetIndex()->AddItem(r->GetREFID(), r->GetPOS(), record_num);
        }
    };
    
    if (!parser.ParseFile(datafile, header_callback, record_callback)) {
        printf("Failed to parse SAM file %s\n", datafile);
        delete f;
        return;
    }
    
    tree->Print();
    tree->Write();
    
    RAMRecord::GetRnameRefs()->Print();
    RAMRecord::GetRnextRefs()->Print();
    RAMRecord::WriteAllRefs();
    
    if (index) {
        RAMRecord::WriteIndex();
    }
    
    delete f;
    
    printf("\nProcessed %zu SAM headers\n", parser.GetLinesProcessed() - parser.GetRecordsProcessed());
    printf("Processed %zu SAM records\n\n", parser.GetRecordsProcessed());
    
    stopwatch.Print();
}

