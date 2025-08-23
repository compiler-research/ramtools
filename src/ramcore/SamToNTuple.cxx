#include "ramcore/SamToNTuple.h"
#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TStopwatch.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>

void samtoramntuple(const char *datafile,
                    const char *treefile,
                    bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy)
{
    TStopwatch stopwatch;
    stopwatch.Start();
    
    auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
    if (!rootFile || !rootFile->IsOpen()) {
        printf("Failed to create RAM file %s\n", treefile);
        return;
    }
    
    RAMNTupleRecord::InitializeRefs();
    
    auto model = RAMNTupleRecord::MakeModel();
    
    ROOT::Experimental::RNTupleWriteOptions writeOptions;
    writeOptions.SetCompression(compression_algorithm);
    writeOptions.SetMaxUnzippedPageSize(64000);
    
    auto writer = ROOT::Experimental::RNTupleWriter::Append(
        std::move(model), "RAM", *rootFile, writeOptions);
    auto defaultEntry = writer->CreateEntry();
    auto recordPtr = defaultEntry->GetPtr<RAMNTupleRecord>("record");
    
    TList headers;
    headers.SetName("headers");
    
    ramcore::SamParser parser;
    
    auto header_callback = [&headers](const std::string& tag, const std::string& content) {
        headers.Add(new TNamed(tag.c_str(), content.c_str()));
        
        if (tag == "@SQ") {
            size_t sn_pos = content.find("SN:");
            if (sn_pos != std::string::npos) {
                sn_pos += 3;
                size_t tab_pos = content.find('\t', sn_pos);
                std::string ref_name = content.substr(sn_pos, 
                    tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
                RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name.c_str());
            }
        }
    };
    
    auto record_callback = [&](const ramcore::SamRecord& sam_record, size_t record_num) {
        RAMNTupleRecord r;
        r.SetBit(quality_policy);
        
        r.SetQNAME(sam_record.qname.c_str());
        r.SetFLAG(sam_record.flag);
        r.SetREFID(sam_record.rname.c_str());
        r.SetPOS(sam_record.pos);
        r.SetMAPQ(sam_record.mapq);
        r.SetCIGAR(sam_record.cigar.c_str());
        r.SetREFNEXT(sam_record.rnext.c_str());
        r.SetPNEXT(sam_record.pnext);
        r.SetTLEN(sam_record.tlen);
        r.SetSEQ(sam_record.seq.c_str());
        r.SetQUAL(sam_record.qual.c_str());
        
        r.ResetNOPT();
        for (const auto& opt : sam_record.optional_fields) {
            r.SetOPT(opt.c_str());
        }
        
        *recordPtr = std::move(r);
        writer->Fill(*defaultEntry);
        
        if (index && record_num % 1000 == 0) {
            RAMNTupleRecord::GetIndex()->AddItem(recordPtr->GetREFID(), 
                                                 recordPtr->GetPOS() - 1, record_num);
        }
    };
    
    if (!parser.ParseFile(datafile, header_callback, record_callback)) {
        printf("Failed to parse SAM file %s\n", datafile);
        return;
    }
    
    writer.reset();
    
    if (index) {
        RAMNTupleRecord::WriteIndex(*rootFile);
    }
    RAMNTupleRecord::WriteAllRefs(*rootFile);
    
    headers.Write();
    rootFile->Close();
    
    printf("\nRAM file created: %s\n", treefile);
    printf("Number of entries: %zu\n", parser.GetRecordsProcessed());
    
    RAMNTupleRecord::GetRnameRefs()->Print();
    RAMNTupleRecord::GetRnextRefs()->Print();
    
    if (index) {
        printf("\nIndex entries: %zu\n", RAMNTupleRecord::GetIndex()->Size());
    }
    
    printf("\nProcessed %zu SAM headers\n", 
           parser.GetLinesProcessed() - parser.GetRecordsProcessed());
    printf("Processed %zu SAM records\n\n", parser.GetRecordsProcessed());
    
    stopwatch.Print();
}

