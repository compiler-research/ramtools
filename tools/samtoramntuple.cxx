// samtoramntuple.C
// SAM to RAM converteR

#include "../src/rntuple/RAMNTupleRecord.cxx"
#include "../inc/rntuple/RAMNTupleRecord.h"
#include "../inc/ttree/Utils.h"
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TStopwatch.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <filesystem>

void samtoramntuple(const char *datafile = "samexample.sam",
              const char *treefile = "ramexample.root",
              bool index = true, bool split = true, bool cache = true,
              int compression_algorithm = 505,  // ZSTD level 5
              uint32_t quality_policy = RAMNTupleRecord::kPhred33)
{
    // Convert a SAM file into a RAM file using RNTuple
    
    TStopwatch stopwatch;
    stopwatch.Start();
    
    auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
    if (!rootFile || !rootFile->IsOpen()) {
        printf("Failed to create RAM file %s\n", treefile);
        return;
    }
    
    FILE *fp = fopen(datafile, "r");
    if (!fp) {
        printf("file %s not found\n", datafile);
        return;
    }
    
    RAMNTupleRecord::InitializeRefs();
    
    auto model = RAMNTupleRecord::MakeModel();
    
    ROOT::RNTupleWriteOptions writeOptions;
    writeOptions.SetCompression(compression_algorithm);
    writeOptions.SetMaxUnzippedPageSize(64000); 
    
    // Create writer
    auto writer = ROOT::Experimental::RNTupleWriter::Append(
        std::move(model), "RAM", *rootFile, writeOptions);
    auto defaultEntry = writer->CreateEntry();
    auto recordPtr = defaultEntry->GetPtr<RAMNTupleRecord>("record");
    
    // Store SAM headers
    TList headers;
    headers.SetName("headers");
    
    Long64_t nlines = 0;
    Long64_t nrecords = 0;
    
    const int maxl = 10240;
    char line[maxl];
    
    while (fgets(line, maxl, fp)) {
        bool header = false;
        int ntok = 0;
        char *tok;
        
        // Create a RAMNTupleRecord for this line
        RAMNTupleRecord r;
        r.SetBit(quality_policy);
        
        while ((tok = strtok(ntok ? 0 : line, "\t"))) {
            if ((ntok == 0 && tok[0] == '@') || header) {
                
                headers.Add(new TNamed(tok, line+strlen(tok)+1));
                header = true;
                
                if (strncmp(tok, "@SQ", 3) == 0) {
                    char *sn = strstr(line, "SN:");
                    if (sn) {
                        sn += 3;
                        char *end = strchr(sn, '\t');
                        if (end) *end = '\0';
                        RAMNTupleRecord::GetRnameRefs()->GetRefId(sn);
                        if (end) *end = '\t';
                    }
                }
                break;
            } else {
                // Parse SAM record fields
                if (ntok == 0) r.SetQNAME(tok);      // qname
                if (ntok == 1) r.SetFLAG(atoi(tok)); // flag
                if (ntok == 2) r.SetREFID(tok);      // rname
                if (ntok == 3) r.SetPOS(atoi(tok));  // pos
                if (ntok == 4) r.SetMAPQ(atoi(tok)); // mapq
                if (ntok == 5) r.SetCIGAR(tok);      // cigar
                if (ntok == 6) r.SetREFNEXT(tok);    // rnext
                if (ntok == 7) r.SetPNEXT(atoi(tok));// pnext
                if (ntok == 8) r.SetTLEN(atoi(tok)); // tlen
                if (ntok == 9) r.SetSEQ(tok);        // seq
                if (ntok == 10) {                     // qual
                    stripcrlf(tok);
                    r.SetQUAL(tok);
                    r.ResetNOPT();
                }
                if (ntok >= 11) {                     // optional fields
                    stripcrlf(tok);
                    r.SetOPT(tok);
                }
            }
            ntok++;
        }
        
        if (!header) {
            
            *recordPtr = std::move(r);
    
            writer->Fill(*defaultEntry);

            if (index && nrecords % 1000 == 0) {
                RAMNTupleRecord::GetIndex()->AddItem(recordPtr->GetREFID(), recordPtr->GetPOS() - 1, nrecords);
            }
            nrecords++;
        }
        nlines++;
    }
    
    writer.reset();
    
    if (index) {
        RAMNTupleRecord::WriteIndex(*rootFile);
    }
    RAMNTupleRecord::WriteAllRefs(*rootFile);
    
    headers.Write();
    rootFile->Close();
    
    printf("\nRAM file created: %s\n", treefile);
    printf("Number of entries: %lld\n", nrecords);
    
    RAMNTupleRecord::GetRnameRefs()->Print();
    RAMNTupleRecord::GetRnextRefs()->Print();
    
    if (index) {
        printf("\nIndex entries: %zu\n", RAMNTupleRecord::GetIndex()->Size());
    }
    
    fclose(fp);
    
    printf("\nProcessed %lld SAM headers\n", nlines-nrecords);
    printf("Processed %lld SAM records\n\n", nrecords);
    
    stopwatch.Print();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.sam> [output.ram]\n";
        std::cout << "Options:\n";
        std::cout << "  -noindex     Disable indexing\n";
        std::cout << "  -illumina    Use Illumina quality binning\n";
        std::cout << "  -dropqual    Drop quality scores\n";
        return 1;
    }
    
    const char* input = argv[1];
    const char* output = (argc > 2 && argv[2][0] != '-') ? argv[2] : nullptr;
    
    std::string outfile;
    if (!output) {
        outfile = std::string(input);
        size_t pos = outfile.rfind(".sam");
        if (pos != std::string::npos) {
            outfile.replace(pos, 4, ".ram");
        } else {
            outfile += ".ram";
        }
        output = outfile.c_str();
    }

    bool do_index = true;
    uint32_t quality_mode = RAMNTupleRecord::kPhred33;
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-noindex") {
            do_index = false;
        } else if (arg == "-illumina") {
            quality_mode = RAMNTupleRecord::kIlluminaBinning;
        } else if (arg == "-dropqual") {
            quality_mode = RAMNTupleRecord::kDrop;
        }
    }
    
    try {
        samtoramntuple(input, output, do_index, true, true, 505, quality_mode);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

