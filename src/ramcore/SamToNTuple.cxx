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
#include <TROOT.h>

#include <map>
#include <memory>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <thread>
#include <vector>
#include <mutex>

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

void samtoramntuple_split_by_chromosome(const char *datafile, 
                                        const char *output_prefix,
                                        int compression_algorithm,
                                        uint32_t quality_policy)
{
    RAMNTupleRecord::InitializeRefs();
    
    std::map<std::string, std::vector<ramcore::SamRecord>> chromosome_records;
    std::vector<std::pair<std::string, std::string>> headers;
    
    ramcore::SamParser parser;
    
    auto header_callback = [&](const std::string& tag, const std::string& content) {
        headers.push_back({tag, content});
        
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
        if (sam_record.rname != "*") {
            chromosome_records[sam_record.rname].push_back(sam_record);
        }
    };
    
    parser.ParseFile(datafile, header_callback, record_callback);
    
    for (auto& [chr, records] : chromosome_records) {
        std::sort(records.begin(), records.end(), 
                  [](const ramcore::SamRecord& a, const ramcore::SamRecord& b) {
                      return a.pos < b.pos;
                  });
    }
    
    for (const auto& [chr, records] : chromosome_records) {
        std::string filename = std::string(output_prefix) + "_" + chr + ".root";
        auto file = TFile::Open(filename.c_str(), "RECREATE");
        
        auto model = RAMNTupleRecord::MakeModel();
        ROOT::Experimental::RNTupleWriteOptions writeOptions;
        writeOptions.SetCompression(compression_algorithm);
        writeOptions.SetMaxUnzippedPageSize(128000);
        
        auto writer = ROOT::Experimental::RNTupleWriter::Append(
            std::move(model), "RAM", *file, writeOptions);
        auto entry = writer->CreateEntry();
        auto recordPtr = entry->GetPtr<RAMNTupleRecord>("record").get();
        
        for (const auto& sam_record : records) {
            recordPtr->SetBit(quality_policy);
            recordPtr->SetQNAME(sam_record.qname.c_str());
            recordPtr->SetFLAG(sam_record.flag);
            recordPtr->SetREFID(sam_record.rname.c_str());
            recordPtr->SetPOS(sam_record.pos);
            recordPtr->SetMAPQ(sam_record.mapq);
            recordPtr->SetCIGAR(sam_record.cigar.c_str());
            recordPtr->SetREFNEXT(sam_record.rnext.c_str());
            recordPtr->SetPNEXT(sam_record.pnext);
            recordPtr->SetTLEN(sam_record.tlen);
            recordPtr->SetSEQ(sam_record.seq.c_str());
            recordPtr->SetQUAL(sam_record.qual.c_str());
            
            recordPtr->ResetNOPT();
            for (const auto& opt : sam_record.optional_fields) {
                recordPtr->SetOPT(opt.c_str());
            }
            
            writer->Fill(*entry);
        }
        
        writer.reset();
        
        TList h;
        h.SetName("headers");
        for (const auto& [tag, content] : headers) {
            h.Add(new TNamed(tag.c_str(), content.c_str()));
        }
        
        RAMNTupleRecord::WriteAllRefs(*file);
        h.Write();
        file->Close();
        delete file;
    }
}

void samtoramntuple_split_by_chromosome_parallel(const char *datafile, 
                                                  const char *output_prefix,
                                                  int compression_algorithm,
                                                  uint32_t quality_policy,
                                                  int num_threads)
{
    ROOT::EnableThreadSafety();
    RAMNTupleRecord::InitializeRefs();
    
    std::map<std::string, std::vector<ramcore::SamRecord>> chromosome_records;
    std::vector<std::pair<std::string, std::string>> headers;
    
    ramcore::SamParser parser;
    
    auto header_callback = [&](const std::string& tag, const std::string& content) {
        headers.push_back({tag, content});
        
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
        if (sam_record.rname != "*") {
            chromosome_records[sam_record.rname].push_back(sam_record);
        }
    };
    
    parser.ParseFile(datafile, header_callback, record_callback);
    
    for (auto& [chr, records] : chromosome_records) {
        std::sort(records.begin(), records.end(), 
                  [](const ramcore::SamRecord& a, const ramcore::SamRecord& b) {
                      return a.pos < b.pos;
                  });
    }
    
    std::vector<std::string> chr_names;
    for (const auto& [chr, records] : chromosome_records) {
        chr_names.push_back(chr);
    }
    
    static std::mutex file_destruction_mutex;
    
    auto write_chromosome = [&](const std::string& chr) {
        const auto& records = chromosome_records[chr];
        
        std::string filename = std::string(output_prefix) + "_" + chr + ".root";
        
        std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "RECREATE"));
        
        auto model = RAMNTupleRecord::MakeModel();
        ROOT::Experimental::RNTupleWriteOptions writeOptions;
        writeOptions.SetCompression(compression_algorithm);
        writeOptions.SetMaxUnzippedPageSize(128000);
        
        auto writer = ROOT::Experimental::RNTupleWriter::Append(
            std::move(model), "RAM", *file, writeOptions);
        auto entry = writer->CreateEntry();
        auto recordPtr = entry->GetPtr<RAMNTupleRecord>("record").get();
        
        for (const auto& sam_record : records) {
            recordPtr->SetBit(quality_policy);
            recordPtr->SetQNAME(sam_record.qname.c_str());
            recordPtr->SetFLAG(sam_record.flag);
            recordPtr->SetREFID(sam_record.rname.c_str());
            recordPtr->SetPOS(sam_record.pos);
            recordPtr->SetMAPQ(sam_record.mapq);
            recordPtr->SetCIGAR(sam_record.cigar.c_str());
            recordPtr->SetREFNEXT(sam_record.rnext.c_str());
            recordPtr->SetPNEXT(sam_record.pnext);
            recordPtr->SetTLEN(sam_record.tlen);
            recordPtr->SetSEQ(sam_record.seq.c_str());
            recordPtr->SetQUAL(sam_record.qual.c_str());
            
            recordPtr->ResetNOPT();
            for (const auto& opt : sam_record.optional_fields) {
                recordPtr->SetOPT(opt.c_str());
            }
            
            writer->Fill(*entry);
        }
        
        writer.reset();
        
        TList h;
        h.SetName("headers");
        for (const auto& [tag, content] : headers) {
            h.Add(new TNamed(tag.c_str(), content.c_str()));
        }
        
        RAMNTupleRecord::WriteAllRefs(*file);
        h.Write();
        
        {
            std::lock_guard<std::mutex> lock(file_destruction_mutex);
            file->Close();
            file.reset();  
        }
    };
    
    size_t chr_idx = 0;
    while (chr_idx < chr_names.size()) {
        std::vector<std::thread> threads;
        
        for (int i = 0; i < num_threads && chr_idx < chr_names.size(); ++i, ++chr_idx) {
            threads.emplace_back(write_chromosome, chr_names[chr_idx]);
        }
        
        for (auto& t : threads) {
            t.join();
        }
    }
}

