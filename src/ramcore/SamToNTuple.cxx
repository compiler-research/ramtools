#include "ramcore/SamToNTuple.h"
#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <ROOT/RNTupleFillContext.hxx>
#include <TStopwatch.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>
#include <TROOT.h>

void samtoramntuple(const char *datafile, const char *treefile, bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy) // NOLINT(misc-use-internal-linkage) - public API function
{
   (void)split;
   (void)cache;
   TStopwatch stopwatch;
   stopwatch.Start();

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (!rootFile || !rootFile->IsOpen()) {
      std::cerr << "Failed to create RAM file " << treefile << std::endl;
      return;
   }

   RAMNTupleRecord::InitializeRefs();

   auto model = RAMNTupleRecord::MakeModel();

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(compression_algorithm);
   writeOptions.SetMaxUnzippedPageSize(64000);

   auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, writeOptions);
   auto defaultEntry = writer->GetModel().CreateEntry();
   auto recordPtr = defaultEntry->GetPtr<RAMNTupleRecord>("record");

   TList headers;
   headers.SetName("headers");

   ramcore::SamParser parser;

   int64_t mapped_count = 0;
   int32_t last_refid = -2;
   int32_t last_indexed_pos = -1000000;
   const int32_t POS_INTERVAL = 10000;
   const int64_t MAPPED_INTERVAL = 100;

   auto header_callback = [&headers](const std::string &tag, const std::string &content) {
      headers.Add(new TNamed(tag.c_str(), content.c_str()));

      if (tag == "@SQ") {
         size_t sn_pos = content.find("SN:");
         if (sn_pos != std::string::npos) {
            sn_pos += 3;
            size_t tab_pos = content.find('\t', sn_pos);
            const std::string ref_name =
               content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
            RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name.c_str());
         }
      }
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t record_num) {
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
      for (const auto &opt : sam_record.optional_fields) {
         recordPtr->SetOPT(opt.c_str());
      }

      writer->Fill(*defaultEntry);

      if (index && !(sam_record.flag & 0x4) && recordPtr->GetREFID() >= 0) {
         mapped_count++;
         bool should_index = false;

         int current_refid = recordPtr->GetREFID();
         int current_pos = recordPtr->GetPOS() - 1;

         if (current_refid != last_refid) {
            should_index = true;
            last_refid = current_refid;
            last_indexed_pos = current_pos;
         }

         else if (current_pos - last_indexed_pos >= POS_INTERVAL) {
            should_index = true;
            last_indexed_pos = current_pos;
         }

         else if (mapped_count % MAPPED_INTERVAL == 0) {
            should_index = true;
         }

         if (should_index) {
            RAMNTupleRecord::GetIndex()->AddItem(current_refid, current_pos, record_num);
         }
      }
   };

   if (!parser.ParseFile(datafile, header_callback, record_callback)) {
      std::cerr << "Failed to parse SAM file " << datafile << std::endl;
      return;
   }

   writer.reset();

   if (index) {
      RAMNTupleRecord::WriteIndex(*rootFile);
   }
   RAMNTupleRecord::WriteAllRefs(*rootFile);

   headers.Write();
   rootFile->Close();

   (void)mapped_count;
   (void)stopwatch;
}
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy, int num_threads,
                                        const std::vector<std::string> &only_chromosomes)
{
   ROOT::EnableThreadSafety();
   RAMNTupleRecord::InitializeRefs();

   std::map<std::string, std::vector<ramcore::SamRecord>> chromosome_records;
   std::vector<std::pair<std::string, std::string>> headers;

   ramcore::SamParser parser;

   std::unordered_set<std::string> chromosome_filter(only_chromosomes.begin(), only_chromosomes.end());
   const bool filter_active = !chromosome_filter.empty();

   auto header_callback = [&](const std::string &tag, const std::string &content) {
      headers.push_back({tag, content});

      if (tag == "@SQ") {
         size_t sn_pos = content.find("SN:");
         if (sn_pos != std::string::npos) {
            sn_pos += 3;
            size_t tab_pos = content.find('\t', sn_pos);
            const std::string ref_name =
               content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
            RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name.c_str());
         }
      }
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t record_num) {
      if (sam_record.rname != "*") {
         if (filter_active && chromosome_filter.find(sam_record.rname) == chromosome_filter.end()) {
            return;
         }
         chromosome_records[sam_record.rname].push_back(sam_record);
      }
   };

   parser.ParseFile(datafile, header_callback, record_callback);

   std::vector<std::string> chr_names;
   for (const auto &[chr, records] : chromosome_records) {
      chr_names.push_back(chr);
   }

   std::vector<std::thread> sort_threads;
   for (const auto &chr : chr_names) {
      sort_threads.emplace_back([&chromosome_records, chr]() {
         auto &records = chromosome_records[chr];
         std::sort(records.begin(), records.end(),
                   [](const ramcore::SamRecord &a, const ramcore::SamRecord &b) { return a.pos < b.pos; });
      });
   }
   for (auto &t : sort_threads) {
      t.join();
   }

   std::mutex global_record_mutex;

   auto write_chromosome_parallel = [&](const std::string &chr, std::mutex &record_mutex) {
      const auto &records = chromosome_records[chr];

      std::string filename = std::string(output_prefix) + "_" + chr + ".root";

      auto model = RAMNTupleRecord::MakeModel();
      ROOT::RNTupleWriteOptions writeOptions;

      writeOptions.SetCompression(ROOT::RCompressionSetting::EAlgorithm::kZSTD, 1);
      writeOptions.SetApproxZippedClusterSize(200 * 1024 * 1024);
      writeOptions.SetMaxUnzippedClusterSize(1024 * 1024 * 1024);
      writeOptions.SetMaxUnzippedPageSize(1024 * 1024);
      writeOptions.SetUseBufferedWrite(true);

      auto parallel_writer =
         ROOT::Experimental::RNTupleParallelWriter::Recreate(std::move(model), "RAM", filename, writeOptions);

      const int contexts_per_file = std::min(4, num_threads);
      const size_t records_per_context = (records.size() + contexts_per_file - 1) / contexts_per_file;

      std::vector<std::thread> write_threads;

      for (int ctx = 0; ctx < contexts_per_file && ctx * records_per_context < records.size(); ++ctx) {
         write_threads.emplace_back([&, ctx]() {
            auto fill_context = parallel_writer->CreateFillContext();
            auto entry = fill_context->GetModel().CreateEntry();
            auto recordPtr = entry->GetPtr<RAMNTupleRecord>("record").get();

            size_t start = ctx * records_per_context;
            size_t end = std::min(start + records_per_context, records.size());

            for (size_t i = start; i < end; ++i) {
               const auto &sam_record = records[i];

               recordPtr->SetBit(quality_policy);
               recordPtr->SetQNAME(sam_record.qname.c_str());
               recordPtr->SetFLAG(sam_record.flag);

               {
                  std::lock_guard<std::mutex> lock(record_mutex);
                  recordPtr->SetREFID(sam_record.rname.c_str());
                  recordPtr->SetREFNEXT(sam_record.rnext.c_str());
               }

               recordPtr->SetPOS(sam_record.pos);
               recordPtr->SetMAPQ(sam_record.mapq);
               recordPtr->SetCIGAR(sam_record.cigar.c_str());
               recordPtr->SetPNEXT(sam_record.pnext);
               recordPtr->SetTLEN(sam_record.tlen);
               recordPtr->SetSEQ(sam_record.seq.c_str());
               recordPtr->SetQUAL(sam_record.qual.c_str());

               recordPtr->ResetNOPT();
               for (const auto &opt : sam_record.optional_fields) {
                  recordPtr->SetOPT(opt.c_str());
               }

               fill_context->Fill(*entry);
            }
         });
      }

      for (auto &t : write_threads) {
         t.join();
      }

      parallel_writer.reset();

      std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "UPDATE"));

      TList h;
      h.SetName("headers");
      for (const auto &[tag, content] : headers) {
         h.Add(new TNamed(tag.c_str(), content.c_str()));
      }

      RAMNTupleRecord::WriteAllRefs(*file);
      h.Write();

      file->Close();
   };

   size_t chr_idx = 0;
   while (chr_idx < chr_names.size()) {
      std::vector<std::thread> threads;

      for (int i = 0; i < num_threads && chr_idx < chr_names.size(); ++i, ++chr_idx) {
         threads.emplace_back(write_chromosome_parallel, chr_names[chr_idx], std::ref(global_record_mutex));
      }

      for (auto &t : threads) {
         t.join();
      }
   }
}

