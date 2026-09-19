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

#include <cstdint>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

void samtoramntuple(const char *datafile, const char *treefile, bool split, bool cache, int compression_algorithm,
                    uint32_t quality_policy)
{
   ROOT::EnableImplicitMT();
   TStopwatch stopwatch;
   stopwatch.Start();

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (!rootFile || !rootFile->IsOpen()) {
      std::cout << "Failed to create RAM file " << treefile << "\n";
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

    auto header_callback = [&headers](const std::string& tag, const std::string& content) {
        headers.Add(new TNamed(tag.c_str(), content.c_str()));

        if (tag == "@SQ") {
            size_t sn_pos = content.find("SN:");
            if (sn_pos != std::string::npos) {
                sn_pos += 3;
                size_t tab_pos = content.find('\t', sn_pos);
                std::string ref_name =
                   content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
                RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name);
            }
        }
    };

    auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
       recordPtr->SetBit(quality_policy);

       recordPtr->SetQNAME(sam_record.qname);
       recordPtr->SetFLAG(sam_record.flag);
       recordPtr->SetREFID(sam_record.rname);
       recordPtr->SetPOS(sam_record.pos);
       recordPtr->SetMAPQ(sam_record.mapq);
       recordPtr->SetCIGAR(sam_record.cigar);
       recordPtr->SetREFNEXT(sam_record.rnext);
       recordPtr->SetPNEXT(sam_record.pnext);
       recordPtr->SetTLEN(sam_record.tlen);
       recordPtr->SetSEQ(sam_record.seq);
       recordPtr->SetQUAL(sam_record.qual);

       recordPtr->ResetNOPT();
       for (const auto &opt : sam_record.optional_fields) {
          recordPtr->SetOPT(opt);
       }

       RAMNTupleRecord::NoteRefSpan(recordPtr->GetRefSpan());
       RAMNTupleRecord::NotePlacement(recordPtr->GetREFID(), recordPtr->GetPOS() - 1);
       writer->Fill(*defaultEntry);
    };

    if (!parser.ParseFile(datafile, header_callback, record_callback)) {
        printf("Failed to parse SAM file %s\n", datafile);
        return;
    }

    writer.reset();

    // Region queries can only seek on a sorted file; the file records which it is.
    if (!RAMNTupleRecord::IsCoordinateSorted())
       fprintf(stderr, "%s is not in coordinate order; region queries will read it in full.\n", datafile);
    RAMNTupleRecord::WriteAllRefs(*rootFile);

    // One key for the list; without kSingleKey every line is written as its own
    // key and no reader can get the header back in order.
    headers.Write("headers", TObject::kSingleKey);
    rootFile->Close();

    printf("\nRAM file created: %s\n", treefile);
    printf("Number of entries: %zu\n", parser.GetRecordsProcessed());

    RAMNTupleRecord::GetRnameRefs()->Print();
    RAMNTupleRecord::GetRnextRefs()->Print();

    printf("\nProcessed %zu SAM headers\n", parser.GetLinesProcessed() - parser.GetRecordsProcessed());
    printf("Processed %zu SAM records\n\n", parser.GetRecordsProcessed());

    stopwatch.Print();
}

namespace {

// One output file, held open while the input streams past it.
struct ChromosomeWriter {
   std::unique_ptr<TFile> file{};
   std::unique_ptr<ROOT::RNTupleWriter> writer{};
   std::unique_ptr<ROOT::REntry> entry{};
   std::shared_ptr<RAMNTupleRecord> record{};
   int64_t rows = 0;
   int32_t last_pos = -1;
   bool sorted = true;
};

} // namespace

void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy)
{
   RAMNTupleRecord::InitializeRefs();

   std::map<std::string, ChromosomeWriter> writers;
   TList headers;
   headers.SetName("headers");
   headers.SetOwner(true);

   auto header_callback = [&](const std::string &tag, const std::string &content) {
      headers.Add(std::make_unique<TNamed>(tag.c_str(), content.c_str()).release());

      if (tag == "@SQ") {
         size_t sn_pos = content.find("SN:");
         if (sn_pos != std::string::npos) {
            sn_pos += 3;
            size_t tab_pos = content.find('\t', sn_pos);
            std::string ref_name =
               content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
            RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name);
         }
      }
   };

   auto open_writer = [&](const std::string &chr) -> ChromosomeWriter & {
      auto it = writers.find(chr);
      if (it != writers.end())
         return it->second;

      ChromosomeWriter &cw = writers[chr];
      std::string filename{output_prefix};
      filename += "_";
      filename += chr;
      filename += ".root";
      cw.file.reset(TFile::Open(filename.c_str(), "RECREATE"));
      if (!cw.file || !cw.file->IsOpen())
         throw std::runtime_error("cannot create " + filename);

      ROOT::RNTupleWriteOptions writeOptions;
      writeOptions.SetCompression(compression_algorithm);
      writeOptions.SetMaxUnzippedPageSize(64000);
      // Every chromosome's buffers are open at once, so keep the clusters small.
      writeOptions.SetApproxZippedClusterSize(8 * 1024 * 1024);

      cw.writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *cw.file, writeOptions);
      cw.entry = cw.writer->GetModel().CreateEntry();
      cw.record = cw.entry->GetPtr<RAMNTupleRecord>("record");
      return cw;
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
      // A record with no reference has no chromosome file to go to.
      if (sam_record.rname == "*")
         return;

      ChromosomeWriter &cw = open_writer(sam_record.rname);
      RAMNTupleRecord &rec = *cw.record;

      rec.SetBit(quality_policy);
      rec.SetQNAME(sam_record.qname);
      rec.SetFLAG(sam_record.flag);
      rec.SetREFID(sam_record.rname);
      rec.SetPOS(sam_record.pos);
      rec.SetMAPQ(sam_record.mapq);
      rec.SetCIGAR(sam_record.cigar);
      rec.SetREFNEXT(sam_record.rnext);
      rec.SetPNEXT(sam_record.pnext);
      rec.SetTLEN(sam_record.tlen);
      rec.SetSEQ(sam_record.seq);
      rec.SetQUAL(sam_record.qual);

      rec.ResetNOPT();
      for (const auto &opt : sam_record.optional_fields)
         rec.SetOPT(opt);

      RAMNTupleRecord::NoteRefSpan(rec.GetRefSpan());
      cw.writer->Fill(*cw.entry);
      cw.rows++;

      // One reference per file, so the order check is on the position alone.
      const int32_t pos = rec.GetPOS() - 1;
      if (pos < cw.last_pos)
         cw.sorted = false;
      cw.last_pos = pos;
   };

   ramcore::SamParser parser;
   if (!parser.ParseFile(datafile, header_callback, record_callback)) {
      printf("Failed to parse SAM file %s\n", datafile);
      return;
   }

   // The reference table and the longest span are only complete once the whole
   // input has been read, so every file is finished here.
   for (auto &[chr, cw] : writers) {
      cw.writer.reset();

      RAMNTupleRecord::SetCoordinateSorted(cw.sorted);
      if (!cw.sorted)
         fprintf(stderr, "%s: %s is not in coordinate order; region queries will read it in full.\n", datafile,
                 chr.c_str());
      RAMNTupleRecord::WriteAllRefs(*cw.file);

      cw.file->cd();
      headers.Write("headers", TObject::kSingleKey);
      cw.file->Close();

      printf("%s_%s.root: %lld records\n", output_prefix, chr.c_str(), static_cast<long long>(cw.rows));
   }
}
