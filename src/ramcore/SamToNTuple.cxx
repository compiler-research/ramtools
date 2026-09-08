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

namespace {

constexpr uint16_t kUnmapped = 0x4;
constexpr int32_t kPositionInterval = 10000;
constexpr int64_t kMappedInterval = 100;

} // namespace

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
    int32_t last_refid = -1;
    int32_t last_indexed_pos = -kPositionInterval;

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

    auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t record_num) {
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
       if (!(sam_record.flag & kUnmapped) && recordPtr->GetREFID() >= 0)
          RAMNTupleRecord::NotePlacement(recordPtr->GetREFID(), recordPtr->GetPOS() - 1);
       writer->Fill(*defaultEntry);

       // Index building: create a sparse lookup table so region queries can jump
       // directly to the relevant records instead of scanning from the beginning.
       //
       // Only mapped reads with a valid chromosome get indexed. The index maps
       // (chromosome, position) → record_number.
       //
       // An entry is created when any of these triggers fire:
       //   1. New chromosome — so queries on that chromosome have a starting point
       //   2. Position gap >= 10kb — limits how far a query must scan forward
       //   3. Every 100th mapped read — density guarantee for pileup regions
       //      where reads are close together and the 10kb gap never triggers
       //
       // Duplicate entries at the same (chromosome, position) are skipped to ensure
       // the index always points to the first record at any given position.
       if (index && !(sam_record.flag & kUnmapped) && recordPtr->GetREFID() >= 0) {
          int32_t current_refid = recordPtr->GetREFID();
          int32_t current_pos = recordPtr->GetPOS() - 1;

          bool new_chrom = (current_refid != last_refid);
          bool far_enough = (current_pos - last_indexed_pos >= kPositionInterval);
          bool periodic = (mapped_count % kMappedInterval == 0);
          bool duplicate = (!new_chrom && current_pos == last_indexed_pos);

          if ((new_chrom || far_enough || periodic) && !duplicate) {
             RAMNTupleRecord::GetIndex()->AddItem(current_refid, current_pos, record_num);
             last_refid = current_refid;
             last_indexed_pos = current_pos;
          }
          mapped_count++;
       }
    };

    if (!parser.ParseFile(datafile, header_callback, record_callback)) {
        printf("Failed to parse SAM file %s\n", datafile);
        return;
    }

    writer.reset();

    // An index is only usable on a sorted file; the file records which it is.
    const bool sorted = RAMNTupleRecord::IsCoordinateSorted();
    if (index && !sorted) {
       fprintf(stderr, "%s is not in coordinate order, so no index was written; region queries will read it in full.\n",
               datafile);
    }
    if (index && sorted) {
       RAMNTupleRecord::WriteIndex(*rootFile);
    }
    RAMNTupleRecord::WriteAllRefs(*rootFile);

    // One key for the list; without kSingleKey every line is written as its own
    // key and no reader can get the header back in order.
    headers.Write("headers", TObject::kSingleKey);
    rootFile->Close();

    printf("\nRAM file created: %s\n", treefile);
    printf("Number of entries: %zu\n", parser.GetRecordsProcessed());

    RAMNTupleRecord::GetRnameRefs()->Print();
    RAMNTupleRecord::GetRnextRefs()->Print();

    if (index && sorted) {
       printf("\nIndex entries: %zu\n", RAMNTupleRecord::GetIndex()->Size());
    }

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
   RAMNTupleIndex index{};
   int64_t rows = 0;
   int64_t mapped = 0;
   int32_t last_pos = -1;
   int32_t last_indexed_pos = -kPositionInterval;
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
      const int64_t row = cw.rows++;

      if (sam_record.flag & kUnmapped)
         return;

      const int32_t pos = rec.GetPOS() - 1;
      if (pos < cw.last_pos)
         cw.sorted = false;
      cw.last_pos = pos;

      // Same sparse-index rule as the single-file writer, kept per file so the
      // rows it records are the rows of this file.
      const bool far_enough = (pos - cw.last_indexed_pos >= kPositionInterval);
      const bool periodic = (cw.mapped % kMappedInterval == 0);
      if ((far_enough || periodic) && pos != cw.last_indexed_pos) {
         cw.index.AddItem(rec.GetREFID(), pos, row);
         cw.last_indexed_pos = pos;
      }
      cw.mapped++;
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
      if (cw.sorted)
         RAMNTupleRecord::WriteIndex(*cw.file, cw.index);
      else
         fprintf(stderr, "%s: %s is not in coordinate order, so no index was written.\n", datafile, chr.c_str());
      RAMNTupleRecord::WriteAllRefs(*cw.file);

      cw.file->cd();
      headers.Write("headers", TObject::kSingleKey);
      cw.file->Close();

      printf("%s_%s.root: %lld records, %zu index entries\n", output_prefix, chr.c_str(),
             static_cast<long long>(cw.rows), cw.sorted ? cw.index.Size() : 0);
   }
}
