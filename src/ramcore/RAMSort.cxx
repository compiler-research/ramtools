#include "ramcore/RAMSort.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TFile.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

int ramsortntuple(const char *inputFile, const char *outputFile, bool byName)
{
   RAMNTupleRecord::ReadAllRefs(inputFile);

   std::unique_ptr<ROOT::RNTupleReader> reader;
   try {
      reader = ROOT::RNTupleReader::Open("RAM", inputFile);
   } catch (const std::exception &e) {
      std::cerr << "Error opening input: " << e.what() << "\n";
      return 1;
   }

   const uint64_t nEntries = reader->GetNEntries();
   if (nEntries == 0) {
      std::cerr << "Input file has no entries.\n";
      return 1;
   }

   auto viewRefId = reader->GetView<int32_t>("record.refid");
   auto viewPos   = reader->GetView<int32_t>("record.pos");
   auto viewQname = reader->GetView<std::string>("record.qname");

   std::vector<uint64_t> order(nEntries);
   std::iota(order.begin(), order.end(), 0);

   std::cout << "Sorting " << nEntries << " records";
   if (byName)
      std::cout << " by QNAME...\n";
   else
      std::cout << " by coordinate (refid, pos)...\n";

   if (byName) {
      std::vector<std::string> qnames(nEntries);
      for (uint64_t i = 0; i < nEntries; ++i)
         qnames[i] = viewQname(i);
      std::stable_sort(order.begin(), order.end(),
                       [&](uint64_t a, uint64_t b) { return qnames[a] < qnames[b]; });
   } else {
      std::vector<int32_t> refids(nEntries), positions(nEntries);
      for (uint64_t i = 0; i < nEntries; ++i) {
         refids[i]    = viewRefId(i);
         positions[i] = viewPos(i);
      }
      std::stable_sort(order.begin(), order.end(), [&](uint64_t a, uint64_t b) {
         if (refids[a] != refids[b])
            return refids[a] < refids[b];
         return positions[a] < positions[b];
      });
   }

   auto viewRecord = reader->GetView<RAMNTupleRecord>("record");

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(outputFile, "RECREATE"));
   if (!rootFile || !rootFile->IsOpen()) {
      std::cerr << "Error: could not create output file " << outputFile << "\n";
      return 1;
   }

   RAMNTupleRecord::InitializeRefs();
   auto model = RAMNTupleRecord::MakeModel();
   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(505);
   auto writer    = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, writeOptions);
   auto entry     = writer->GetModel().CreateEntry();
   auto recordPtr = entry->GetPtr<RAMNTupleRecord>("record");

   for (uint64_t idx : order) {
      *recordPtr = viewRecord(idx);
      writer->Fill(*entry);
   }

   RAMNTupleRecord::WriteAllRefs(*rootFile);
   RAMNTupleRecord::WriteIndex(*rootFile);

   std::cout << "Sorted output written to " << outputFile << "\n";
   return 0;
}
