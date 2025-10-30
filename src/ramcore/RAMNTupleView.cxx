#include "ramcore/RAMNTupleView.h"
#include <iostream>
#include <cstring>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TStopwatch.h>
#include <TString.h>
#include "rntuple/RAMNTupleRecord.h"
#include <Rtypes.h> 

Long64_t ramntupleview(const char *file, const char *query, bool cache, bool perfstats, const char *perfstatsfilename)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      printf("ramntupleview: failed to open file %s\n", file);
      return 0; 
   }

   std::string region = query;
   int chrDelimiterPos = region.find(":");
   if (chrDelimiterPos == std::string::npos) {
      std::cerr << "Invalid region format. Use rname:start-end\n";
      return 0; 
   }
   TString rname = region.substr(0, chrDelimiterPos);
   int rangeDelimiterPos = region.find("-", chrDelimiterPos);
   if (rangeDelimiterPos == std::string::npos) {
      std::cerr << "Invalid region format. Use rname:start-end\n";
      return 0; 
   }
   Int_t range_start = std::stoi(region.substr(chrDelimiterPos + 1, rangeDelimiterPos - chrDelimiterPos - 1));
   Int_t range_end = std::stoi(region.substr(rangeDelimiterPos + 1));

   auto refid = RAMNTupleRecord::GetRnameRefs()->GetRefId(rname.Data());
   auto index = RAMNTupleRecord::GetIndex();

   auto recordView = reader->GetView<RAMNTupleRecord>("record");

   Long64_t count = 0; 

   if (!index || index->Size() == 0) {
      
      for (auto i : reader->GetEntryRange()) {
         const auto &rec = recordView(i);
         if (rec.refid == refid && rec.pos >= range_start - 1 && rec.pos <= range_end - 1) {
            count++;
         }
      }
   } else {
      
      auto start_entry = index->GetRow(refid, range_start);
      auto end_entry = index->GetRow(refid, range_end);
      if (start_entry < 0)
         start_entry = 0;
      if (end_entry < 0)
         end_entry = reader->GetNEntries();

      for (; start_entry < end_entry; start_entry++) {
         const auto &rec = recordView(start_entry);
         int seqlen = rec.GetSEQLEN();
         if (rec.pos + seqlen > range_start - 1) {
            break;
         }
      }

      Long64_t j;
      for (j = start_entry; j < reader->GetNEntries(); j++) {
         const auto &rec = recordView(j);
         
         if (rec.pos >= range_end) {
            break;
         }
         count++;
      }
   }

   stopwatch.Print(); 
   
   return count; 
}

