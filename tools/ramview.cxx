
#include <iostream>
#include <cstring>

#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TTreeIndex.h>
#include <TTreePerfStats.h>
#include <Rtypes.h>

#include "ttree/Utils.h"
#include "ttree/RAMRecord.h"

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root")
{
   TStopwatch stopwatch;
   stopwatch.Start();

   auto f = TFile::Open(file);
   if (!f) {
      printf("ramview: failed to open file %s\n", file);
      return 0;
   }
   auto t = RAMRecord::GetTree(f);

   RAMRecord *r = 0;

   if (!cache)
      t->SetCacheSize(0);

   TTreePerfStats *ps = 0;
   if (perfstats)
      ps = new TTreePerfStats("ioperf", t);

   t->SetBranchAddress("RAMRecord.", &r);

   TBranch *b = t->GetBranch("RAMRecord.");

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

   Int_t range_start = std::stoi(region.substr(chrDelimiterPos + 1, rangeDelimiterPos - chrDelimiterPos));
   Int_t range_end = std::stoi(region.substr(rangeDelimiterPos + 1, region.size() - rangeDelimiterPos));
   Int_t rs0 = range_start > 0 ? range_start - 1 : 0;
   Int_t re0 = range_end > 0 ? range_end - 1 : 0;

   auto refid = RAMRecord::GetRnameRefs()->GetRefId(rname);
   const UInt_t FLAG_FILTER = 0x904;

   auto start_entry = RAMRecord::GetIndex()->GetRow(refid, rs0);
   auto end_entry = RAMRecord::GetIndex()->GetRow(refid, re0);

   printf("ramview: %s:%d (%lld) - %d (%lld)\n", rname.Data(), range_start, start_entry, range_end, end_entry);

   if (b->GetSplitLevel() > 0)
      t->SetBranchStatus("RAMRecord.*", 0);

   if (b->GetSplitLevel() > 0) {
      t->SetBranchStatus("RAMRecord.v_refid", 1);
      t->SetBranchStatus("RAMRecord.v_pos", 1);
      t->SetBranchStatus("RAMRecord.v_lseq", 1);
   }

   for (; start_entry < t->GetEntries(); start_entry++) {
      t->GetEntry(start_entry);
      if (r->GetREFID() < refid)
         continue;
      if (r->GetREFID() > refid)
         break;
      if (r->GetPOS() > re0)
         break;
      if (r->GetPOS() + r->GetSEQLEN() > rs0)
         break;
   }

   if (b->GetSplitLevel() > 0)
      t->SetBranchStatus("RAMRecord.*", 1);

   Long64_t j;
   Long64_t count = 0;

   for (j = start_entry; j < t->GetEntries(); j++) {
      t->GetEntry(j);
      if (r->GetREFID() < refid)
         continue;
      if (r->GetREFID() > refid)
         break;
      if (r->GetFLAG() & FLAG_FILTER)
         continue;
      if (r->GetPOS() > re0)
         break;
      const Int_t read_start = r->GetPOS();
      if (read_start >= rs0) {
         count++;
         continue;
      }
      const Int_t read_end = r->GetSEQLEN() > 0 ? (read_start + r->GetSEQLEN() - 1) : read_start;
      if (read_end >= rs0)
         count++;
   }

   stopwatch.Print();

   std::string region_str =
      std::string(rname.Data()) + ":" + std::to_string(range_start) + "-" + std::to_string(range_end);
   std::cout << "Found " << static_cast<long long>(count) << " records in region " << region_str << std::endl;

   if (perfstats) {
      ps->SaveAs(perfstatsfilename);
      delete ps;
   }

   return count;
}

