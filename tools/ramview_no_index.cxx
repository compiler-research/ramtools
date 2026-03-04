//
// View a region of a RAM file, don't use index, but brute force scan.
//
// Author: Jose Javier Gonzalez Ortiz, 5/7/2017
//

#include <iostream>
#include <string>

#include <TBranch.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TTreePerfStats.h>
#include "ttree/RAMRecord.h"

void ramview_no_index(const char *file, const char *query, bool cache = false, bool perfstats = false,
                      const char *perfstatsfilename = "perf.root")
{
   TStopwatch stopwatch;
   stopwatch.Start();

   // Open the file and load tree and reader
   auto f = TFile::Open(file);
   auto t = RAMRecord::GetTree(f);

   RAMRecord *r = 0;

   if (!cache) {
      t->SetCacheSize(0);
   }

   TTreePerfStats *ps = 0;

   if (perfstats) {
      ps = new TTreePerfStats("ioperf", t);
   }

   t->SetBranchAddress("RAMRecord.", &r);

   TBranch *b = t->GetBranch("RAMRecord.");

   // Parse queried region string
   std::string region = query;
   int chrDelimiterPos = region.find(":");
   TString rname = region.substr(0, chrDelimiterPos);

   int rangeDelimiterPos = region.find("-");

   UInt_t rangeStart = std::stoi(region.substr(chrDelimiterPos + 1, rangeDelimiterPos - chrDelimiterPos));
   UInt_t rangeEnd = std::stoi(region.substr(rangeDelimiterPos + 1, region.size() - rangeDelimiterPos));

   // Default values to ensure correctness
   int rnameStart = -1;
   int posStart = -1;

   // Assume RNAME are chunked together
   // We look only at the RNAME column
   // We can only do this when there are columns
   if (b->GetSplitLevel() > 0) {
      t->SetBranchStatus("RAMRecord.*", 0);
      t->SetBranchStatus("RAMRecord.v_rname", 1);
   }

   for (int i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (rname.EqualTo(r->GetRNAME())) {
         rnameStart = i;
         break;
      }
   }

   // If the RNAME was found
   if (rnameStart >= 0) {

      // We need to look both at the leftmost position (v_pos)
      // as well as the length of sequence (v_lseq)
      if (b->GetSplitLevel() > 0) {
         t->SetBranchStatus("RAMRecord.v_pos", 1);
         t->SetBranchStatus("RAMRecord.v_lseq", 1);
      }

      for (int i = rnameStart; i < t->GetEntries(); i++) {
         t->GetEntry(i);

         // If the RNAME region ends
         if (!rname.EqualTo(r->GetRNAME())) {
            break;
         } else {
            // GetPOS() returns 1-based SAM coordinate, convert to internal 0-based
            if (r->GetPOS() - 1 + r->GetSEQLEN() > rangeStart) {
               // Register first valid position for printing
               posStart = i;
               break;
            }
         }
      }

      // If the position was found
      if (posStart >= 0) {

         // Enable all fields for printing
         if (b->GetSplitLevel() > 0) {
            t->SetBranchStatus("RAMRecord.*", 1);
         }
         for (int i = posStart; i < t->GetEntries(); i++) {
            t->GetEntry(i);

            // If the RNAME region ends
            if (!rname.EqualTo(r->GetRNAME())) {
               break;
            } else {
               // Within the region
               // GetPOS() returns 1-based SAM coordinate, convert to internal 0-based
               if (r->GetPOS() - 1 <= rangeEnd) {
                  r->Print();
               } else {
                  break;
               }
            }
         }
      }
   }

   stopwatch.Print();

   if (perfstats) {
      ps->SaveAs(perfstatsfilename);
      delete ps;
      printf("Reading %lld bytes in %d transactions\n", f->GetBytesRead(), f->GetReadCalls());
   }
}

