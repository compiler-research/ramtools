//
// Class in order to access and create a random subset of a given BAM file of given size
// Implemented in order to test vs bamtools random function.
//
// Author: Taghi Aliyev, Date: 17/07/2017
//

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <iostream>

using namespace std;

#include "ttree/RAMRecord.h"

void ramrandom(const char *file = "/eos/genome/local/14007a/realigned_SAM/6148.root", const char *outfile = "out.out",
               int n = 10)
{
   auto f = TFile::Open(file);
   auto t = RAMRecord::GetTree(f);
   if (!t) {
      ::Error("ramrandom", "file %s, not found or open", file);
      return;
   }

   RAMRecord *r = 0;

   t->SetBranchAddress("RAMRecord.", &r);

   if (n > t->GetEntries()) {
      cout << "Error : n is larger than number of entries!" << endl;
      return;
   }

   cout << "There are : " << t->GetEntries() << " entries" << endl;
   UInt_t numberOfSamples = t->GetEntries();

   // Random access loop
   for (int i = 0; i < n; i++) {
      int index = gRandom->Integer(numberOfSamples);
      t->GetEvent(index);
      cout << "Accessing : " << index << endl;
      r->Print();
   }
}
