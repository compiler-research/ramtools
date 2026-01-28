//
// Read a RAM file.
//
// Author: Fons Rademakers, 29/6/2017
//

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>

#include "ttree/RAMRecord.h"

void ramreader(const char *file = "ramexample.root")
{
   auto f = TFile::Open(file);
   if (!f) {
      printf("ramreader: failed to open file %s\n", file);
      return;
   }
   auto t = RAMRecord::GetTree(f);

   RAMRecord *r = 0;

   t->SetBranchAddress("RAMRecord.", &r);

   printf("The file contains %lld RAMRecords\n\n", t->GetEntries());

   t->SetBranchStatus("RAMRecord.*", 0);
   t->SetBranchStatus("RAMRecord.TObject.fBits", 1);
   t->SetBranchStatus("RAMRecord.v_qname", 1);
   t->SetBranchStatus("RAMRecord.v_lseq", 1);
   t->SetBranchStatus("RAMRecord.v_seq", 1);
   t->SetBranchStatus("RAMRecord.v_qual", 1);

   // access sequentially first 10 records
   printf("Sequentially access the first 10 records from the file:\n");
   for (int i = 0; i < 10; i++) {
      t->GetEvent(i);
      printf("%2d QNAME: %s\n", i, r->GetQNAME());
      printf("%2d SEQ:   %s\n", i, r->GetSEQ());
      printf("%2d QUAL:  %s\n", i, r->GetQUAL());
   }

   printf("\nFull print of RAMRecord 10:\n");
   r->Print();

   // no need anymore for QNAME
   t->SetBranchStatus("RAMRecord.v_qname", 0);

   // Randomly access 10 records
   printf("\nRandomly access 10 records from the file:\n");
   for (int i = 0; i < 10; i++) {
      int n = gRandom->Rndm() * 1000.;
      t->GetEvent(n);
      printf("%2d SEQ:  %s\n", n, r->GetSEQ());
      printf("%2d QUAL: %s\n", n, r->GetQUAL());
   }

   // Get full last RAMRecord, turn on all branches
   t->SetBranchStatus("RAMRecord.*", 1);
   t->GetEvent(t->GetEntries() - 1);

   RAMRecord r2 = *r;
   printf("\nFull print of copied last RAMRecord:\n");
   r2.Print();

   printf("\nPrint RnameRefs:\n");
   RAMRecord::GetRnameRefs()->Print();
   printf("\nPrint RnextRefs:\n");
   RAMRecord::GetRnextRefs()->Print();
   printf("\nPrint Index:\n");
   RAMRecord::GetIndex()->Print();
}
