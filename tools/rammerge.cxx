//
// Merge two or more RAM files into another
//
// Author: Jose Javier Gonzalez Ortiz, 31/7/2017
//

#include <TList.h>
#include <TFile.h>
#include <TTree.h>

#include "RAMRecord.h"

void rammerge(const char *outfile, const char *infile1, const char *infile2)
{

   printf("*** Needs to be fixed to use new index mechanism. ***\n");
   return;

   TList *list = new TList;

   auto *in1 = TFile::Open(infile1);
   auto t1 = RAMRecord::GetTree(in1);
   list->Add(t1);

   auto *in2 = TFile::Open(infile2);
   auto t2 = RAMRecord::GetTree(in2);
   list->Add(t2);

   TFile *out = new TFile(outfile, "RECREATE");
   TTree *newtree = TTree::MergeTrees(list);
   newtree->SetMaxTreeSize(500000000000LL); // Default is 100GB, change to 500GB
   newtree->SetName("RAM");
   newtree->BuildIndex("v_refid", "v_pos");
   newtree->Write();

   out->Close();
   in1->Close();
   in2->Close();

   delete list;
}

