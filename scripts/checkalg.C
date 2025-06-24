//
// Read a RAM file.
//
// Author: Fons Rademakers, 29/6/2017
//

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <iostream>
#include <Compression.h>
#include "../lib/RAMRecord.C"

void checkalg(const char *file = "ramexample.root")
{
   auto f = TFile::Open(file);
   auto t = RAMRecord::GetTree(f);

   RAMRecord *r = 0;

   t->SetBranchAddress("RAMRecord.", &r);

   // t->Print();
   std::cout << "LZMA = " << ROOT::kLZMA << std::endl;
   std::cout << "GZIP = " << ROOT::kZLIB << std::endl;
   std::cout << "COMPALG = " << f->GetCompressionAlgorithm() << std::endl;
   std::cout << "COMPLEV = " << f->GetCompressionLevel() << std::endl;

   TBranch *b = t->GetBranch("RAMRecord.");

   std::cout << "SPLIT = " << b->GetSplitLevel() << std::endl;
}

