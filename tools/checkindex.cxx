#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TTreeIndex.h>
#include <TFile.h>
#include <TString.h>
#include <cstring>
#include "ttree/Utils.h"
#include "RAMRecord.h"

void checkindex(const char *file = "ramexample.root")
{
   auto f = TFile::Open(file);
   auto t = RAMRecord::GetTree(f);

   RAMRecord *r = 0;
   t->SetBranchAddress("RAMRecord.", &r);

   // Check for TTreeIndex in TTree, same file and external file (in that order)

   auto i = t->GetTreeIndex();
   if (i) {
      std::cout << "Index present in TTree" << std::endl;
   } else {
      std::string indexfile = file;
      indexfile.append(".rai");
      std::FILE *fp = std::fopen(indexfile.c_str(), "r");
      if (fp) {
         auto f2 = TFile::Open(indexfile.c_str());
         i = (TTreeIndex *)f2->Get("INDEX");
         if (i) {
            std::cout << "HOLA" << std::endl;
            std::cout << "Index present in separate TFile" << std::endl;
         }
         delete fp;
         delete f2;
      }
   }

   if (i) {
      std::cout << "Majorname is " << i->GetMajorName() << std::endl;
      std::cout << "Minorname is " << i->GetMinorName() << std::endl;
   } else {
      std::cout << "Index Missing" << std::endl;
   }

   delete i;
   delete t;
   delete f;
}

