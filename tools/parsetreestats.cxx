//
// Parse results from TTreePerfStats
//
// Author: Jose Javier Gonzalez, 02/08/2017
//

#include <TCanvas.h>
#include <TImage.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreePerfStats.h>
#include <iostream>
#include "ttree/RAMRecord.h"

void parsetreestats(const char *file = "treestats.root", bool saveimage = false, const char *imagefile = "canvas.png")
{
   auto f = TFile::Open(file);
   auto ioperf = (TTreePerfStats *)f->Get("ioperf");
   ioperf->Print();

   // Save the image

   if (saveimage) {
      TCanvas *c = new TCanvas;
      TImage *img = TImage::Create();
      ioperf->Draw();
      img->FromPad(c);
      img->WriteImage(imagefile);
      delete img;
      delete c;
   } else {
      ioperf->Draw();
   }
}

