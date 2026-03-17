#include "ramcore/RAMNTupleView.h"
#include "ramcore/RDF_RAMNTupleView.h"
#include <iostream>
#include <stdio.h>
#include <Rtypes.h>
//#define RDF_IMPL
int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <file.root> [rname:start-end]\n";
      std::cerr << "Example: " << argv[0] << " output.root chr1:1000-2000\n";
      return 1;
   }

   const char *file = argv[1];
   const char *region_str = (argc > 2) ? argv[2] : "";
//#ifdef RDF_IMPL
   ULong64_t s = rdf_ramntupleview(file, region_str);
//#else
 Long64_t read_count = ramntupleview(file, region_str);
//#endif
   printf("Found %lld records in region %s [single thread]\n Found %lld records in region %s [multi-thread]", read_count, region_str, s, region_str);

   return 0;
}
