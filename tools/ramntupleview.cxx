#include "ramcore/RAMNTupleView.h"
#include "ramcore/RDF_RAMNTupleView.h"
#include <iostream>
#include <stdio.h>
#include <Rtypes.h>

int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <file.root> [rname:start-end]\n";
      std::cerr << "Example: " << argv[0] << " output.root chr1:1000-2000\n";
      return 1;
   }

   const char *file = argv[1];
   const char *region_str = (argc > 2) ? argv[2] : "";

   ULong64_t s = rdf_ramntupleview(file);
   Long64_t read_count = ramntupleview(file, region_str);

   printf("Found %lld records in region %s\n", read_count, region_str);

   return 0;
}
