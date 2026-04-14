#include "ramcore/RAMNTupleView.h"
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
   bool mt_set = false;
   int num_threads = 0;
   if (argc > 3) {
      std::string mt_flag = argv[3];
      if (mt_flag.find("-m") != std::string::npos) {
         mt_set = true;
         num_threads = std::atoi(mt_flag.substr(2).data());
      }
   }

   const char *file = argv[1];
   const char *region_str = (argc > 2) ? argv[2] : "";
   ULong64_t read_count;
   if (mt_set) {
      read_count = mt_ramntupleview(num_threads, file, region_str);
   } else {
      read_count = ramntupleview(file, region_str);
   }
   printf("Found %lld records in region %s", read_count, region_str);

   return 0;
}
