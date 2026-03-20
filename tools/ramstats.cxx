#include "ramcore/RAMStats.h"
#include <cstring>
#include <iostream>

int main(int argc, char **argv)
{
   if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <input.root> [--verbose]\n";
      std::cerr << "  Compute summary statistics for a RAM (RNTuple) file.\n";
      std::cerr << "  --verbose  Also print per-chromosome read counts.\n";
      return 1;
   }
   bool verbose = false;
   for (int i = 2; i < argc; ++i) {
      if (std::strcmp(argv[i], "--verbose") == 0)
         verbose = true;
   }
   return ramcore::RunRamStats(argv[1], verbose);
}
