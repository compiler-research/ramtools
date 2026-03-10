/// \file ramsort.cxx
/// \brief Command-line tool to sort a RAM (RNTuple) file by coordinate or name.
///
/// Usage:
///   ./tools/ramsort <input.root> <output.root> [--by-name]

#include "ramcore/RAMSort.h"
#include <cstring>
#include <iostream>

int main(int argc, char **argv)
{
   if (argc < 3) {
      std::cerr << "Usage: " << argv[0] << " <input.root> <output.root> [--by-name]\n";
      std::cerr << "  Sort a RAM (RNTuple) file by genomic coordinate (refid, pos).\n";
      std::cerr << "  --by-name  Sort by QNAME instead.\n";
      return 1;
   }

   bool byName = false;
   for (int i = 3; i < argc; ++i) {
      if (std::strcmp(argv[i], "--by-name") == 0)
         byName = true;
   }

   return ramsortntuple(argv[1], argv[2], byName);
}
