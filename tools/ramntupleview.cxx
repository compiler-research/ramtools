#include "ramcore/RAMNTupleView.h"
#include <iostream>

int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <file.root> [rname:start-end]\n";
      std::cerr << "Example: " << argv[0] << " output.root chr1:1000-2000\n";
      return 1;
   }

   const char *file = argv[1];
   const char *region = (argc > 2) ? argv[2] : "";

   ramntupleview(file, region);
   return 0;
}
