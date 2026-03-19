#include "ramcore/BamtoNTuple.h"

#include "rntuple/RAMNTupleRecord.h"

#include <cstdint>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cout << "Usage: " << argv[0] << " <input.bam> [output]\n"
                << "Options:\n"
                << "  -noindex     Disable indexing\n"
                << "  -illumina    Use Illumina quality binning\n"
                << "  -dropqual    Drop quality scores\n";
      return 1;
   }

   const char *input = argv[1];
   const char *output = nullptr;

   bool do_index = true;
   uint32_t quality_mode = RAMNTupleRecord::kPhred33;

   for (int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "-noindex")
         do_index = false;
      else if (arg == "-illumina")
         quality_mode = RAMNTupleRecord::kIlluminaBinning;
      else if (arg == "-dropqual")
         quality_mode = RAMNTupleRecord::kDrop;
      else if (arg[0] != '-')
         output = argv[i];
   }

   std::string outfile;
   if (output == nullptr) {
      outfile = input;
      const auto pos = outfile.rfind(".bam");
      if (pos != std::string::npos)
         outfile.erase(pos);
      output = outfile.c_str();
   }

   std::string ramfile = output;
   if (ramfile.find(".root") == std::string::npos && ramfile.find(".ram") == std::string::npos)
      ramfile += ".ram";

   bamtoramntuple(input, ramfile.c_str(),
                  /*index=*/do_index, /*split=*/false, /*cache=*/true,
                  /*compression_algorithm=*/505, /*quality_policy=*/quality_mode);

   return 0;
}