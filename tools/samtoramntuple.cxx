#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <algorithm>
#include <cstring>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
   std::vector<std::string> args;
   args.reserve(static_cast<size_t>(argc));
   std::for_each_n(argv, static_cast<size_t>(argc), [&](char *arg) { args.emplace_back(arg); });

   if (args.size() < 2) {
      std::cout << "Usage: " << args[0] << " <input.sam> [output]\n";
      std::cout << "Options:\n";
      std::cout << "  -split             Split by chromosome\n";
      std::cout << "  -only <chr>        When used with -split, emit only the specified chromosome (repeatable)\n";
      std::cout << "  -noindex           Disable indexing\n";
      std::cout << "  -illumina          Use Illumina quality binning\n";
      std::cout << "  -dropqual          Drop quality scores\n";
      return 1;
   }

   const char *input = args[1].c_str();
   const char *output = nullptr;

   bool do_split = false;
   bool do_index = true;
   uint32_t quality_mode = RAMNTupleRecord::kPhred33;
   std::vector<std::string> only_chromosomes;

   for (int i = 2; i < static_cast<int>(args.size()); i++) {
      std::string arg = args[static_cast<size_t>(i)];
      if (arg == "-split") {
         do_split = true;
      } else if (arg == "-only") {
         if (i + 1 < static_cast<int>(args.size())) {

            only_chromosomes.emplace_back(args[static_cast<size_t>(++i)]);
         } else {
            std::cerr << "-only expects a chromosome name\n";
            return 1;
         }
      } else if (arg == "-noindex") {
         do_index = false;
      } else if (arg == "-illumina" || arg == "-dropqual") {
         quality_mode = (arg == "-illumina") ? RAMNTupleRecord::kIlluminaBinning : RAMNTupleRecord::kDrop;
      } else if (!arg.empty() && arg[0] != '-') {
         output = args[static_cast<size_t>(i)].c_str();
      }
   }

   std::string outfile{};
   if (!output) {
      outfile = std::string(input);
      size_t pos = outfile.rfind(".sam");
      if (pos != std::string::npos) {
         outfile.erase(pos);
      }
      output = outfile.c_str();
   }

   try {
      if (do_split) {
         samtoramntuple_split_by_chromosome(input, output, 505, quality_mode, 4, only_chromosomes);
      } else {
         std::string ramfile = std::string(output);
         if (ramfile.find(".root") == std::string::npos && ramfile.find(".ram") == std::string::npos) {
            ramfile += ".ram";
         }
         samtoramntuple(input, ramfile.c_str(), do_index, true, true, 505, quality_mode);
      }
   } catch (const std::exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return 1;
   }

   return 0;
}
