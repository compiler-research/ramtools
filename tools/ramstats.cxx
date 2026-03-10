#include "ramcore/RAMStats.h"
#include <iostream>
#include <cstring>

int main(int argc, char **argv)
{
   if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " <input.root> [--verbose]\n";
      std::cerr << "  Compute summary statistics for a RAM (RNTuple) file.\n";
      std::cerr << "  --verbose  Also print per-chromosome read counts.\n";
      return 1;
   }

   const char *inputFile = argv[1];
   bool verbose = false;
   for (int i = 2; i < argc; ++i) {
      if (std::strcmp(argv[i], "--verbose") == 0)
         verbose = true;
   }

   auto result = ramcore::ComputeStats(inputFile);

   if (!result.ok) {
      std::cerr << "Error: " << result.error_message << "\n";
      return 1;
   }

   // Computation and I/O are separate — tool controls printing
   result.stats.Print();

   // Per-chromosome breakdown only with --verbose
   if (verbose && !result.stats.reads_per_chromosome.empty()) {
      std::cout << "\n--- Reads per Chromosome ---\n";
      for (const auto &[chrom, count] : result.stats.reads_per_chromosome) {
         double pct = 100.0 * count / result.stats.total_reads;
         std::cout << "  " << chrom << ": " << count
                   << "  (" << pct << "%)\n";
      }
      std::cout << "\n";
   }

   return 0;
}
