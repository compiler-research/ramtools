#ifndef RAMCORE_RAMSTATS_H
#define RAMCORE_RAMSTATS_H

#include <cstdint>
#include <string>
#include <map>

namespace ramcore {

/// Holds summary statistics computed from a RAM (RNTuple) file.
struct RAMStats {
   uint64_t total_reads = 0;
   uint64_t mapped_reads = 0;
   uint64_t unmapped_reads = 0;
   uint64_t duplicate_reads = 0;
   uint64_t paired_reads = 0;
   uint64_t properly_paired_reads = 0;
   uint64_t forward_strand = 0;
   uint64_t reverse_strand = 0;

   double mean_mapping_quality = 0.0;
   double mean_read_length = 0.0;

   uint64_t total_bases = 0;

   std::map<std::string, uint64_t> reads_per_chromosome;

   void Print() const;
};

/// Compute statistics from an RNTuple-based RAM file.
/// \param filename  Path to the .root file produced by samtoramntuple.
/// \param verbose   If true, print per-chromosome breakdown.
/// \return Populated RAMStats struct.
RAMStats ComputeStats(const char *filename, bool verbose = false);

} // namespace ramcore

#endif // RAMCORE_RAMSTATS_H
