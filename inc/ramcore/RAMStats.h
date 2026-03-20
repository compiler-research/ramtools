#ifndef RAMCORE_RAMSTATS_H
#define RAMCORE_RAMSTATS_H

#include <cstdint>
#include <string>
#include <map>

namespace ramcore {

struct RAMStats {
   uint64_t total_reads         = 0;
   uint64_t mapped_reads        = 0;
   uint64_t unmapped_reads      = 0;
   uint64_t duplicate_reads     = 0;
   uint64_t paired_reads        = 0;
   uint64_t properly_paired_reads = 0;
   uint64_t forward_strand      = 0;
   uint64_t reverse_strand      = 0;

   double mean_mapping_quality  = 0.0;
   double mean_read_length      = 0.0;
   uint64_t total_bases         = 0;

   std::map<std::string, uint64_t> reads_per_chromosome;

   void Print() const;
};

/// Result wrapper so callers can distinguish errors from empty files
struct RAMStatsResult {
   RAMStats stats;
   bool     ok    = false;
   std::string error_message;
};

/// Compute statistics from an RNTuple-based RAM file.
/// Returns RAMStatsResult — check .ok before using .stats.
RAMStatsResult ComputeStats(const char *filename);

/// Run the full ramstats pipeline: compute stats and print results.
/// \param filename  Path to input .root RAM file
/// \param verbose   If true, print per-chromosome breakdown
/// \return 0 on success, 1 on error
int RunRamStats(const char *filename, bool verbose = false);

} // namespace ramcore

#endif // RAMCORE_RAMSTATS_H
