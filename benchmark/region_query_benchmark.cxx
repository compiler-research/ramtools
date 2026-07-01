#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include <Rtypes.h>
#include <numeric>
#include <string>
#include <vector>

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root");
Long64_t ramntupleview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");

// A scale ladder of regions, from single-position to whole-chromosome. Indices into
// this list are what --regions selects (default 0,3,6,9 keeps the historical set).
static const std::vector<std::string> kRegions = {"chr1:1000000-1001000",
                                                  "chr2:5000000-5010000",
                                                  "chrX:100000-150000",
                                                  "chr1:1000000-2000000",
                                                  "chr5:10000000-15000000",
                                                  "chr10:50000000-60000000",
                                                  "chr1:1-50000000",
                                                  "chr2:1-100000000",
                                                  "chr7:50000000-150000000",
                                                  "chr21:1-48129895",
                                                  "chrM:1-16571",
                                                  "chrY:2600000-2700000",
                                                  "GL000227.1:1-100000",
                                                  "chr1:1-1000",
                                                  "chr1:249250621-249250621",
                                                  "chr22:51304566-51304566",
                                                  "chr17:41196312-41277500",
                                                  "chr13:32889611-32973805"};

static void
RunQuery(benchmark::State &state, const std::string &file, const std::string &region, int region_idx, bool rntuple)
{
   int64_t total_reads_processed = 0;
   Long64_t reads_in_this_run = 0;

   {
      // Suppress the per-call chatter (stopwatch + "Found N records") for the whole run;
      // constructed outside the timed loop so it adds no measured overhead.
      benchutil::ScopedStdoutSuppressor quiet;
      for (auto _ : state) {
         reads_in_this_run = rntuple ? ramntupleview(file.c_str(), region.c_str(), true, false, "perf.root")
                                     : ramview(file.c_str(), region.c_str(), true, false, "perf.root");
         total_reads_processed += reads_in_this_run;
      }
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

static std::vector<int> SelectedRegions(const benchutil::BenchmarkConfig &cfg)
{
   if (cfg.allRegions) {
      std::vector<int> v(kRegions.size());
      std::iota(v.begin(), v.end(), 0);
      return v;
   }
   if (!cfg.regions.empty())
      return cfg.regions;
   return {0, 3, 6, 9};
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);

   // Resolve (or convert once) the files the queries run against.
   const std::string ttree_root = cfg.EnsureTTreeRoot();
   const std::string rntuple_root = cfg.EnsureRNTupleRoot();

   for (int idx : SelectedRegions(cfg)) {
      const std::string region = kRegions[idx % kRegions.size()];

      benchmark::RegisterBenchmark("RegionQuery/TTree/r" + std::to_string(idx), [ttree_root, region,
                                                                                 idx](benchmark::State &state) {
         RunQuery(state, ttree_root, region, idx, false);
      })->Unit(benchmark::kSecond);

      benchmark::RegisterBenchmark("RegionQuery/RNTuple/r" + std::to_string(idx), [rntuple_root, region,
                                                                                   idx](benchmark::State &state) {
         RunQuery(state, rntuple_root, region, idx, true);
      })->Unit(benchmark::kSecond);
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
