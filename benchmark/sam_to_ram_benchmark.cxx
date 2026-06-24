#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/SamToTTree.h"
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

// Convert one on-disk SAM to both TTree and RNTuple, comparing output sizes.
static void CompareConversion(benchmark::State &state, const std::string &sam_file, int compression,
                              unsigned int quality)
{
   const std::string ttree_file = "cmp_ttree.root";
   const std::string rntuple_file = "cmp_rntuple.root";

   for (auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet;
         samtoram(sam_file.c_str(), ttree_file.c_str(), true, true, true, 1, quality);
         samtoramntuple(sam_file.c_str(), rntuple_file.c_str(), true, true, true, compression, quality);
      }

      if (std::filesystem::exists(ttree_file))
         state.counters["ttree_size_mb"] = std::filesystem::file_size(ttree_file) / (1024.0 * 1024.0);
      if (std::filesystem::exists(rntuple_file))
         state.counters["rntuple_size_mb"] = std::filesystem::file_size(rntuple_file) / (1024.0 * 1024.0);
      if (state.counters["ttree_size_mb"] > 0 && state.counters["rntuple_size_mb"] > 0)
         state.counters["compression_ratio"] = state.counters["ttree_size_mb"] / state.counters["rntuple_size_mb"];

      std::remove(ttree_file.c_str());
      std::remove(rntuple_file.c_str());
   }
}

// Generated-data sweep: build a synthetic SAM of state.range(0) reads, then compare.
static void BM_Generated(benchmark::State &state, int compression, unsigned int quality)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "cmp_gen_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   CompareConversion(state, sam, compression, quality);
   std::remove(sam.c_str());
   state.counters["reads"] = num_reads;
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);

   std::cout << "SAM to RAM Conversion Benchmark" << std::endl;
   std::cout << "Comparing TTree vs RNTuple performance and file sizes" << std::endl;
   std::cout << std::endl;

   const int compression = cfg.compression;
   const unsigned int quality = cfg.quality;

   if (cfg.HasRealDataset()) {
      const std::string sam = cfg.sam;
      benchmark::RegisterBenchmark("BM_SamToRamComparison/real",
                                   [sam, compression, quality](benchmark::State &state) {
                                      CompareConversion(state, sam, compression, quality);
                                   })
         ->Unit(benchmark::kMillisecond);
   } else {
      benchmark::RegisterBenchmark("BM_SamToRamComparison",
                                   [compression, quality](benchmark::State &state) {
                                      BM_Generated(state, compression, quality);
                                   })
         ->Arg(1000)
         ->Arg(10000)
         ->Arg(100000)
         ->Unit(benchmark::kMillisecond);
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
