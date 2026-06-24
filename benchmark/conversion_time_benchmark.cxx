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

// Time SAM->TTree conversion of an on-disk SAM, recording the output size.
static void TimeTTreeConversion(benchmark::State &state, const std::string &sam_file)
{
   const std::string out = "conv_ttree_out.root";
   for (auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet;
         samtoram(sam_file.c_str(), out.c_str(), true, true, true, 1, 0);
      }
      if (std::filesystem::exists(out))
         state.counters["file_size_mb"] = std::filesystem::file_size(out) / (1024.0 * 1024.0);
      std::remove(out.c_str());
   }
}

// Time SAM->RNTuple conversion of an on-disk SAM, recording the output size.
static void TimeRNTupleConversion(benchmark::State &state, const std::string &sam_file, int compression,
                                  unsigned int quality)
{
   const std::string out = "conv_rntuple_out.root";
   for (auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet;
         samtoramntuple(sam_file.c_str(), out.c_str(), true, true, true, compression, quality);
      }
      if (std::filesystem::exists(out))
         state.counters["file_size_mb"] = std::filesystem::file_size(out) / (1024.0 * 1024.0);
      std::remove(out.c_str());
   }
}

// Generated-data sweep: build a synthetic SAM of state.range(0) reads, then time conversion.
static void BM_TTreeGenerated(benchmark::State &state)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "conv_gen_ttree_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   TimeTTreeConversion(state, sam);
   std::remove(sam.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_RNTupleGenerated(benchmark::State &state, int compression, unsigned int quality)
{
   const int num_reads = static_cast<int>(state.range(0));
   const std::string sam = "conv_gen_rntuple_" + std::to_string(num_reads) + ".sam";
   GenerateSAMFile(sam, num_reads);
   TimeRNTupleConversion(state, sam, compression, quality);
   std::remove(sam.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);

   std::cout << "Individual Conversion Time Benchmark" << std::endl;
   std::cout << "====================================" << std::endl;
   std::cout << "Measuring TTree and RNTuple conversion times separately" << std::endl;
   std::cout << std::endl;

   const int compression = cfg.compression;
   const unsigned int quality = cfg.quality;

   if (cfg.HasRealDataset()) {
      const std::string sam = cfg.sam;
      benchmark::RegisterBenchmark("TTree_Conversion/real",
                                   [sam](benchmark::State &state) { TimeTTreeConversion(state, sam); })
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("RNTuple_Conversion/real",
                                   [sam, compression, quality](benchmark::State &state) {
                                      TimeRNTupleConversion(state, sam, compression, quality);
                                   })
         ->Unit(benchmark::kMillisecond);
   } else {
      benchmark::RegisterBenchmark("TTree_Conversion", BM_TTreeGenerated)
         ->Arg(1000)
         ->Arg(10000)
         ->Arg(100000)
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("RNTuple_Conversion",
                                   [compression, quality](benchmark::State &state) {
                                      BM_RNTupleGenerated(state, compression, quality);
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
