#include <benchmark/benchmark.h>
#include "generate_sam_benchmark.h"
#include "benchmark_utils.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <filesystem>
#include <iostream>
#include <cstdio>
#include <cstring>

static void BM_SamToRamComparison(benchmark::State &state)
{
   int num_reads = state.range(0);

   std::string sam_file = "benchmark_test.sam";
   std::string ttree_file = "benchmark_ttree.root";
   std::string rntuple_file = "benchmark_rntuple.root";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {

      FILE *original_stdout = stdout;
      stdout = fopen(NULL_DEVICE, "w");

      samtoram(sam_file.c_str(), ttree_file.c_str(), true, true, true, 1, 0);

      samtoramntuple(sam_file.c_str(), rntuple_file.c_str(), true, true, true, 505, 0);

      fclose(stdout);
      stdout = original_stdout;

      if (std::filesystem::exists(ttree_file)) {
         auto ttree_size = std::filesystem::file_size(ttree_file);
         state.counters["ttree_size_mb"] = ttree_size / (1024.0 * 1024.0);
      }

      if (std::filesystem::exists(rntuple_file)) {
         auto rntuple_size = std::filesystem::file_size(rntuple_file);
         state.counters["rntuple_size_mb"] = rntuple_size / (1024.0 * 1024.0);
      }

      if (state.counters["ttree_size_mb"] > 0 && state.counters["rntuple_size_mb"] > 0) {
         state.counters["compression_ratio"] = state.counters["ttree_size_mb"] / state.counters["rntuple_size_mb"];
      }
   }

   std::remove(sam_file.c_str());
   std::remove(ttree_file.c_str());
   std::remove(rntuple_file.c_str());

   state.counters["reads"] = num_reads;
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   std::cout << "SAM to RAM Conversion Benchmark" << std::endl;
   std::cout << "Comparing TTree vs RNTuple performance and file sizes" << std::endl;
   std::cout << std::endl;

   ::benchmark::RegisterBenchmark("BM_SamToRamComparison", BM_SamToRamComparison)
      ->Args({1000})
      ->Args({10000})
      ->Args({100000})
      ->Unit(benchmark::kMillisecond);

   ::benchmark::Initialize(&argc, argv);
   ::benchmark::RunSpecifiedBenchmarks();
   return 0;
}
