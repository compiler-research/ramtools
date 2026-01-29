#include <benchmark/benchmark.h>
#include "generate_sam_benchmark.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <filesystem>
#include <iostream>
#include <cstdio>

#ifdef _WIN32
#define NULL_DEVICE "NUL"
#else
#define NULL_DEVICE "/dev/null"
#endif

static void BM_SamToTTree(benchmark::State &state)
{
   int num_reads = state.range(0);
   std::string sam_file = "ttree_test.sam";
   std::string ttree_file = "ttree_output.root";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {

      FILE *original_stdout = stdout;
      stdout = fopen(NULL_DEVICE, "w");

      samtoram(sam_file.c_str(), ttree_file.c_str(), true, true, true, 1, 0);

      fclose(stdout);
      stdout = original_stdout;

      if (std::filesystem::exists(ttree_file)) {
         auto file_size = std::filesystem::file_size(ttree_file);
         state.counters["file_size_mb"] = file_size / (1024.0 * 1024.0);
      }

      std::remove(ttree_file.c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_SamToRNTuple(benchmark::State &state)
{
   int num_reads = state.range(0);
   std::string sam_file = "rntuple_test.sam";
   std::string rntuple_file = "rntuple_output.root";

   GenerateSAMFile(sam_file, num_reads);

   for (auto _ : state) {

      FILE *original_stdout = stdout;
      stdout = fopen(NULL_DEVICE, "w");

      samtoramntuple(sam_file.c_str(), rntuple_file.c_str(), true, true, true, 505, 0);

      fclose(stdout);
      stdout = original_stdout;

      if (std::filesystem::exists(rntuple_file)) {
         auto file_size = std::filesystem::file_size(rntuple_file);
         state.counters["file_size_mb"] = file_size / (1024.0 * 1024.0);
      }

      std::remove(rntuple_file.c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   std::cout << "Individual Conversion Time Benchmark" << std::endl;
   std::cout << "====================================" << std::endl;
   std::cout << "Measuring TTree and RNTuple conversion times separately" << std::endl;
   std::cout << std::endl;

   ::benchmark::RegisterBenchmark("TTree_Conversion", BM_SamToTTree)
      ->Args({1000})
      ->Args({10000})
      ->Args({100000})
      ->Unit(benchmark::kMillisecond);

   ::benchmark::RegisterBenchmark("RNTuple_Conversion", BM_SamToRNTuple)
      ->Args({1000})
      ->Args({10000})
      ->Args({100000})
      ->Unit(benchmark::kMillisecond);

   ::benchmark::Initialize(&argc, argv);
   ::benchmark::RunSpecifiedBenchmarks();
   return 0;
}
