#include <benchmark/benchmark.h>
#include "generate_sam_benchmark.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <filesystem>
#include <iostream>
#include <cstdio>
#include <vector>
#include <string>

#ifdef _WIN32
#define NULL_DEVICE "NUL"
#else
#define NULL_DEVICE "/dev/null"
#endif

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::vector<std::string> input_files;

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

static void BM_SamToRNTupleRealData(benchmark::State &state)
{
   auto file_index = static_cast<std::size_t>(state.range(0));
   std::string sam_file = input_files[file_index];
   std::string rntuple_file = sam_file + ".rntuple.root";

   // NOLINTNEXTLINE(clang-analyzer-deadcode.DeadStores)
   for (auto _ : state) {

      FILE *original_stdout = stdout;
      stdout = fopen(NULL_DEVICE, "w");

      samtoramntuple(sam_file.c_str(), rntuple_file.c_str(),
                     /*index=*/true,
                     /*split=*/true,
                     /*cache=*/true,
                     /*compression_algorithm=*/505,
                     /*quality_policy=*/0);

      // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
      fclose(stdout);
      stdout = original_stdout;

      if (std::filesystem::exists(rntuple_file)) {
         auto file_size = std::filesystem::file_size(rntuple_file);
         double file_size_mb = static_cast<double>(file_size) / (1024.0 * 1024.0);

         state.counters["file_size_mb"] = file_size_mb;
         state.counters["MB_per_sec"] = benchmark::Counter(file_size_mb, benchmark::Counter::kIsRate);
      }

      std::remove(rntuple_file.c_str());
   }
}

int main(int argc, char **argv)
{
   if (argc > 1) {
      std::cout << "Input Files Conversion Benchmark\n";

      for (int i = 1; i < argc; i++) {
         // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
         input_files.emplace_back(argv[i]);
      }

      for (int i = 0; i < input_files.size(); i++) {
         ::benchmark::RegisterBenchmark((input_files[i]), BM_SamToRNTupleRealData)
            ->Args({i})
            ->Unit(benchmark::kMillisecond);
      }

      ::benchmark::Initialize(&argc, argv);
      ::benchmark::RunSpecifiedBenchmarks();

      return 0;
   }

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
