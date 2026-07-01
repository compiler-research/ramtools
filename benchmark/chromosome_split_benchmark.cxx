#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

// Empty => generate synthetic data per benchmark arg; non-empty => split this real SAM.
static std::string g_realSam;

static std::vector<std::string> GetChromosomes(const std::string &sam_file)
{
   std::vector<std::string> chroms;
   std::ifstream f(sam_file);
   std::string line;

   while (std::getline(f, line) && line[0] == '@') {
      if (line.find("@SQ\tSN:") == 0) {
         size_t start = 7;
         size_t end = line.find('\t', start);
         chroms.push_back(line.substr(start, end - start));
      }
   }
   return chroms;
}

// Resolve the SAM to operate on: the real dataset, or a freshly generated synthetic file.
// Sets `generated` so the caller knows whether to delete it afterwards.
static std::string PrepareSam(int num_reads, const std::string &gen_name, bool &generated)
{
   if (!g_realSam.empty()) {
      generated = false;
      return g_realSam;
   }
   generated = true;
   GenerateSAMFile(gen_name, num_reads);
   return gen_name;
}

static void BM_SamtoolsSplit(benchmark::State &state)
{
   int num_reads = static_cast<int>(state.range(0));
   bool generated = false;
   std::string sam_file = PrepareSam(num_reads, "bench_st_" + std::to_string(num_reads) + ".sam", generated);
   auto chromosomes = GetChromosomes(sam_file);

   for (auto _ : state) {
      std::string bam_file = "bench_st_tmp.bam";
      std::string sorted_bam = "bench_st_sorted.bam";

      system(("samtools view -bS " + sam_file + " -o " + bam_file + " 2>/dev/null").c_str());

      system(("samtools sort " + bam_file + " -o " + sorted_bam + " 2>/dev/null").c_str());

      system(("samtools index " + sorted_bam + " 2>/dev/null").c_str());

      for (const auto &chr : chromosomes) {
         std::string cmd = "samtools view -b " + sorted_bam + " " + chr + " > bench_st_" + chr + ".bam 2>/dev/null";
         system(cmd.c_str());
      }

      state.counters["size_MB"] = benchutil::GetTotalFileSize("bench_st_chr") / (1024.0 * 1024.0);

      benchutil::CleanupFiles("bench_st_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   if (generated)
      std::remove(sam_file.c_str());
   if (g_realSam.empty())
      state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_SamtoolsSplitThreaded(benchmark::State &state)
{
   int num_reads = static_cast<int>(state.range(0));
   int num_threads = static_cast<int>(state.range(1));
   bool generated = false;
   std::string sam_file = PrepareSam(num_reads, "bench_st_mt_" + std::to_string(num_reads) + ".sam", generated);
   auto chromosomes = GetChromosomes(sam_file);

   for (auto _ : state) {
      std::string bam_file = "bench_st_mt_tmp.bam";
      std::string sorted_bam = "bench_st_mt_sorted.bam";

      std::string cmd =
         "samtools view -@ " + std::to_string(num_threads) + " -bS " + sam_file + " -o " + bam_file + " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools sort -@ " + std::to_string(num_threads) + " -m 1G " + bam_file + " -o " + sorted_bam +
            " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools index -@ " + std::to_string(num_threads) + " " + sorted_bam + " 2>/dev/null";
      system(cmd.c_str());

      std::vector<std::thread> threads;
      for (const auto &chr : chromosomes) {
         threads.emplace_back([&sorted_bam, &chr]() {
            std::string cmd =
               "samtools view -@ 2 -b " + sorted_bam + " " + chr + " > bench_st_mt_" + chr + ".bam 2>/dev/null";
            system(cmd.c_str());
         });

         if (threads.size() >= static_cast<size_t>(num_threads)) {
            for (auto &t : threads) {
               t.join();
            }
            threads.clear();
         }
      }

      for (auto &t : threads) {
         t.join();
      }

      state.counters["size_MB"] = benchutil::GetTotalFileSize("bench_st_mt_chr") / (1024.0 * 1024.0);
      state.counters["threads"] = num_threads;

      benchutil::CleanupFiles("bench_st_mt_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   if (generated)
      std::remove(sam_file.c_str());
   if (g_realSam.empty())
      state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_ChromosomeSplitThreads(benchmark::State &state)
{
   int num_reads = static_cast<int>(state.range(0));
   int num_threads = static_cast<int>(state.range(1));
   bool generated = false;
   std::string sam_file = PrepareSam(num_reads, "bench_split_par_" + std::to_string(num_reads) + ".sam", generated);

   for (auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet(true);
         samtoramntuple_split_by_chromosome(sam_file.c_str(), "bench_split_par_out", 505, 1, num_threads);
      }

      state.counters["size_MB"] = benchutil::GetTotalFileSize("bench_split_par_out_") / (1024.0 * 1024.0);
      state.counters["threads"] = num_threads;
      benchutil::CleanupFiles("bench_split_par_out_");
   }

   if (generated)
      std::remove(sam_file.c_str());
   if (g_realSam.empty())
      state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);
   g_realSam = cfg.sam;

   if (!g_realSam.empty()) {
      // Real dataset: vary only the thread count (read count is fixed by the file).
      benchmark::RegisterBenchmark("BM_SamtoolsSplit/real", BM_SamtoolsSplit)->Arg(0)->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BM_SamtoolsSplitThreaded/real", BM_SamtoolsSplitThreaded)
         ->Args({0, 2})
         ->Args({0, 4})
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BM_ChromosomeSplitThreads/real", BM_ChromosomeSplitThreads)
         ->Args({0, 2})
         ->Args({0, 4})
         ->Unit(benchmark::kMillisecond);
   } else {
      benchmark::RegisterBenchmark("BM_SamtoolsSplit", BM_SamtoolsSplit)
         ->Arg(100000)
         ->Arg(500000)
         ->Arg(1000000)
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BM_SamtoolsSplitThreaded", BM_SamtoolsSplitThreaded)
         ->Args({100000, 2})
         ->Args({100000, 4})
         ->Args({500000, 2})
         ->Args({500000, 4})
         ->Args({1000000, 2})
         ->Args({1000000, 4})
         ->Unit(benchmark::kMillisecond);
      benchmark::RegisterBenchmark("BM_ChromosomeSplitThreads", BM_ChromosomeSplitThreads)
         ->Args({100000, 2})
         ->Args({100000, 4})
         ->Args({500000, 2})
         ->Args({500000, 4})
         ->Args({1000000, 2})
         ->Args({1000000, 4})
         ->Unit(benchmark::kMillisecond);
   }

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
