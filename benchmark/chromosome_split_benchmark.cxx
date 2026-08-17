#include "benchmark_config.h"
#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include <benchmark/benchmark.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>
// Empty => generate synthetic data per benchmark arg; non-empty => split this real SAM.
static std::string g_realSam; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

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

   for ([[maybe_unused]] auto _ : state) {
      std::string bam_file = "bench_st_tmp.bam";
      std::string sorted_bam = "bench_st_sorted.bam";

      std::string cmd = "samtools view -bS ";
      cmd += sam_file;
      cmd += " -o ";
      cmd += bam_file;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools sort ";
      cmd += bam_file;
      cmd += " -o ";
      cmd += sorted_bam;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools index ";
      cmd += sorted_bam;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      for (const auto &chr : chromosomes) {
         cmd = "samtools view -b ";
         cmd += sorted_bam;
         cmd += " ";
         cmd += chr;
         cmd += " > bench_st_";
         cmd += chr;
         cmd += ".bam 2>/dev/null";
         system(cmd.c_str());
      }

      state.counters["size_MB"] = static_cast<double>(benchutil::GetTotalFileSize("bench_st_chr")) / (1024.0 * 1024.0);

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

   for ([[maybe_unused]] auto _ : state) {
      std::string bam_file = "bench_st_mt_tmp.bam";
      std::string sorted_bam = "bench_st_mt_sorted.bam";

      std::string cmd = "samtools view -@ ";
      cmd += std::to_string(num_threads);
      cmd += " -bS ";
      cmd += sam_file;
      cmd += " -o ";
      cmd += bam_file;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools sort -@ ";
      cmd += std::to_string(num_threads);
      cmd += " -m 1G ";
      cmd += bam_file;
      cmd += " -o ";
      cmd += sorted_bam;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools index -@ ";
      cmd += std::to_string(num_threads);
      cmd += " ";
      cmd += sorted_bam;
      cmd += " 2>/dev/null";
      system(cmd.c_str());

      std::vector<std::thread> threads;
      for (const auto &chr : chromosomes) {
         threads.emplace_back([&sorted_bam, &chr]() {
            std::string thread_cmd = "samtools view -@ 2 -b ";
            thread_cmd += sorted_bam;
            thread_cmd += " ";
            thread_cmd += chr;
            thread_cmd += " > bench_st_mt_";
            thread_cmd += chr;
            thread_cmd += ".bam 2>/dev/null";
            system(thread_cmd.c_str());
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

      state.counters["size_MB"] =
         static_cast<double>(benchutil::GetTotalFileSize("bench_st_mt_chr")) / (1024.0 * 1024.0);
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

   for ([[maybe_unused]] auto _ : state) {
      {
         benchutil::ScopedStdoutSuppressor quiet(/*suppress_stderr=*/true);
         samtoramntuple_split_by_chromosome(sam_file.c_str(), /*output_prefix=*/"bench_split_par_out",
                                            /*compression_algorithm=*/505, /*quality_policy=*/1, num_threads);
      }

      state.counters["size_MB"] =
         static_cast<double>(benchutil::GetTotalFileSize("bench_split_par_out_")) / (1024.0 * 1024.0);
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
