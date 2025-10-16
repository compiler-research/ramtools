#include <benchmark/benchmark.h>
#include "ramcore/SamToNTuple.h"
#include "generate_sam_benchmark.h"
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <thread>

static void CleanupFiles(const std::string &pattern)
{
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos) {
         std::remove(entry.path().c_str());
      }
   }
}

static size_t GetTotalFileSize(const std::string &pattern)
{
   size_t total = 0;
   for (const auto &entry : std::filesystem::directory_iterator(".")) {
      if (entry.path().filename().string().find(pattern) != std::string::npos) {
         total += std::filesystem::file_size(entry.path());
      }
   }
   return total;
}

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

static void BM_SamtoolsSplit(benchmark::State &state)
{
   int num_reads = state.range(0);
   std::string sam_file = "bench_st_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);
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

      state.counters["size_MB"] = GetTotalFileSize("bench_st_chr") / (1024.0 * 1024.0);

      CleanupFiles("bench_st_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_SamtoolsSplitThreaded(benchmark::State &state)
{
   int num_reads = state.range(0);
   int num_threads = state.range(1);
   std::string sam_file = "bench_st_mt_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);
   auto chromosomes = GetChromosomes(sam_file);

   for (auto _ : state) {
      std::string bam_file = "bench_st_mt_tmp.bam";
      std::string sorted_bam = "bench_st_mt_sorted.bam";

      std::string cmd = "samtools view -@ " + std::to_string(num_threads) +
                       " -bS " + sam_file + " -o " + bam_file + " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools sort -@ " + std::to_string(num_threads) +
            " -m 1G " + bam_file + " -o " + sorted_bam + " 2>/dev/null";
      system(cmd.c_str());

      cmd = "samtools index -@ " + std::to_string(num_threads) +
            " " + sorted_bam + " 2>/dev/null";
      system(cmd.c_str());

      std::vector<std::thread> threads;
      for (const auto &chr : chromosomes) {
         threads.emplace_back([&sorted_bam, &chr]() {
            
            std::string cmd = "samtools view -@ 2 -b " + sorted_bam + " " + chr +
                            " > bench_st_mt_" + chr + ".bam 2>/dev/null";
            system(cmd.c_str());
         });
         
         if (threads.size() >= num_threads) {
            for (auto &t : threads) {
               t.join();
            }
            threads.clear();
         }
      }
      
      for (auto &t : threads) {
         t.join();
      }

      state.counters["size_MB"] = GetTotalFileSize("bench_st_mt_chr") / (1024.0 * 1024.0);
      state.counters["threads"] = num_threads;

      CleanupFiles("bench_st_mt_chr");
      std::remove(bam_file.c_str());
      std::remove(sorted_bam.c_str());
      std::remove((sorted_bam + ".bai").c_str());
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

static void BM_ChromosomeSplitThreads(benchmark::State &state)
{
   int num_reads = state.range(0);
   int num_threads = state.range(1);
   std::string sam_file = "bench_split_par_" + std::to_string(num_reads) + ".sam";

   GenerateSAMFile(sam_file, num_reads);

   FILE *old_stdout = stdout;
   FILE *old_stderr = stderr;

   for (auto _ : state) {
      stdout = fopen("/dev/null", "w");
      stderr = fopen("/dev/null", "w");

      samtoramntuple_split_by_chromosome(sam_file.c_str(), "bench_split_par_out", 505, 1, num_threads);

      fclose(stdout);
      fclose(stderr);
      stdout = old_stdout;
      stderr = old_stderr;

      state.counters["size_MB"] = GetTotalFileSize("bench_split_par_out_") / (1024.0 * 1024.0);
      state.counters["threads"] = num_threads;
      CleanupFiles("bench_split_par_out_");
   }

   std::remove(sam_file.c_str());
   state.counters["reads/s"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
}

BENCHMARK(BM_SamtoolsSplit)
   ->Arg(100000)
   ->Arg(500000)
   ->Arg(1000000)
   ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_SamtoolsSplitThreaded)
   ->Args({100000, 2})
   ->Args({100000, 4})
   ->Args({500000, 2})
   ->Args({500000, 4})
   ->Args({1000000, 2})
   ->Args({1000000, 4})
   ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_ChromosomeSplitThreads)
   ->Args({100000, 2})
   ->Args({100000, 4})
   ->Args({500000, 2})
   ->Args({500000, 4})
   ->Args({1000000, 2})
   ->Args({1000000, 4})
   ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

