#include <benchmark/benchmark.h>
#include "generate_sam_benchmark.h"
#include <cstdio> // For std::remove

static void BM_GenerateSAM(benchmark::State &state)
{
   int num_reads = state.range(0);
   for (auto _ : state) {
      GenerateSAMFile("benchmark_temp.sam", num_reads);
      std::remove("benchmark_temp.sam");
   }

   state.counters["reads_per_second"] = benchmark::Counter(num_reads, benchmark::Counter::kIsRate);
   state.counters["bytes_per_second"] = benchmark::Counter(num_reads * 200, benchmark::Counter::kIsRate);
}

BENCHMARK(BM_GenerateSAM)->Range(100, 10000);
BENCHMARK_MAIN();
