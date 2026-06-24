// Columnar read benchmark: RNTuple's differentiator.
//
// Reading a single field (e.g. FLAG or MAPQ) from a columnar store touches only that
// column's pages, whereas a full-record read must decode every column (qname, seq, qual,
// cigar, tags, ...). This benchmark contrasts the two on the same RNTuple so the gap is
// explicit -- the genomics analogue of the HEP "read one branch" pattern.

#include <benchmark/benchmark.h>
#include "benchmark_config.h"
#include "rntuple/RAMNTupleRecord.h"
#include <Rtypes.h>
#include <cstdint>
#include <string>

// Scan only the FLAG column.
static void BM_ColumnarFlagOnly(benchmark::State &state, const std::string &file)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      state.SkipWithError("could not open RNTuple file");
      return;
   }
   const Long64_t n = reader->GetNEntries();

   for (auto _ : state) {
      auto flagView = reader->GetView<uint16_t>("record.flag");
      uint64_t sum = 0;
      for (Long64_t i = 0; i < n; ++i)
         sum += flagView(i);
      benchmark::DoNotOptimize(sum);
   }

   state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * n);
   state.SetLabel(std::to_string(n) + " rows");
}

// Scan only the MAPQ column.
static void BM_ColumnarMapqOnly(benchmark::State &state, const std::string &file)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      state.SkipWithError("could not open RNTuple file");
      return;
   }
   const Long64_t n = reader->GetNEntries();

   for (auto _ : state) {
      auto mapqView = reader->GetView<uint8_t>("record.mapq");
      uint64_t sum = 0;
      for (Long64_t i = 0; i < n; ++i)
         sum += mapqView(i);
      benchmark::DoNotOptimize(sum);
   }

   state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * n);
   state.SetLabel(std::to_string(n) + " rows");
}

// Baseline: materialize the whole record (decodes every column) and read FLAG from it.
static void BM_FullRecordRead(benchmark::State &state, const std::string &file)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      state.SkipWithError("could not open RNTuple file");
      return;
   }
   const Long64_t n = reader->GetNEntries();

   for (auto _ : state) {
      auto recordView = reader->GetView<RAMNTupleRecord>("record");
      uint64_t sum = 0;
      for (Long64_t i = 0; i < n; ++i)
         sum += recordView(i).GetFLAG();
      benchmark::DoNotOptimize(sum);
   }

   state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * n);
   state.SetLabel(std::to_string(n) + " rows");
}

int main(int argc, char **argv)
{
   benchutil::BenchmarkConfig cfg = benchutil::BenchmarkConfig::FromArgs(&argc, argv);
   const std::string rntuple_root = cfg.EnsureRNTupleRoot();

   benchmark::RegisterBenchmark("Columnar/FlagOnly", [rntuple_root](benchmark::State &state) {
      BM_ColumnarFlagOnly(state, rntuple_root);
   })->Unit(benchmark::kMillisecond);
   benchmark::RegisterBenchmark("Columnar/MapqOnly", [rntuple_root](benchmark::State &state) {
      BM_ColumnarMapqOnly(state, rntuple_root);
   })->Unit(benchmark::kMillisecond);
   benchmark::RegisterBenchmark("Columnar/FullRecord", [rntuple_root](benchmark::State &state) {
      BM_FullRecordRead(state, rntuple_root);
   })->Unit(benchmark::kMillisecond);

   benchmark::Initialize(&argc, argv);
   benchmark::RunSpecifiedBenchmarks();
   benchmark::Shutdown();
   return 0;
}
