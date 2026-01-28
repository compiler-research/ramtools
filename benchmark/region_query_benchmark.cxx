#include <benchmark/benchmark.h>
#include "benchmark_utils.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <string>
#include <cstdio>
#include <vector>
#include <Rtypes.h>

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root");
Long64_t ramntupleview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");

class RegionQueryFixture : public benchmark::Fixture {
public:
   void SetUp(const benchmark::State &state) override
   {
      region_idx_ = static_cast<int>(state.range(0));

      sam_file_ = "/media/aditya/213e0e46-6f86-4288-8b79-74851c34314f/output_big.sam";
      ttree_root_file_ = "/media/aditya/213e0e46-6f86-4288-8b79-74851c34314f/output_big_lzma.root";
      rntuple_root_file_ = "/home/aditya/ramtools/build/rntuple.root";
   }

   void TearDown(const benchmark::State &) override {}

protected:
   int region_idx_;
   std::string sam_file_;
   std::string ttree_root_file_;
   std::string rntuple_root_file_;

   static const std::vector<std::string> regions;

   void suppress_output() { freopen(NULL_DEVICE, "w", stdout); }
   void restore_output() { freopen("/dev/tty", "w", stdout); }

   [[nodiscard]] const char *get_current_region() const { return regions[region_idx_ % regions.size()].c_str(); }
};

const std::vector<std::string> RegionQueryFixture::regions = {"1:1000000-1001000",
                                                              "2:5000000-5010000",
                                                              "X:100000-150000",
                                                              "1:1000000-2000000",
                                                              "5:10000000-15000000",
                                                              "10:50000000-60000000",
                                                              "1:1-50000000",
                                                              "2:1-100000000",
                                                              "7:50000000-150000000",
                                                              "21:1-48129895",
                                                              "MT:1-16571",
                                                              "Y:2600000-2700000",
                                                              "GL000227.1:1-100000",
                                                              "1:1-1000",
                                                              "1:249250621-249250621",
                                                              "22:51304566-51304566",
                                                              "17:41196312-41277500",
                                                              "13:32889611-32973805"};

BENCHMARK_DEFINE_F(RegionQueryFixture, TTree)(benchmark::State &state)
{
   const char *region = get_current_region();
   int64_t total_reads_processed = 0;
   Long64_t reads_in_this_run = 0;

   for (auto _ : state) {
      suppress_output();
      reads_in_this_run = ramview(ttree_root_file_.c_str(), region, true, false, "perf.root");
      restore_output();

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

BENCHMARK_DEFINE_F(RegionQueryFixture, RNTuple)(benchmark::State &state)
{
   const char *region = get_current_region();
   int64_t total_reads_processed = 0;
   Long64_t reads_in_this_run = 0;

   for (auto _ : state) {
      suppress_output();
      reads_in_this_run = ramntupleview(rntuple_root_file_.c_str(), region, true, false, "perf.root");
      restore_output();

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

BENCHMARK_REGISTER_F(RegionQueryFixture, TTree)->DenseRange(0, 17, 1)->Unit(benchmark::kSecond);
BENCHMARK_REGISTER_F(RegionQueryFixture, RNTuple)->DenseRange(0, 17, 1)->Unit(benchmark::kSecond);

BENCHMARK_MAIN();
