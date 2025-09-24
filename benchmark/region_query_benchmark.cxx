#include <benchmark/benchmark.h>
#include "generate_sam_benchmark.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <string>
#include <cstdio>

#ifdef _WIN32
#define NULL_DEVICE "NUL"
#else
#define NULL_DEVICE "/dev/null"
#endif

// Forward declarations
void ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
             const char *perfstatsfilename = "perf.root");
void ramntupleview(const char *file, const char *query, bool cache = true,
                   bool perfstats = false, const char *perfstatsfilename = "perf.root");

class RegionQueryFixture : public benchmark::Fixture {
public:
    void SetUp(const benchmark::State& state) override {
        num_reads_ = static_cast<int>(state.range(0));
        sam_file_ = "region_query_test_" + std::to_string(num_reads_) + ".sam";
        
        GenerateSAMFile(sam_file_, num_reads_);
    }
    
    void TearDown(const benchmark::State&) override {
        std::remove(sam_file_.c_str());
    }

protected:
    int num_reads_;
    std::string sam_file_;
    static constexpr const char* region_ = "chr1:1-100000000";
    
    void suppress_output() {
        freopen(NULL_DEVICE, "w", stdout);
    }
    
    void restore_output() {
        freopen("/dev/tty", "w", stdout);
    }
};

BENCHMARK_DEFINE_F(RegionQueryFixture, TTree)(benchmark::State& state) {
    std::string root_file = "ttree_" + std::to_string(num_reads_) + ".root";
    
    // Convert once before timing
    suppress_output();
    samtoram(sam_file_.c_str(), root_file.c_str(), true, true, true, 1, 0);
    restore_output();
    
    // Benchmark only the query performance
    for (auto _ : state) {
        suppress_output();
        ramview(root_file.c_str(), region_, true, false, "perf.root");
        restore_output();
    }
    
    std::remove(root_file.c_str());
    
    state.counters["reads_per_sec"] = benchmark::Counter(num_reads_, benchmark::Counter::kIsRate);
}

BENCHMARK_DEFINE_F(RegionQueryFixture, RNTuple)(benchmark::State& state) {
    std::string root_file = "rntuple_" + std::to_string(num_reads_) + ".root";
    
    // Convert once before timing
    suppress_output();
    samtoramntuple(sam_file_.c_str(), root_file.c_str(), true, true, true, 505, 0);
    restore_output();
    
    // Benchmark only the query performance
    for (auto _ : state) {
        suppress_output();
        ramntupleview(root_file.c_str(), region_, true, false, "perf.root");
        restore_output();
    }
    
    std::remove(root_file.c_str());
    
    state.counters["reads_per_sec"] = benchmark::Counter(num_reads_, benchmark::Counter::kIsRate);
}

BENCHMARK_REGISTER_F(RegionQueryFixture, TTree)
    ->Args({1000})
    ->Args({10000})
    ->Args({100000})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_REGISTER_F(RegionQueryFixture, RNTuple)
    ->Args({1000})
    ->Args({10000})
    ->Args({100000})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();

