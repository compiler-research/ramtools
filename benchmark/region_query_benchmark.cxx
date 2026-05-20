#include <benchmark/benchmark.h>
#include "benchmark_utils.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <Rtypes.h>

static std::string env_or(const char *name, const char *fallback)
{
   const char *v = std::getenv(name);
   return (v && *v) ? std::string(v) : std::string(fallback);
}

// Count records overlapping a region in a BAM or CRAM file, using htslib
// directly (no shell-out). For CRAM, pass a non-empty reference FASTA path.
// Matches `samtools view -c <region>` semantics. Returns -1 on error.
static Long64_t htslib_count_region(const char *file, const char *region,
                                    const char *reference)
{
   samFile *fp = sam_open(file, "r");
   if (!fp) return -1;
   if (reference && *reference) hts_set_fai_filename(fp, reference);
   sam_hdr_t *hdr = sam_hdr_read(fp);
   if (!hdr) { sam_close(fp); return -1; }
   hts_idx_t *idx = sam_index_load(fp, file);
   if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); return -1; }
   hts_itr_t *iter = sam_itr_querys(idx, hdr, region);
   Long64_t n = 0;
   if (iter) {
      bam1_t *b = bam_init1();
      while (sam_itr_next(fp, iter, b) >= 0) ++n;
      bam_destroy1(b);
      hts_itr_destroy(iter);
   }
   hts_idx_destroy(idx);
   sam_hdr_destroy(hdr);
   sam_close(fp);
   return n;
}

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root");
Long64_t ramntupleview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                       const char *perfstatsfilename = "perf.root");

class RegionQueryFixture : public benchmark::Fixture {
public:
   void SetUp(const benchmark::State &state) override
   {
      region_idx_ = static_cast<int>(state.range(0));
   }

   void TearDown(const benchmark::State &) override {}

protected:
   int region_idx_;
   static const std::vector<std::string> regions_;

   void suppress_output() { freopen(NULL_DEVICE, "w", stdout); }
   void restore_output() { freopen("/dev/tty", "w", stdout); }

   const char *get_current_region() const { return regions_[region_idx_ % regions_.size()].c_str(); }
};

// Regions for the complete 1000 Genomes HG00097 chromosome-11 low-coverage
// sample. GRCh37 names contigs without a "chr" prefix; chr11 is 135,006,516 bp.
// Four increasing size tiers: Small (~100 bp), 1 Mb, 10 Mb, full chr11.
const std::vector<std::string> RegionQueryFixture::regions_ = {"11:5000000-5000100",
                                                               "11:5000000-6000000",
                                                               "11:60000000-70000000",
                                                               "11:1-135006516"};

// ---------- TTree variants (per ROOT compression codec) ------------------
// Each variant reads its file path from a distinct env var so all codecs
// can be measured in a single binary run. If the env var is unset or empty,
// the fixture skips with a clear error rather than crashing.
#define DEFINE_TTREE_VARIANT(NAME, ENV_NAME)                                    \
   BENCHMARK_DEFINE_F(RegionQueryFixture, TTree_##NAME)(benchmark::State &state)\
   {                                                                            \
      std::string file = env_or(ENV_NAME, "");                                  \
      if (file.empty()) { state.SkipWithError(ENV_NAME " not set"); return; }  \
      const char *region = get_current_region();                                \
      int64_t total = 0; Long64_t n = 0;                                        \
      for (auto _ : state) {                                                    \
         suppress_output();                                                     \
         n = ramview(file.c_str(), region, true, false, "perf.root");           \
         restore_output();                                                      \
         total += n;                                                            \
      }                                                                         \
      state.SetItemsProcessed(total);                                           \
      state.counters["region_idx"] = region_idx_;                               \
      state.SetLabel(std::to_string(n) + " reads");                             \
   }

DEFINE_TTREE_VARIANT(LZMA, "RAMBENCH_TTREE_LZMA")
DEFINE_TTREE_VARIANT(LZ4,  "RAMBENCH_TTREE_LZ4")
DEFINE_TTREE_VARIANT(ZLIB, "RAMBENCH_TTREE_ZLIB")
DEFINE_TTREE_VARIANT(ZSTD, "RAMBENCH_TTREE_ZSTD")

// ---------- RNTuple variants (per ROOT compression codec) ----------------
#define DEFINE_RNTUPLE_VARIANT(NAME, ENV_NAME)                                     \
   BENCHMARK_DEFINE_F(RegionQueryFixture, RNTuple_##NAME)(benchmark::State &state) \
   {                                                                               \
      std::string file = env_or(ENV_NAME, "");                                     \
      if (file.empty()) { state.SkipWithError(ENV_NAME " not set"); return; }     \
      const char *region = get_current_region();                                   \
      int64_t total = 0; Long64_t n = 0;                                           \
      for (auto _ : state) {                                                       \
         suppress_output();                                                        \
         n = ramntupleview(file.c_str(), region, true, false, "perf.root");        \
         restore_output();                                                         \
         total += n;                                                               \
      }                                                                            \
      state.SetItemsProcessed(total);                                              \
      state.counters["region_idx"] = region_idx_;                                  \
      state.SetLabel(std::to_string(n) + " reads");                                \
   }

DEFINE_RNTUPLE_VARIANT(LZMA, "RAMBENCH_RNTUPLE_LZMA")
DEFINE_RNTUPLE_VARIANT(LZ4,  "RAMBENCH_RNTUPLE_LZ4")
DEFINE_RNTUPLE_VARIANT(ZLIB, "RAMBENCH_RNTUPLE_ZLIB")
DEFINE_RNTUPLE_VARIANT(ZSTD, "RAMBENCH_RNTUPLE_ZSTD")

// ---------- BAM (single config, BGZF default) ----------------------------
BENCHMARK_DEFINE_F(RegionQueryFixture, BAM)(benchmark::State &state)
{
   std::string file = env_or("RAMBENCH_BAM", "data/input.bam");
   if (file.empty()) { state.SkipWithError("RAMBENCH_BAM not set"); return; }
   const char *region = get_current_region();
   int64_t total = 0; Long64_t n = 0;
   for (auto _ : state) {
      n = htslib_count_region(file.c_str(), region, nullptr);
      if (n < 0) { state.SkipWithError("BAM read failed"); break; }
      total += n;
   }
   state.SetItemsProcessed(total);
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(n) + " reads");
}

// ---------- CRAM variants (two presets: 3.1 archive, 3.0 default) --------
#define DEFINE_CRAM_VARIANT(NAME, ENV_NAME)                                     \
   BENCHMARK_DEFINE_F(RegionQueryFixture, CRAM_##NAME)(benchmark::State &state) \
   {                                                                            \
      std::string file = env_or(ENV_NAME, "");                                  \
      std::string ref = env_or("RAMBENCH_REFERENCE", "data/reference.fa");      \
      if (file.empty()) { state.SkipWithError(ENV_NAME " not set"); return; }  \
      const char *region = get_current_region();                                \
      int64_t total = 0; Long64_t n = 0;                                        \
      for (auto _ : state) {                                                    \
         n = htslib_count_region(file.c_str(), region, ref.c_str());            \
         if (n < 0) { state.SkipWithError("CRAM read failed"); break; }         \
         total += n;                                                            \
      }                                                                         \
      state.SetItemsProcessed(total);                                           \
      state.counters["region_idx"] = region_idx_;                               \
      state.SetLabel(std::to_string(n) + " reads");                             \
   }

DEFINE_CRAM_VARIANT(archive, "RAMBENCH_CRAM_ARCHIVE")  // CRAM 3.1 archive preset
DEFINE_CRAM_VARIANT(v30,     "RAMBENCH_CRAM_V30")      // CRAM 3.0 default preset

// ---------- Registration -------------------------------------------------
#define REG(NAME) BENCHMARK_REGISTER_F(RegionQueryFixture, NAME) \
   ->Args({0})->Args({1})->Args({2})->Args({3})->Unit(benchmark::kSecond)

REG(TTree_LZMA);
REG(TTree_LZ4);
REG(TTree_ZLIB);
REG(TTree_ZSTD);
REG(RNTuple_LZMA);
REG(RNTuple_LZ4);
REG(RNTuple_ZLIB);
REG(RNTuple_ZSTD);
REG(BAM);
REG(CRAM_archive);
REG(CRAM_v30);

BENCHMARK_MAIN();
