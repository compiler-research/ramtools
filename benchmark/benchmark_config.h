#pragma once

// Shared command-line configuration for the RAMtools benchmark suite.
//
// Lets every benchmark binary run on a user-supplied dataset instead of only on
// synthetic data or hardcoded paths. Recognized flags (all optional):
//
//   --sam=PATH            SAM file to benchmark (default: generate synthetic)
//   --ttree-root=PATH     pre-converted TTree .root (default: derive from --sam)
//   --rntuple-root=PATH   pre-converted RNTuple .root (default: derive from --sam)
//   --compression=N       ROOT compression setting for RNTuple (default 505 = ZSTD-5)
//   --quality=N           quality policy bitmask (default 0)
//   --threads=N           worker threads for split benchmarks (default 4)
//   --reads=N             synthetic read count when no --sam given (default 100000)
//   --regions=all|i,j,k   region indices for region_query (default: 0,3,6,9)
//
// Parse with BenchmarkConfig::FromArgs(&argc, argv) BEFORE benchmark::Initialize so
// these flags are stripped and only --benchmark_* flags reach Google Benchmark.
// When no dataset path is supplied, EnsureSam()/EnsureTTreeRoot()/EnsureRNTupleRoot()
// fall back to GenerateSAMFile() and convert on demand, cleaning up temp files on
// destruction.

#include "benchmark_utils.h"
#include "generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/SamToTTree.h"

#include <cstdlib>
#include <string>
#include <vector>

namespace benchutil {

class BenchmarkConfig {
public:
   std::string sam;         // empty => generate synthetic
   std::string ttreeRoot;   // empty => derive by converting sam
   std::string rntupleRoot; // empty => derive by converting sam
   int compression = 505;   // ROOT compression setting for RNTuple (ZSTD-5)
   unsigned int quality = 0;
   int threads = 4;
   int reads = 100000;       // synthetic read count for fallback
   std::vector<int> regions; // explicit region indices for region_query
   bool allRegions = false;  // --regions=all

   // Parse and strip recognized "--key=value" flags from argv (compacting it in place
   // and lowering *argc), leaving --benchmark_* and other flags for the caller.
   static BenchmarkConfig FromArgs(int *argc, char **argv)
   {
      BenchmarkConfig cfg;
      int out = 1;
      for (int i = 1; i < *argc; ++i) {
         if (!cfg.Consume(argv[i]))
            argv[out++] = argv[i];
      }
      *argc = out;
      return cfg;
   }

   // Path to a SAM file: the user-supplied one, or a freshly generated synthetic file.
   const std::string &EnsureSam()
   {
      if (!sam.empty())
         return sam;
      if (m_genSam.empty()) {
         m_genSam = "bench_gen_" + std::to_string(reads) + ".sam";
         GenerateSAMFile(m_genSam, reads);
         m_created.push_back(m_genSam);
      }
      return m_genSam;
   }

   // Path to an RNTuple .root: the user-supplied one, or one converted from the SAM.
   const std::string &EnsureRNTupleRoot()
   {
      if (!rntupleRoot.empty())
         return rntupleRoot;
      if (m_genRntuple.empty()) {
         m_genRntuple = "bench_gen_rntuple.root";
         const std::string &s = EnsureSam();
         ScopedStdoutSuppressor quiet(true);
         samtoramntuple(s.c_str(), m_genRntuple.c_str(), true, true, true, compression, quality);
         m_created.push_back(m_genRntuple);
      }
      return m_genRntuple;
   }

   // Path to a TTree .root: the user-supplied one, or one converted from the SAM.
   const std::string &EnsureTTreeRoot()
   {
      if (!ttreeRoot.empty())
         return ttreeRoot;
      if (m_genTtree.empty()) {
         m_genTtree = "bench_gen_ttree.root";
         const std::string &s = EnsureSam();
         ScopedStdoutSuppressor quiet(true);
         samtoram(s.c_str(), m_genTtree.c_str(), true, true, true, 1, quality);
         m_created.push_back(m_genTtree);
      }
      return m_genTtree;
   }

   // True when the user supplied a real SAM (vs. falling back to synthetic data).
   bool HasRealDataset() const { return !sam.empty(); }

   void Cleanup()
   {
      for (const auto &f : m_created)
         std::remove(f.c_str());
      m_created.clear();
   }

   ~BenchmarkConfig() { Cleanup(); }

private:
   std::vector<std::string> m_created;
   std::string m_genSam;
   std::string m_genRntuple;
   std::string m_genTtree;

   static bool Match(const std::string &arg, const char *key, std::string &val)
   {
      const std::string prefix = std::string("--") + key + "=";
      if (arg.rfind(prefix, 0) == 0) {
         val = arg.substr(prefix.size());
         return true;
      }
      return false;
   }

   bool Consume(const std::string &arg)
   {
      std::string v;
      if (Match(arg, "sam", v)) {
         sam = v;
         return true;
      }
      if (Match(arg, "ttree-root", v)) {
         ttreeRoot = v;
         return true;
      }
      if (Match(arg, "rntuple-root", v)) {
         rntupleRoot = v;
         return true;
      }
      if (Match(arg, "compression", v)) {
         compression = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, "quality", v)) {
         quality = static_cast<unsigned int>(std::strtoul(v.c_str(), nullptr, 10));
         return true;
      }
      if (Match(arg, "threads", v)) {
         threads = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, "reads", v)) {
         reads = std::atoi(v.c_str());
         return true;
      }
      if (Match(arg, "regions", v)) {
         ParseRegions(v);
         return true;
      }
      return false;
   }

   void ParseRegions(const std::string &v)
   {
      regions.clear();
      allRegions = false;
      if (v == "all") {
         allRegions = true;
         return;
      }
      std::size_t start = 0;
      while (start <= v.size()) {
         const std::size_t comma = v.find(',', start);
         const std::string tok = v.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
         if (!tok.empty())
            regions.push_back(std::atoi(tok.c_str()));
         if (comma == std::string::npos)
            break;
         start = comma + 1;
      }
   }
};

} // namespace benchutil
