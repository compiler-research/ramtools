#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <TFile.h>
#include <cstdio>
#include <filesystem>

#include "../benchmark/generate_sam_benchmark.h"
#include "ramcore/SamToNTuple.h"
#include "ramcore/RAMStats.h"

namespace {

class RAMStatsTest : public ::testing::Test {
protected:
   static constexpr int kNumReads = 100;
   const char *kSamFile    = "stats_test.sam";
   const char *kRootFile   = "stats_test_rntuple.root";

   void SetUp() override
   {
      GenerateSAMFile(kSamFile, kNumReads);
      std::remove(kRootFile);
      // compression=505 (LZMA), quality_policy=0 — matches existing tests
      samtoramntuple(kSamFile, kRootFile, true, true, true, 505, 0);
   }

   void TearDown() override
   {
      std::remove(kSamFile);
      std::remove(kRootFile);
   }
};

/// Total reads reported must equal the number of entries in the RNTuple.
TEST_F(RAMStatsTest, TotalReadsMatchesNTupleEntries)
{
   auto stats = ramcore::ComputeStats(kRootFile);

   auto reader = ROOT::RNTupleReader::Open("RAM", kRootFile);
   ASSERT_TRUE(reader != nullptr);

   EXPECT_EQ(stats.total_reads, static_cast<uint64_t>(reader->GetNEntries()));
   EXPECT_EQ(stats.total_reads, static_cast<uint64_t>(kNumReads));
}

/// Mapped + unmapped must equal total reads (no reads are lost).
TEST_F(RAMStatsTest, MappedPlusUnmappedEqualsTotal)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   EXPECT_EQ(stats.mapped_reads + stats.unmapped_reads, stats.total_reads);
}

/// Forward + reverse strand counts must equal total reads.
TEST_F(RAMStatsTest, StrandCountsEqualTotal)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   EXPECT_EQ(stats.forward_strand + stats.reverse_strand, stats.total_reads);
}

/// Mean read length must be positive for a non-empty file.
TEST_F(RAMStatsTest, MeanReadLengthIsPositive)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   EXPECT_GT(stats.mean_read_length, 0.0);
}

/// Total bases must be consistent with mean read length * total reads.
TEST_F(RAMStatsTest, TotalBasesConsistentWithMeanLength)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   double expected_bases = stats.mean_read_length * stats.total_reads;
   // Allow small floating-point rounding tolerance
   EXPECT_NEAR(static_cast<double>(stats.total_bases), expected_bases, 1.0);
}

/// Mean mapping quality must be in valid SAM range [0, 255].
TEST_F(RAMStatsTest, MeanMappingQualityInValidRange)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   EXPECT_GE(stats.mean_mapping_quality, 0.0);
   EXPECT_LE(stats.mean_mapping_quality, 255.0);
}

/// Per-chromosome map must not be empty for a valid SAM file.
TEST_F(RAMStatsTest, ReadsPerChromosomeNotEmpty)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   EXPECT_FALSE(stats.reads_per_chromosome.empty());
}

/// Chromosome read counts must sum to at most total_reads
/// (unmapped reads with rname="*" are excluded from the map).
TEST_F(RAMStatsTest, ChromosomeCountsSumToAtMostTotal)
{
   auto stats = ramcore::ComputeStats(kRootFile);
   uint64_t chrom_sum = 0;
   for (const auto &[chrom, count] : stats.reads_per_chromosome) {
      chrom_sum += count;
   }
   EXPECT_LE(chrom_sum, stats.total_reads);
}

/// ComputeStats on a non-existent file must return zero total reads.
TEST(RAMStatsEdgeCases, NonExistentFileReturnsEmpty)
{
   auto stats = ramcore::ComputeStats("nonexistent_file.root");
   EXPECT_EQ(stats.total_reads, 0u);

}

/// ComputeStats on a non-existent file must return zero total reads.
} // namespace
