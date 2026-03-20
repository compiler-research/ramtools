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
   auto result = ramcore::ComputeStats(kRootFile);

   auto reader = ROOT::RNTupleReader::Open("RAM", kRootFile);
   ASSERT_TRUE(reader != nullptr);

   EXPECT_EQ(result.stats.total_reads, static_cast<uint64_t>(reader->GetNEntries()));
   EXPECT_EQ(result.stats.total_reads, static_cast<uint64_t>(kNumReads));
}

/// Mapped + unmapped must equal total reads (no reads are lost).
TEST_F(RAMStatsTest, MappedPlusUnmappedEqualsTotal)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_EQ(result.stats.mapped_reads + result.stats.unmapped_reads, result.stats.total_reads);
}

/// Forward + reverse strand counts must equal total reads.
TEST_F(RAMStatsTest, StrandCountsEqualTotal)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_EQ(result.stats.forward_strand + result.stats.reverse_strand, result.stats.total_reads);
}

/// Mean read length must be positive for a non-empty file.
TEST_F(RAMStatsTest, MeanReadLengthIsPositive)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_GT(result.stats.mean_read_length, 0.0);
}

/// Total bases must be consistent with mean read length * total reads.
TEST_F(RAMStatsTest, TotalBasesConsistentWithMeanLength)
{
   auto result = ramcore::ComputeStats(kRootFile);
   double expected_bases = result.stats.mean_read_length * result.stats.total_reads;
   // Allow small floating-point rounding tolerance
   EXPECT_NEAR(static_cast<double>(result.stats.total_bases), expected_bases, 1.0);
}

/// Mean mapping quality must be in valid SAM range [0, 255].
TEST_F(RAMStatsTest, MeanMappingQualityInValidRange)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_GE(result.stats.mean_mapping_quality, 0.0);
   EXPECT_LE(result.stats.mean_mapping_quality, 255.0);
}

/// Per-chromosome map must not be empty for a valid SAM file.
TEST_F(RAMStatsTest, ReadsPerChromosomeNotEmpty)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_FALSE(result.stats.reads_per_chromosome.empty());
}

/// Chromosome read counts must sum to at most total_reads
/// (unmapped reads with rname="*" are excluded from the map).
TEST_F(RAMStatsTest, ChromosomeCountsSumToAtMostTotal)
{
   auto result = ramcore::ComputeStats(kRootFile);
   uint64_t chrom_sum = 0;
   for (const auto &[chrom, count] : result.stats.reads_per_chromosome) {
      chrom_sum += count;
   }
   EXPECT_LE(chrom_sum, result.stats.total_reads);
}

/// Print() must produce non-empty output containing expected headers.
TEST_F(RAMStatsTest, PrintProducesOutput)
{
   auto result = ramcore::ComputeStats(kRootFile);
   ASSERT_TRUE(result.ok);
   testing::internal::CaptureStdout();
   result.stats.Print();
   std::string output = testing::internal::GetCapturedStdout();
   EXPECT_FALSE(output.empty());
   EXPECT_NE(output.find("Total reads:"), std::string::npos);
   EXPECT_NE(output.find("Mapped reads:"), std::string::npos);
   EXPECT_NE(output.find("Mean mapping quality:"), std::string::npos);
}

/// ComputeStats on a non-existent file must return ok=false with error message.
TEST(RAMStatsEdgeCases, BadFileReturnsErrorResult)
{
   auto result = ramcore::ComputeStats("nonexistent_file.root");
   EXPECT_FALSE(result.ok);
   EXPECT_FALSE(result.error_message.empty());
}

/// ComputeStats result ok must be true for a valid file.
TEST_F(RAMStatsTest, ValidFileReturnsOk)
{
   auto result = ramcore::ComputeStats(kRootFile);
   EXPECT_TRUE(result.ok);
   EXPECT_TRUE(result.error_message.empty());
}

/// RunRamStats on a valid file must return 0.
TEST_F(RAMStatsTest, RunRamStatsValidFileReturnsZero)
{
   testing::internal::CaptureStdout();
   int ret = ramcore::RunRamStats(kRootFile);
   std::string output = testing::internal::GetCapturedStdout();
   EXPECT_EQ(ret, 0);
   EXPECT_FALSE(output.empty());
}

/// RunRamStats with verbose must print chromosome breakdown.
TEST_F(RAMStatsTest, RunRamStatsVerbosePrintsChromosomes)
{
   testing::internal::CaptureStdout();
   int ret = ramcore::RunRamStats(kRootFile, true);
   std::string output = testing::internal::GetCapturedStdout();
   EXPECT_EQ(ret, 0);
   EXPECT_NE(output.find("Reads per Chromosome"), std::string::npos);
}

/// RunRamStats on a bad file must return 1.
TEST(RAMStatsEdgeCases, RunRamStatsBadFileReturnsOne)
{
   testing::internal::CaptureStderr();
   int ret = ramcore::RunRamStats("nonexistent.root");
   testing::internal::GetCapturedStderr();
   EXPECT_EQ(ret, 1);
}

} // namespace
