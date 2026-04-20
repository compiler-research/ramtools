#include "ramcore/BamtoNTuple.h"

#include "generate_sam_benchmark.h"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"

#include <gtest/gtest.h>

#include <ROOT/RNTupleReader.hxx>

#include <htslib/sam.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>

namespace {

void GenerateBAMFile(const char *bam_path, int num_reads)
{
   const char *sam_path = "temp_for_bam.sam";
   GenerateSAMFile(sam_path, num_reads);

   samFile *in = sam_open(sam_path, "r");
   samFile *out = sam_open(bam_path, "wb");
   sam_hdr_t *hdr = sam_hdr_read(in);

   if (sam_hdr_write(out, hdr) < 0) {
      throw std::runtime_error("Failed to write BAM header");
   }

   bam1_t *rec = bam_init1();
   while (sam_read1(in, hdr, rec) >= 0) {
      if (sam_write1(out, hdr, rec) < 0) {
         throw std::runtime_error("Failed to write BAM record");
      }
   }

   bam_destroy1(rec);
   sam_hdr_destroy(hdr);
   sam_close(in);
   sam_close(out);
   std::remove(sam_path);
}

void GenerateRichBAMFile(const char *bam_path)
{
   samFile *out = sam_open(bam_path, "wb");

   sam_hdr_t *hdr = sam_hdr_init();

   // Bypassed varargs entirely by passing the structured header string
   const std::string header_text = "@HD\tVN:1.6\tSO:coordinate\n"
                                   "@SQ\tSN:chr1\tLN:100000\n"
                                   "@SQ\tSN:chr2\tLN:80000\n";

   if (sam_hdr_add_lines(hdr, header_text.c_str(), header_text.length()) < 0) {
      throw std::runtime_error("Failed to add rich BAM header lines");
   }

   if (sam_hdr_write(out, hdr) < 0) {
      throw std::runtime_error("Failed to write rich BAM header");
   }

   bam1_t *rec = bam_init1();

   // Mapped read with all aux tag types: A, c, C, s, S, i, I, f, Z, H
   {
      const char *qname = "read_tags";
      const char *seq = "ACGTACGTACGT";
      const int len = static_cast<int>(strlen(seq));
      std::string qual(len, 'I');

      std::array<uint32_t, 1> cigar = {static_cast<uint32_t>(len) << 4 | 0};

      bam_set1(rec, strlen(qname), qname, 0, 0, 1000, 60, 1, cigar.data(), -1, -1, 0, len, seq, qual.c_str(), 128);

      int8_t val_c = -5;
      bam_aux_append(rec, "Xc", 'c', sizeof(val_c), reinterpret_cast<const uint8_t *>(&val_c));

      uint8_t val_C = 200;
      bam_aux_append(rec, "XC", 'C', sizeof(val_C), &val_C);

      int16_t val_s = -1000;
      bam_aux_append(rec, "Xs", 's', sizeof(val_s), reinterpret_cast<const uint8_t *>(&val_s));

      uint16_t val_S = 50000;
      bam_aux_append(rec, "XS", 'S', sizeof(val_S), reinterpret_cast<const uint8_t *>(&val_S));

      int32_t val_i = -100000;
      bam_aux_append(rec, "Xi", 'i', sizeof(val_i), reinterpret_cast<const uint8_t *>(&val_i));

      uint32_t val_I = 3000000;
      bam_aux_append(rec, "XI", 'I', sizeof(val_I), reinterpret_cast<const uint8_t *>(&val_I));

      float val_f = 3.14F;
      bam_aux_append(rec, "Xf", 'f', sizeof(val_f), reinterpret_cast<const uint8_t *>(&val_f));

      uint8_t val_A = 'Q';
      bam_aux_append(rec, "XA", 'A', 1, &val_A);

      const char *val_Z = "hello";
      bam_aux_append(rec, "XZ", 'Z', static_cast<int>(strlen(val_Z) + 1), reinterpret_cast<const uint8_t *>(val_Z));

      const char *val_H = "1AE301";
      bam_aux_append(rec, "XH", 'H', static_cast<int>(strlen(val_H) + 1), reinterpret_cast<const uint8_t *>(val_H));

      if (sam_write1(out, hdr, rec) < 0)
         throw std::runtime_error("Failed to write record 1");
   }

   // Unmapped read, no CIGAR, no quality
   {
      const char *qname = "read_unmapped";
      const char *seq = "NNNNNN";
      const int len = static_cast<int>(strlen(seq));

      bam_set1(rec, strlen(qname), qname, 4, -1, -1, 0, 0, nullptr, -1, -1, 0, len, seq, nullptr, 0);

      if (sam_write1(out, hdr, rec) < 0)
         throw std::runtime_error("Failed to write record 2");
   }

   // Mate on different chromosome
   {
      const char *qname = "read_diffchrom";
      const char *seq = "ACGTACGT";
      const int len = static_cast<int>(strlen(seq));
      std::string qual(len, 'F');

      std::array<uint32_t, 1> cigar = {static_cast<uint32_t>(len) << 4 | 0};

      bam_set1(rec, strlen(qname), qname, 1, 0, 500, 30, 1, cigar.data(), 1, 2000, 0, len, seq, qual.c_str(), 32);

      if (sam_write1(out, hdr, rec) < 0)
         throw std::runtime_error("Failed to write record 3");
   }

   // Mate on same chromosome, multi-op CIGAR (4M4S)
   {
      const char *qname = "read_samechrom";
      const char *seq = "TTTTGGGG";
      const int len = static_cast<int>(strlen(seq));
      std::string qual(len, 'A');

      std::array<uint32_t, 2> cigar = {4U << 4 | 0, 4U << 4 | 4};

      bam_set1(rec, strlen(qname), qname, 1, 0, 600, 40, 2, cigar.data(), 0, 800, 200, len, seq, qual.c_str(), 32);

      if (sam_write1(out, hdr, rec) < 0)
         throw std::runtime_error("Failed to write record 4");
   }

   // Bulk reads for index exercising
   for (int r = 0; r < 1001; ++r) {
      std::string qname = "bulk_" + std::to_string(r);
      const char *seq = "ACGT";
      const int len = 4;
      const char *qual = "IIII";

      std::array<uint32_t, 1> cigar = {4U << 4 | 0};

      bam_set1(rec, qname.size(), qname.c_str(), 0, 0, 100 + r, 50, 1, cigar.data(), -1, -1, 0, len, seq, qual, 0);

      if (sam_write1(out, hdr, rec) < 0)
         throw std::runtime_error("Failed to write bulk record");
   }

   bam_destroy1(rec);
   sam_hdr_destroy(hdr);
   sam_close(out);
}

class BamToNTupleTest : public ::testing::Test {
protected:
   void SetUp() override
   {
      GenerateBAMFile(/*bam_path=*/"test.bam", /*num_reads=*/100);
      GenerateSAMFile(/*sam_path=*/"test_compare.sam", /*num_reads=*/100);
      std::remove("test_bam.ram");
      std::remove("test_sam.ram");
   }

   void TearDown() override
   {
      std::remove("test_bam.ram");
      std::remove("test_sam.ram");
      std::remove("test.bam");
      std::remove("test_compare.sam");
      std::remove("test_rich.bam");
      std::remove("test_rich.ram");
   }
};

TEST_F(BamToNTupleTest, ConversionProducesEntries)
{
   bamtoramntuple("test.bam", "test_bam.ram", true, false, true, 505, 0U);

   ASSERT_TRUE(std::filesystem::exists("test_bam.ram"));

   auto reader = ROOT::RNTupleReader::Open("RAM", "test_bam.ram");
   ASSERT_NE(reader, nullptr);
   EXPECT_GT(reader->GetNEntries(), 0);
}

TEST_F(BamToNTupleTest, SameEntryCountAsSAMPath)
{
   bamtoramntuple("test.bam", "test_bam.ram", true, false, true, 505, 0U);

   samtoramntuple("test_compare.sam", "test_sam.ram", true, true, true, 505, 0U);

   auto bamReader = ROOT::RNTupleReader::Open("RAM", "test_bam.ram");
   auto samReader = ROOT::RNTupleReader::Open("RAM", "test_sam.ram");

   ASSERT_NE(bamReader, nullptr);
   ASSERT_NE(samReader, nullptr);
   EXPECT_EQ(bamReader->GetNEntries(), samReader->GetNEntries());
}

TEST_F(BamToNTupleTest, RegionQueryWorks)
{
   bamtoramntuple("test.bam", "test_bam.ram", true, false, true, 505, 0U);

   const std::array<const char *, 3> regions = {"chr1:100-200", "chr2:500-1000", "chr5:1000-5000"};

   for (const char *region : regions) {
      testing::internal::CaptureStdout();
      const Long64_t count = ramntupleview("test_bam.ram", region, true, false, nullptr);
      testing::internal::GetCapturedStdout();
      EXPECT_GE(count, 0);
   }
}

TEST_F(BamToNTupleTest, MandatoryFieldsPresent)
{
   bamtoramntuple("test.bam", "test_bam.ram", true, false, true, 505, 0U);

   auto reader = ROOT::RNTupleReader::Open("RAM", "test_bam.ram");
   ASSERT_NE(reader, nullptr);

   const auto &desc = reader->GetDescriptor();
   EXPECT_NE(desc.FindFieldId("record"), ROOT::kInvalidDescriptorId);
   EXPECT_GT(reader->GetNEntries(), 0);
}

TEST_F(BamToNTupleTest, RichBAMCoversAllTagTypes)
{
   GenerateRichBAMFile(/*bam_path=*/"test_rich.bam");

   bamtoramntuple("test_rich.bam", "test_rich.ram", true, false, true, 505, 0U);

   auto reader = ROOT::RNTupleReader::Open("RAM", "test_rich.ram");
   ASSERT_NE(reader, nullptr);
   EXPECT_EQ(reader->GetNEntries(), 1005);
}

TEST_F(BamToNTupleTest, InvalidFileReturnsGracefully)
{
   testing::internal::CaptureStderr();
   bamtoramntuple("nonexistent.bam", "out.ram", true, false, true, 505, 0U);
   std::string err = testing::internal::GetCapturedStderr();
   EXPECT_NE(err.find("Cannot open BAM"), std::string::npos);
   EXPECT_FALSE(std::filesystem::exists("out.ram"));
}

TEST_F(BamToNTupleTest, ConversionWithoutIndex)
{
   bamtoramntuple("test.bam", "test_bam.ram", false, false, true, 505, 0U);

   auto reader = ROOT::RNTupleReader::Open("RAM", "test_bam.ram");
   ASSERT_NE(reader, nullptr);
   EXPECT_GT(reader->GetNEntries(), 0);
}

TEST_F(BamToNTupleTest, UnmappedAndMissingQualityHandled)
{
   GenerateRichBAMFile(/*bam_path=*/"test_rich.bam");

   bamtoramntuple("test_rich.bam", "test_rich.ram", true, false, true, 505, 0U);

   testing::internal::CaptureStdout();
   const Long64_t count = ramntupleview("test_rich.ram", "chr1:100-2000", true, false, nullptr);
   testing::internal::GetCapturedStdout();
   EXPECT_GE(count, 0);
}

TEST_F(BamToNTupleTest, DifferentChromMateHandled)
{
   GenerateRichBAMFile(/*bam_path=*/"test_rich.bam");

   bamtoramntuple("test_rich.bam", "test_rich.ram", true, false, true, 505, 0U);

   testing::internal::CaptureStdout();
   const Long64_t count = ramntupleview("test_rich.ram", "chr2:400-600", true, false, nullptr);
   testing::internal::GetCapturedStdout();
   EXPECT_GE(count, 0);
}

} // namespace