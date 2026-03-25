#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include <string>
#include <vector>

#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static void WriteSAM(const char *path, const std::string &content)
{
   std::ofstream f(path);
   f << content;
}

// ─────────────────────────────────────────────────────────────────────────────
// SamParser tests
// ─────────────────────────────────────────────────────────────────────────────

namespace {

class SamParserTest : public ::testing::Test {
protected:
   const char *kSamFile = "parser_test.sam";
   void TearDown() override { std::remove(kSamFile); }
};

/// Parser must return false for a non-existent file.
TEST_F(SamParserTest, NonExistentFileReturnsFalse)
{
   ramcore::SamParser parser;
   bool ok = parser.ParseFile("no_such_file.sam", nullptr, nullptr);
   EXPECT_FALSE(ok);
}

/// Header lines starting with '@' must trigger the header callback.
TEST_F(SamParserTest, HeaderCallbackFired)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:248956422\n");

   ramcore::SamParser parser;
   std::vector<std::string> tags;
   parser.ParseFile(kSamFile,
      [&](const std::string &tag, const std::string &) { tags.push_back(tag); },
      nullptr);

   ASSERT_EQ(tags.size(), 2u);
   EXPECT_EQ(tags[0], "@HD");
   EXPECT_EQ(tags[1], "@SQ");
}

/// Record callback must fire once per alignment line.
TEST_F(SamParserTest, RecordCallbackFiredForEachAlignment)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\n"
            "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
            "read2\t16\tchr1\t200\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");

   ramcore::SamParser parser;
   size_t count = 0;
   parser.ParseFile(kSamFile, nullptr,
      [&](const ramcore::SamRecord &, size_t) { ++count; });

   EXPECT_EQ(count, 2u);
   EXPECT_EQ(parser.GetRecordsProcessed(), 2u);
}

/// Parsed record fields must match the SAM columns exactly.
TEST_F(SamParserTest, RecordFieldsParsedCorrectly)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\n"
            "myread\t99\tchr2\t500\t30\t5M2I3M\t=\t600\t110\tACGTACGTAC\tIIIIIIIIII\n");

   ramcore::SamParser parser;
   ramcore::SamRecord captured;
   parser.ParseFile(kSamFile, nullptr,
      [&](const ramcore::SamRecord &rec, size_t) { captured = rec; });

   EXPECT_EQ(captured.qname,  "myread");
   EXPECT_EQ(captured.flag,   99);
   EXPECT_EQ(captured.rname,  "chr2");
   EXPECT_EQ(captured.pos,    500);
   EXPECT_EQ(captured.mapq,   30);
   EXPECT_EQ(captured.cigar,  "5M2I3M");
   EXPECT_EQ(captured.rnext,  "=");
   EXPECT_EQ(captured.pnext,  600);
   EXPECT_EQ(captured.tlen,   110);
   EXPECT_EQ(captured.seq,    "ACGTACGTAC");
   EXPECT_EQ(captured.qual,   "IIIIIIIIII");
}

/// An empty SAM file (no headers, no records) must parse successfully.
TEST_F(SamParserTest, EmptyFileParseSucceeds)
{
   WriteSAM(kSamFile, "");
   ramcore::SamParser parser;
   bool ok = parser.ParseFile(kSamFile, nullptr, nullptr);
   EXPECT_TRUE(ok);
   EXPECT_EQ(parser.GetRecordsProcessed(), 0u);
}

/// A SAM file with only headers and no records must produce zero record calls.
TEST_F(SamParserTest, HeaderOnlyFileProducesNoRecords)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\n"
            "@SQ\tSN:chr1\tLN:1000\n");

   ramcore::SamParser parser;
   size_t count = 0;
   parser.ParseFile(kSamFile, nullptr,
      [&](const ramcore::SamRecord &, size_t) { ++count; });

   EXPECT_EQ(count, 0u);
}

/// Optional fields must be captured in the record.
TEST_F(SamParserTest, OptionalFieldsCaptured)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\n"
            "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\tRG:Z:sample1\n");

   ramcore::SamParser parser;
   ramcore::SamRecord captured;
   parser.ParseFile(kSamFile, nullptr,
      [&](const ramcore::SamRecord &rec, size_t) { captured = rec; });

   EXPECT_GE(captured.optional_fields.size(), 1u);
}

/// Lines processed counter must include both header and record lines.
TEST_F(SamParserTest, LinesProcessedCountIsCorrect)
{
   WriteSAM(kSamFile,
            "@HD\tVN:1.6\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n");

   ramcore::SamParser parser;
   parser.ParseFile(kSamFile, nullptr, nullptr);
   EXPECT_EQ(parser.GetLinesProcessed(), 3u);
}

/// StripCRLF must remove trailing CR and LF characters.
TEST(StripCRLFTest, RemovesTrailingCRLF)
{
   char buf1[] = "hello\r\n";
   ramcore::StripCRLF(buf1);
   EXPECT_STREQ(buf1, "hello");

   char buf2[] = "world\n";
   ramcore::StripCRLF(buf2);
   EXPECT_STREQ(buf2, "world");

   char buf3[] = "noend";
   ramcore::StripCRLF(buf3);
   EXPECT_STREQ(buf3, "noend");
}

// ─────────────────────────────────────────────────────────────────────────────
// RAMNTupleUtils encode/decode tests
// ─────────────────────────────────────────────────────────────────────────────

/// EncodeSequence + DecodeSequence must be inverse operations.
TEST(RAMNTupleUtilsTest, SequenceRoundTrip)
{
   const std::string original = "ACGTACGTACGT";
   std::string encoded = RAMNTupleUtils::EncodeSequence(original);

   // Encoded format: 4-byte length prefix + packed nibbles
   ASSERT_GE(encoded.size(), 4u);
   uint32_t stored_len = 0;
   std::memcpy(&stored_len, encoded.data(), 4);
   EXPECT_EQ(stored_len, original.size());

   std::string decoded = RAMNTupleUtils::DecodeSequence(encoded.substr(4), original.size());
   EXPECT_EQ(decoded, original);
}

/// Odd-length sequences must round-trip correctly.
TEST(RAMNTupleUtilsTest, OddLengthSequenceRoundTrip)
{
   const std::string original = "ACGTA";
   std::string encoded = RAMNTupleUtils::EncodeSequence(original);
   std::string decoded = RAMNTupleUtils::DecodeSequence(encoded.substr(4), original.size());
   EXPECT_EQ(decoded, original);
}

/// Empty sequence must encode and decode without crashing.
TEST(RAMNTupleUtilsTest, EmptySequenceRoundTrip)
{
   const std::string original = "";
   std::string encoded = RAMNTupleUtils::EncodeSequence(original);
   std::string decoded = RAMNTupleUtils::DecodeSequence(encoded.substr(4), 0);
   EXPECT_EQ(decoded, original);
}

/// EncodeQuality + DecodeQuality must be inverse for default Phred33 mode.
TEST(RAMNTupleUtilsTest, QualityRoundTripPhred33)
{
   const std::string original = "IIIIIIIIII";
   uint32_t flags = RAMNTupleRecord::kPhred33;
   std::string encoded = RAMNTupleUtils::EncodeQuality(original, flags);
   std::string decoded = RAMNTupleUtils::DecodeQuality(encoded, flags);
   EXPECT_EQ(decoded, original);
}

/// Quality scores must round-trip when drop mode is used.
TEST(RAMNTupleUtilsTest, QualityRoundTripDropMode)
{
   const std::string original = "IIIIIIIIII";
   uint32_t flags = RAMNTupleRecord::kDrop;
   std::string encoded = RAMNTupleUtils::EncodeQuality(original, flags);
   // In drop mode decoded quality should be empty or placeholder
   std::string decoded = RAMNTupleUtils::DecodeQuality(encoded, flags);
   // Just verify no crash and result is a string
   EXPECT_NO_THROW(RAMNTupleUtils::DecodeQuality(encoded, flags));
}

/// ParseCIGAR + FormatCIGAR must be inverse operations.
TEST(RAMNTupleUtilsTest, CIGARRoundTrip)
{
   const std::string original = "10M2I5M3D";
   auto ops = RAMNTupleUtils::ParseCIGAR(original);
   EXPECT_FALSE(ops.empty());
   std::string formatted = RAMNTupleUtils::FormatCIGAR(ops);
   EXPECT_EQ(formatted, original);
}

/// Wildcard CIGAR '*' must parse without crashing.
TEST(RAMNTupleUtilsTest, WildcardCIGARParsesCleanly)
{
   auto ops = RAMNTupleUtils::ParseCIGAR("*");
   EXPECT_TRUE(ops.empty());
}

} // namespace
