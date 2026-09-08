#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <Rtypes.h>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include "../benchmark/generate_sam_benchmark.h"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamParser.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
namespace {

const RAMNTupleViewOpts opts = {true, false, ""};
class ramcoreTest : public ::testing::Test {
protected:
   static constexpr const char *kParserTestFile = "test_sam_parser_validation.sam";

   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void SetUp() override
   {
      GenerateSAMFile("samexample.sam", 100);
      std::remove("test_rntuple.root");
   }

   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override
   {
      std::remove("test_rntuple.root");
      std::remove("samexample.sam");
      std::remove(kParserTestFile);
      RAMNTupleRecord::GetIndex()->Clear();
   }

   template <typename OnRecord>
   static size_t ParseSamRecords(std::initializer_list<const char *> records, OnRecord on_record)
   {
      {
         std::ofstream sam(kParserTestFile);
         sam << "@HD\tVN:1.6\n";
         sam << "@SQ\tSN:chr1\tLN:1000\n";
         for (const auto &record : records) {
            sam << record << '\n';
         }
      }

      size_t count = 0;
      ramcore::SamParser parser;
      const bool parsed = parser.ParseFile(
         kParserTestFile, [](const std::string &, const std::string &) {},
         [&](const ramcore::SamRecord &record, size_t) {
            ++count;
            on_record(record);
         });

      EXPECT_TRUE(parsed);
      return count;
   }

   static size_t ParseSamRecords(std::initializer_list<const char *> records)
   {
      return ParseSamRecords(records, [](const ramcore::SamRecord &) {});
   }
};

TEST_F(ramcoreTest, ConversionProducesExpectedEntries)
{
   const char *samFile = "samexample.sam";
   const char *rntupleFile = "test_rntuple.root";

   samtoramntuple(samFile, rntupleFile, true, true, true, 505, 0);

   auto reader = ROOT::RNTupleReader::Open("RAM", rntupleFile);
   ASSERT_NE(reader, nullptr);
   EXPECT_EQ(reader->GetNEntries(), 100);
}

TEST_F(ramcoreTest, RNTupleViewRegionQueries)
{
   const char *rntupleFile = "test_rntuple.root";
   samtoramntuple("samexample.sam", rntupleFile, true, true, true, 505, 0);
   Long64_t hit = ramntupleview(rntupleFile, "chr1:1-1000000", opts);
   EXPECT_GE(hit, 0);

   Long64_t miss = ramntupleview(rntupleFile, "chrNonExistent:1-100", opts);
   EXPECT_EQ(miss, 0);

   Long64_t wildcard = ramntupleview(rntupleFile, "*", opts);
   EXPECT_EQ(wildcard, 100);

   Long64_t empty = ramntupleview(rntupleFile, "", opts);
   EXPECT_EQ(empty, 100);

   Long64_t null = ramntupleview(rntupleFile, nullptr, opts);
   EXPECT_EQ(null, 100);

   Long64_t whole = ramntupleview(rntupleFile, "chr1", opts);
   EXPECT_GE(whole, 0);

   Long64_t single = ramntupleview(rntupleFile, "chr1:500", opts);
   EXPECT_GE(single, 0);

   Long64_t invalid = ramntupleview(rntupleFile, "chr1:abc-def", opts);
   EXPECT_EQ(invalid, 0);

   Long64_t lateChr = ramntupleview(rntupleFile, "chrX:1-100", opts);
   EXPECT_GE(lateChr, 0);

   Long64_t zeroStart = ramntupleview(rntupleFile, "chr1:0-100", opts);
   EXPECT_GE(zeroStart, 0);
}

TEST_F(ramcoreTest, RNTupleViewOpenFailure)
{
   Long64_t count = ramntupleview("nonexistent_file.root", "chr1:1-100", opts);
   EXPECT_EQ(count, 0);
}

TEST_F(ramcoreTest, RNTupleDataIntegrity)
{
   const char *rntupleFile = "test_rntuple.root";
   samtoramntuple("samexample.sam", rntupleFile, true, true, true, 505, 0);

   auto reader = ROOT::RNTupleReader::Open("RAM", rntupleFile);
   ASSERT_NE(reader, nullptr);

   auto viewPos = reader->GetView<int>("record.pos");
   auto viewSeq = reader->GetView<std::string>("record.seq");

   ASSERT_GT(reader->GetNEntries(), 0);

   int firstPos = viewPos(0);
   std::string firstSeq = viewSeq(0);

   EXPECT_GT(firstPos, 0);

   size_t expectedSize = 18 + 4;
   EXPECT_EQ(firstSeq.size(), expectedSize);

   uint32_t storedLen = 0;
   std::memcpy(&storedLen, firstSeq.data(), sizeof(storedLen));
   EXPECT_EQ(storedLen, 36);

   std::cout << "[   INFO   ] Data Integrity - Pos: " << firstPos << ", Encoded Seq Size: " << firstSeq.size()
             << " (Header: " << storedLen << ")" << std::endl;
}

TEST_F(ramcoreTest, RNTupleViewCigarOverlap)
{
   const char *customSam = "test_cigar.sam";
   const char *rntupleFile = "test_cigar.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      sam << "@SQ\tSN:chr1\tLN:10000\n";
      sam << "cigar_read\t0\tchr1\t100\t60\t10M80D10M\t*\t0\t0\t"
          << "AAAAAAAAAAAAAAAAAAAA\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, false, false, false, 505, 0);

   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:140-160", opts), 1);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:201-210", opts), 0);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:50-110", opts), 1);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:1-50", opts), 0);

   std::remove(customSam);
   std::remove(rntupleFile);
}

// `samtools view <region>` returns every alignment overlapping the region,
// secondary and supplementary included -- they cover the locus as much as the
// primary does. Silently returning a subset makes every count from this format
// disagree with samtools.
TEST_F(ramcoreTest, RegionQueryReturnsAllOverlappingAlignments)
{
   const char *customSam = "test_flags.sam";
   const char *rntupleFile = "test_flags.root";

   const std::string seq(50, 'A');

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000\n";
      sam << "primary\t0\tchr1\t100\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      sam << "unmapped\t4\tchr1\t100\t0\t*\t*\t0\t0\t" << std::string(10, 'A') << "\t*\n";
      sam << "secondary\t256\tchr1\t150\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      sam << "supplementary\t2048\tchr1\t200\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, false, false, false, 505, 0);

   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:1-500", opts), 4)
      << "every alignment placed over the region must be reported";

   std::remove(customSam);
   std::remove(rntupleFile);
}

// The index is sparse, so the entry at or before a region start only
// approximates where the region's reads begin. A read starting earlier and
// running long still overlaps; seeking straight to the anchor would step over
// it. Long reads and spliced RNA-seq alignments make this routine.
TEST_F(ramcoreTest, RegionQueryFindsReadsStartingBeforeTheIndexAnchor)
{
   const char *customSam = "test_span.sam";
   const char *rntupleFile = "test_span.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000000\n";

      // Starts at 100000 and reaches 300009 through a long intron.
      sam << "spanning\t0\tchr1\t100000\t60\t10M199990N10M\t*\t0\t0\t" << std::string(20, 'A') << "\t*\n";

      // Enough short reads to put index anchors between it and the query region.
      const std::string seq(50, 'C');
      for (int i = 0; i < 300; ++i) {
         sam << "short" << i << "\t0\tchr1\t" << (290000 + i * 10) << "\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      }
      sam << "inside\t0\tchr1\t300050\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
      // An anchor past the region, so the lookup lands on the anchor before the
      // region instead of falling back to a scan from the first record.
      sam << "after\t0\tchr1\t310000\t60\t50M\t*\t0\t0\t" << seq << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, /*index=*/true, false, false, 505, 0);

   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:300000-300100", opts), 2)
      << "the spanning read overlaps the region and must not be skipped";

   std::remove(customSam);
   std::remove(rntupleFile);
}

// ramntuplescan hands every overlapping row to the callback in file order and
// returns the same count ramntupleview reports; a null callback only counts.
TEST_F(ramcoreTest, ScanVisitsEveryOverlappingRowInOrder)
{
   const char *customSam = "test_scan.sam";
   const char *rntupleFile = "test_scan.root";
   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n";
      sam << "a\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\t*\n";
      sam << "b\t0\tchr1\t150\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "c\t0\tchr1\t300\t60\t4M\t*\t0\t0\tACGT\t*\n";
      sam << "d\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\t*\n";
   }
   samtoramntuple(customSam, rntupleFile, /*index=*/true, false, false, 505, 0);

   auto reader = RAMNTupleRecord::OpenRAMFile(rntupleFile);
   ASSERT_NE(reader, nullptr);
   std::vector<Long64_t> rows;
   auto collect = [&](Long64_t row) { rows.push_back(row); };

   EXPECT_EQ(ramntuplescan(*reader, "chr1:150-200", collect), 1);
   EXPECT_EQ(rows, std::vector<Long64_t>{1});

   rows.clear();
   EXPECT_EQ(ramntuplescan(*reader, "", collect), 4) << "no region means every record";
   EXPECT_EQ(rows, (std::vector<Long64_t>{0, 1, 2, 3}));

   EXPECT_EQ(ramntuplescan(*reader, "chr1", nullptr), 3);
   EXPECT_EQ(ramntuplescan(*reader, "chr1", nullptr), ramntupleview(rntupleFile, "chr1", opts));

   std::remove(customSam);
   std::remove(rntupleFile);
}

// A region query seeks to an index entry and stops at the first record past
// the region. Both assume coordinate order, so a file without it is read end
// to end instead, and gets no index.
TEST_F(ramcoreTest, UnsortedFileIsReadInFullAndNotIndexed)
{
   const char *customSam = "test_unsorted.sam";
   const char *rntupleFile = "test_unsorted.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      // Only a and d overlap chr1:1000-1100. A query that stopped at the first
      // record past the region would report a and never reach d.
      sam << "a\t0\tchr1\t1000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "b\t0\tchr1\t90000\t60\t50M\t*\t0\t0\t" << std::string(50, 'C') << "\t*\n";
      sam << "c\t0\tchr1\t50000\t60\t50M\t*\t0\t0\t" << std::string(50, 'G') << "\t*\n";
      sam << "d\t0\tchr1\t1050\t60\t50M\t*\t0\t0\t" << std::string(50, 'T') << "\t*\n";
   }

   testing::internal::CaptureStderr();
   samtoramntuple(customSam, rntupleFile, /*index=*/true, false, false, 505, 0);
   EXPECT_NE(testing::internal::GetCapturedStderr().find("not in coordinate order"), std::string::npos);

   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:1000-1100", opts), 2);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:50000-50010", opts), 1);
   EXPECT_EQ(ramntupleview(rntupleFile, "chr1:2000-3000", opts), 0);
   EXPECT_FALSE(RAMNTupleRecord::IsCoordinateSorted());
   EXPECT_EQ(RAMNTupleRecord::GetIndex()->Size(), 0U);

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, SortedFileIsIndexedAndMarkedSorted)
{
   const char *customSam = "test_sorted.sam";
   const char *rntupleFile = "test_sorted.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      for (int i = 0; i < 300; ++i)
         sam << "r" << i << "\t0\tchr1\t" << (1000 + i * 100) << "\t60\t50M\t*\t0\t0\t" << std::string(50, 'A')
             << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, /*index=*/true, false, false, 505, 0);

   auto reader = RAMNTupleRecord::OpenRAMFile(rntupleFile);
   ASSERT_NE(reader, nullptr);
   EXPECT_TRUE(RAMNTupleRecord::IsCoordinateSorted());
   EXPECT_GT(RAMNTupleRecord::GetIndex()->Size(), 0U);

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, IndexGetRowsInRange)
{
   RAMNTupleRecord::InitializeRefs();
   auto *index = RAMNTupleRecord::GetIndex();

   //  test with hard coded entries
   index->AddItem(/*refid=*/0, /*pos=*/100, /*row=*/0);
   index->AddItem(/*refid=*/0, /*pos=*/200, /*row=*/1);
   index->AddItem(/*refid=*/0, /*pos=*/300, /*row=*/2);
   index->AddItem(/*refid=*/1, /*pos=*/150, /*row=*/3);

   auto rows = index->GetRowsInRange(/*refid=*/0, /*start=*/150, /*end=*/250);
   ASSERT_EQ(rows.size(), 1U);
   EXPECT_EQ(rows[0], 1);

   auto all = index->GetRowsInRange(/*refid=*/0, /*start=*/0, /*end=*/400);
   EXPECT_EQ(all.size(), 3U);

   auto none = index->GetRowsInRange(/*refid=*/0, /*start=*/400, /*end=*/500);
   EXPECT_TRUE(none.empty());

   auto otherChrom = index->GetRowsInRange(/*refid=*/1, /*start=*/100, /*end=*/200);
   ASSERT_EQ(otherChrom.size(), 1U);
   EXPECT_EQ(otherChrom[0], 3);

   // test with generated entries; the input has to be sorted to get an index
   const char *mockSam = "test_mock_index.sam";
   const char *mockFile = "test_mock_index.root";
   {
      std::ofstream sam(mockSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000000\n";
      for (int i = 0; i < 100; ++i)
         sam << "r" << i << "\t0\tchr1\t" << (1 + i * 1000) << "\t60\t36M\t*\t0\t0\t" << std::string(36, 'A')
             << "\t*\n";
   }
   samtoramntuple(mockSam, mockFile, /*index=*/true, /*split=*/true, /*cache=*/true, /*compression_algorithm=*/505,
                  /*quality_policy=*/0);

   auto reader = RAMNTupleRecord::OpenRAMFile(mockFile);
   ASSERT_NE(reader, nullptr);
   EXPECT_GT(index->Size(), 0U);

   int chr1_refid = RAMNTupleRecord::GetRnameRefs()->GetRefId("chr1");
   EXPECT_GE(chr1_refid, 0);

   auto wideRows = index->GetRowsInRange(/*refid=*/chr1_refid, /*start=*/0, /*end=*/1000000000);
   for (int64_t row : wideRows) {
      EXPECT_GE(row, 0);
      EXPECT_LT(row, 100);
   }

   auto invalidRows = index->GetRowsInRange(/*refid=*/-1, /*start=*/0, /*end=*/1000000000);
   EXPECT_TRUE(invalidRows.empty());

   std::remove(mockSam);
   std::remove(mockFile);
}

TEST_F(ramcoreTest, RecordGetters)
{
   RAMNTupleRecord record;

   record.SetRNEXT("chr1");
   EXPECT_EQ(record.GetRNEXT(), "chr1");
   record.SetRNEXT("=");
   EXPECT_EQ(record.GetRNEXT(), "=");
   record.SetRNEXT("*");
   EXPECT_EQ(record.GetRNEXT(), "*");

   // all 9 CIGAR operations (M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8)
   record.SetCIGAR("1M1I1D1N1S1H1P1=1X");
   EXPECT_EQ(record.GetCIGAR(), "1M1I1D1N1S1H1P1=1X");
   EXPECT_EQ(record.GetNCIGAROP(), 9U);

   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/0), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/1), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/2), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/3), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/4), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/5), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/6), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/7), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/8), 1);
   EXPECT_EQ(record.GetCIGAROPLEN(/*idx=*/9), 0);

   EXPECT_EQ(record.GetCIGAROP(/*idx=*/0), 0);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/1), 1);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/2), 2);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/3), 3);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/4), 4);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/5), 5);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/6), 6);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/7), 7);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/8), 8);
   EXPECT_EQ(record.GetCIGAROP(/*idx=*/9), 0);

   // all 15 IUPAC bases
   record.SetSEQ("ATT");
   EXPECT_EQ(record.GetSEQ(), "ATT");
   record.SetSEQ("=ACMGRSVTWYHKDBN");
   EXPECT_EQ(record.GetSEQ(), "=ACMGRSVTWYHKDBN");
   record.SetSEQ("");
   EXPECT_EQ(record.GetSEQ(), "");

   // kPhred33 (default), returns as-is
   record.SetQUAL("IIIII");
   EXPECT_EQ(record.GetQUAL(), "IIIII");

   // kDrop, always returns *
   RAMNTupleRecord dropRecord;
   dropRecord.SetBit(RAMNTupleRecord::kDrop);
   dropRecord.SetQUAL("IIIII");
   EXPECT_EQ(dropRecord.GetQUAL(), "*");

   // kIlluminaBinning maps a Phred VALUE to one of 0,1,6,15,22,27,33,37,40.
   // SAM writes quality as Phred+33 ASCII, so the encoder must subtract 33
   // before the lookup. These expectations are stated in Phred space and
   // converted, so they cannot silently drift back to indexing by ASCII.
   RAMNTupleRecord binRecord;
   binRecord.SetBit(RAMNTupleRecord::kIlluminaBinning);

   // Two characters, not one: Phred 9 encodes to ASCII 42, which is '*'. A
   // one-base read whose quality is "*" is ambiguous in SAM itself (sentinel
   // vs. Q9), so single-character quality strings are a bad test vector.
   auto phred = [](int q) { return std::string(2, static_cast<char>(q + 33)); };
   auto roundTrip = [&](int q) {
      binRecord.SetQUAL(phred(q));
      const std::string out = binRecord.GetQUAL();
      EXPECT_EQ(out.size(), 2U);
      return static_cast<int>(static_cast<unsigned char>(out[0])) - 33;
   };

   struct Bin {
      int in;
      int want;
   };
   const std::array<Bin, 16> kBins = {{
      {0, 0},
      {1, 1},
      {2, 6},
      {9, 6},
      {10, 15},
      {19, 15},
      {20, 22},
      {24, 22},
      {25, 27},
      {29, 27},
      {30, 33},
      {34, 33},
      {35, 37},
      {39, 37},
      {40, 40},
      {93, 40},
   }};
   for (const auto &c : kBins)
      EXPECT_EQ(roundTrip(c.in), c.want) << "Q" << c.in << " should bin to Q" << c.want;

   // NOTE: binning is NOT monotonically downward -- Illumina maps each bin to
   // a representative value near its middle, so Q2 legitimately becomes Q6.
   // "never raises a quality" is therefore the wrong invariant. What the
   // ASCII-indexing bug actually violated is captured below.

   // Q0 means "no usable base". It must survive as Q0; the old code rewrote it
   // as Q33 (0.05% error), fabricating confidence a variant caller would trust.
   EXPECT_EQ(roundTrip(0), 0) << "a zero-quality base must not be upgraded";

   int previous = -1;
   for (int q = 0; q <= 93; ++q) {
      const int got = roundTrip(q);
      EXPECT_GE(got, previous) << "binning must be monotonic; broke at Q" << q;
      previous = got;
      EXPECT_TRUE(got == 0 || got == 1 || got == 6 || got == 15 || got == 22 || got == 27 || got == 33 || got == 37 ||
                  got == 40)
         << "Q" << q << " produced Q" << got << ", not a legal Illumina bin";
   }

   // "*" means "no quality available"; it is a sentinel, not a Phred string,
   // and must survive rather than be run through the table. Feeding it through
   // mapped '*' (ASCII 42) to bin 40, producing a one-character QUAL against a
   // full-length SEQ -- a malformed record, not merely a wrong one.
   binRecord.SetQUAL("*");
   EXPECT_EQ(binRecord.GetQUAL(), "*");

   // The lossless path must pass the sentinel through untouched too.
   RAMNTupleRecord losslessRecord;
   losslessRecord.SetQUAL("*");
   EXPECT_EQ(losslessRecord.GetQUAL(), "*");

   // Length must be preserved for real quality strings.
   binRecord.SetQUAL("IIIIIIIIII");
   EXPECT_EQ(binRecord.GetQUAL().size(), 10U);
}

} // namespace

TEST_F(ramcoreTest, SmartIndexSkipsUnmappedReads)
{
   const char *customSam = "test_unmapped_index.sam";
   const char *rntupleFile = "test_unmapped_index.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      sam << "mapped1\t0\tchr1\t1000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      sam << "mapped2\t0\tchr1\t2000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);

   Long64_t count = ramntupleview(rntupleFile, "chr1:900-2100", opts);
   EXPECT_EQ(count, 2) << "Both mapped reads should be queryable";

   Long64_t unmapped = ramntupleview(rntupleFile, "*:0-100", opts);
   EXPECT_EQ(unmapped, 0) << "Unmapped reads should not appear in index queries";

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, SmartIndexCreatesEntryAtChromosomeBoundary)
{
   const char *customSam = "test_chrom_boundary.sam";
   const char *rntupleFile = "test_chrom_boundary.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:100000\n";
      sam << "@SQ\tSN:chr2\tLN:100000\n";
      for (int i = 0; i < 50; ++i)
         sam << "chr1_r" << i << "\t0\tchr1\t" << (1000 + i * 100) << "\t60\t50M\t*\t0\t0\t" << std::string(50, 'A')
             << "\t*\n";
      for (int i = 0; i < 50; ++i)
         sam << "chr2_r" << i << "\t0\tchr2\t" << (500 + i * 100) << "\t60\t50M\t*\t0\t0\t" << std::string(50, 'A')
             << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);

   Long64_t chr1_hits = ramntupleview(rntupleFile, "chr1:1000-6000", opts);
   EXPECT_GT(chr1_hits, 0) << "chr1 reads should be queryable";

   Long64_t chr2_hits = ramntupleview(rntupleFile, "chr2:500-5500", opts);
   EXPECT_GT(chr2_hits, 0) << "chr2 reads should be immediately queryable at boundary";

   std::remove(customSam);
   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, SmartIndexRespectsPositionInterval)
{
   const char *customSam = "test_pos_interval.sam";
   const char *rntupleFile = "test_pos_interval.root";

   {
      std::ofstream sam(customSam);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000000\n";
      // Cluster of 200 reads at position 1000 (same position, should not generate 200 index entries)
      for (int i = 0; i < 200; ++i)
         sam << "cluster_r" << i << "\t0\tchr1\t1000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
      // Read far away at 50000 (>10kb gap, should trigger position interval index)
      sam << "far_read\t0\tchr1\t50000\t60\t50M\t*\t0\t0\t" << std::string(50, 'A') << "\t*\n";
   }

   samtoramntuple(customSam, rntupleFile, true, false, false, 505, 0);
   Long64_t cluster = ramntupleview(rntupleFile, "chr1:900-1100", opts);
   EXPECT_EQ(cluster, 200);

   Long64_t far = ramntupleview(rntupleFile, "chr1:49900-50100", opts);
   EXPECT_EQ(far, 1) << "Distant read should be indexed via position interval";

   std::remove(customSam);
   std::remove(rntupleFile);
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, SamParserRejectsMalformedIntegerFields)
{
   EXPECT_EQ(ParseSamRecords({"bad_flag\tabc\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\tIIII",
                              "bad_flag_space\t 1\tchr1\t100\t60\t10M\t*\t0\t0\tACGT\tIIII",
                              "bad_pos\t0\tchr1\t12x\t60\t10M\t*\t0\t0\tACGT\tIIII",
                              "bad_mapq\t0\tchr1\t100\t256\t10M\t*\t0\t0\tACGT\tIIII",
                              "bad_pnext\t0\tchr1\t100\t60\t10M\t*\t-1\t0\tACGT\tIIII",
                              "bad_tlen_min\t0\tchr1\t100\t60\t10M\t*\t0\t-2147483648\tACGT\tIIII",
                              "bad_tlen\t0\tchr1\t100\t60\t10M\t*\t0\t999999999999999999999\tACGT\tIIII",
                              "good\t0\tchr1\t200\t60\t10M\t*\t0\t0\tACGT\tIIII"}),
             1U);
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, SamParserRejectsMalformedCigar)
{
   // An unknown operator, no length, a trailing length, a length past BAM's
   // 28-bit limit, a fractional length, and a second operator with no length.
   for (const char *cigar : {"10Q", "M", "4M2", "268435456M", "1.5M", "4MI"}) {
      const std::string record = std::string("r\t0\tchr1\t100\t60\t") + cigar + "\t*\t0\t0\tACGT\tIIII";
      testing::internal::CaptureStderr();
      EXPECT_EQ(ParseSamRecords({record.c_str()}), 0U) << "CIGAR '" << cigar << "' was accepted";
      EXPECT_NE(testing::internal::GetCapturedStderr().find("malformed CIGAR"), std::string::npos)
         << "no warning for CIGAR '" << cigar << "'";
   }
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, SamParserKeepsValidCigar)
{
   for (const char *cigar : {"*", "4M", "2S4M1I3M1D2N5=1X2H", "268435455M"}) {
      const std::string record = std::string("r\t0\tchr1\t100\t60\t") + cigar + "\t*\t0\t0\tACGT\tIIII";
      std::string seen;
      EXPECT_EQ(ParseSamRecords({record.c_str()}, [&](const ramcore::SamRecord &r) { seen = r.cigar; }), 1U)
         << "CIGAR '" << cigar << "' was rejected";
      EXPECT_EQ(seen, cigar);
   }
}

// A rejected record must not reach the file, and the records around it must.
// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, MalformedCigarRecordIsSkippedByConversion)
{
   const char *samFile = "test_bad_cigar.sam";
   const char *rntupleFile = "test_bad_cigar.root";
   {
      std::ofstream sam(samFile);
      sam << "@HD\tVN:1.6\tSO:coordinate\n";
      sam << "@SQ\tSN:chr1\tLN:1000\n";
      sam << "a\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
      sam << "b\t0\tchr1\t200\t60\t4Q\t*\t0\t0\tACGT\tIIII\n";
      sam << "c\t0\tchr1\t300\t60\t2M2I\t*\t0\t0\tACGT\tIIII\n";
   }

   testing::internal::CaptureStderr();
   samtoramntuple(samFile, rntupleFile, /*index=*/true, false, false, 505, 0);
   EXPECT_NE(testing::internal::GetCapturedStderr().find("malformed CIGAR '4Q' at line 4"), std::string::npos);

   auto reader = RAMNTupleRecord::OpenRAMFile(rntupleFile);
   ASSERT_NE(reader, nullptr);
   ASSERT_EQ(reader->GetNEntries(), 2U);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   EXPECT_EQ(view(0).GetQNAME(), "a");
   EXPECT_EQ(view(0).GetCIGAR(), "4M");
   EXPECT_EQ(view(1).GetQNAME(), "c");
   EXPECT_EQ(view(1).GetCIGAR(), "2M2I");

   std::remove(samFile);
   std::remove(rntupleFile);
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, SamParserParsesValidIntegerBoundaries)
{
   ramcore::SamRecord parsed;

   ASSERT_EQ(ParseSamRecords({"boundary\t65535\tchr1\t0\t255\t10M\t*\t0\t-2147483647\tACGT\tIIII"},
                             [&](const ramcore::SamRecord &record) { parsed = record; }),
             1U);
   EXPECT_EQ(parsed.flag, 65535);
   EXPECT_EQ(parsed.pos, 0);
   EXPECT_EQ(parsed.mapq, 255);
   EXPECT_EQ(parsed.pnext, 0);
   EXPECT_EQ(parsed.tlen, std::numeric_limits<int>::min() + 1);
}

// Resolving a region's reference must be a lookup, never an insert: a query for a
// contig the file does not contain used to append it to fRefVec, giving every
// later record's refid a different meaning than the one it was written with.
// NOLINTNEXTLINE(misc-use-internal-linkage)
TEST_F(ramcoreTest, InvalidChromosomeDoesNotPolluteFRefVec)
{
   const char *samFile = "samexample.sam";
   const char *rntupleFile = "test_rntuple.root";

   samtoramntuple(samFile, rntupleFile, true, true, true, 505, 0);

   const auto refsBefore = RAMNTupleRecord::GetRnameRefs()->Size();

   testing::internal::CaptureStdout();
   testing::internal::CaptureStderr();
   const Long64_t found = ramntupleview(rntupleFile, "chrINVALID:100-200", opts);
   testing::internal::GetCapturedStdout();
   testing::internal::GetCapturedStderr();

   EXPECT_EQ(found, 0) << "a reference that is not in the file cannot have records";
   EXPECT_EQ(refsBefore, RAMNTupleRecord::GetRnameRefs()->Size())
      << "Invalid chromosome 'chrINVALID' was inserted into fRefVec (regression of issue #23)";

   std::remove(rntupleFile);
}

TEST_F(ramcoreTest, QUALEncodingDecodingModes)
{
   const char *samFile = "test_qual.sam";
   const char *ramFile = "test_qual.root";

   std::string seq = "AAAAAAAAAA";
   std::string qual = "@@@FBIEDH!";

   {
      std::ofstream sam(samFile);
      sam << "@HD\tVN:1.6\tSO:unsorted\n";
      sam << "@SQ\tSN:chr1\tLN:1000\n";
      sam << "read1\t0\tchr1\t100\t60\t10M\t*\t0\t0\t" << seq << "\t" << qual << "\n";
   }

   // No compression of QUAL field, stored as it is
   samtoramntuple(samFile, ramFile, /*index=*/true, /*split=*/false, /*cache=*/false, /*compression_algorithm=*/505,
                  /*quality_policy=*/0);

   {
      auto reader = ROOT::RNTupleReader::Open("RAM", ramFile);
      ASSERT_NE(reader, nullptr);

      auto view = reader->GetView<RAMNTupleRecord>("record");
      const auto &rec = view(0);

      std::string decoded = rec.GetQUAL();

      // MUST be identical
      EXPECT_EQ(decoded, qual) << "QUAL should remain unchanged without compression";
   }

   std::remove(ramFile);

   // Drop the quailty field , should stored as "*"
   samtoramntuple(samFile, ramFile, /*index=*/true, /*split=*/false, /*cache=*/false, /*compression_algorithm=*/505,
                  /*quality_policy=*/RAMNTupleRecord::kDrop);

   {
      auto reader = ROOT::RNTupleReader::Open("RAM", ramFile);
      ASSERT_NE(reader, nullptr);

      auto view = reader->GetView<RAMNTupleRecord>("record");
      const auto &rec = view(0);

      EXPECT_EQ(rec.GetQUAL(), "*") << "QUAL should be dropped";
   }

   std::remove(ramFile);

   // Illumina binning
   samtoramntuple(samFile, ramFile, /*index=*/true, /*split=*/false, /*cache=*/false, /*compression_algorithm=*/505,
                  /*quality_policy=*/RAMNTupleRecord::kIlluminaBinning);

   {
      auto reader = ROOT::RNTupleReader::Open("RAM", ramFile);
      ASSERT_NE(reader, nullptr);

      auto view = reader->GetView<RAMNTupleRecord>("record");
      const auto &rec = view(0);

      std::string decoded = rec.GetQUAL();

      // Must NOT be same (lossy compression done during encoding)
      EXPECT_NE(decoded, qual) << "Binning must change quality values";

      // Length must stay same
      EXPECT_EQ(decoded.size(), qual.size());

      // Must still be valid printable ASCII
      for (char c : decoded) {
         EXPECT_GE(c, 33);
         EXPECT_LE(c, 126);
      }
   }

   std::remove(samFile);
   std::remove(ramFile);
}
