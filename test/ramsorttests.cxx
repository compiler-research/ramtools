#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <cstdio>

#include "../benchmark/generate_sam_benchmark.h"
#include "ramcore/RAMSort.h"
#include "ramcore/SamToNTuple.h"

namespace {

class RAMSortTest : public ::testing::Test {
protected:
   static constexpr int kNumReads = 200;
   const char *kSamFile      = "sort_test.sam";
   const char *kUnsortedFile = "sort_test_unsorted.root";
   const char *kSortedFile   = "sort_test_sorted.root";
   const char *kNameSortFile = "sort_test_namesort.root";

   void SetUp() override
   {
      GenerateSAMFile(kSamFile, kNumReads);
      std::remove(kUnsortedFile);
      std::remove(kSortedFile);
      std::remove(kNameSortFile);
      samtoramntuple(kSamFile, kUnsortedFile, true, true, true, 505, 0);
   }

   void TearDown() override
   {
      std::remove(kSamFile);
      std::remove(kUnsortedFile);
      std::remove(kSortedFile);
      std::remove(kNameSortFile);
   }
};

TEST_F(RAMSortTest, EntryCountPreserved)
{
   ASSERT_EQ(ramsortntuple(kUnsortedFile, kSortedFile), 0);
   auto readerIn  = ROOT::RNTupleReader::Open("RAM", kUnsortedFile);
   auto readerOut = ROOT::RNTupleReader::Open("RAM", kSortedFile);
   ASSERT_NE(readerIn,  nullptr);
   ASSERT_NE(readerOut, nullptr);
   EXPECT_EQ(readerIn->GetNEntries(), readerOut->GetNEntries());
}

TEST_F(RAMSortTest, CoordinateSortOrder)
{
   ASSERT_EQ(ramsortntuple(kUnsortedFile, kSortedFile), 0);
   auto reader    = ROOT::RNTupleReader::Open("RAM", kSortedFile);
   ASSERT_NE(reader, nullptr);
   auto viewRefId = reader->GetView<int32_t>("record.refid");
   auto viewPos   = reader->GetView<int32_t>("record.pos");
   int32_t prevRefId = -1, prevPos = -1;
   for (uint64_t i = 0; i < reader->GetNEntries(); ++i) {
      int32_t refid = viewRefId(i);
      int32_t pos   = viewPos(i);
      if (refid == prevRefId) {
         EXPECT_GE(pos, prevPos) << "pos out of order at entry " << i;
      } else {
         EXPECT_GE(refid, prevRefId) << "refid out of order at entry " << i;
      }
      prevRefId = refid;
      prevPos   = pos;
   }
}

TEST_F(RAMSortTest, NameSortOrder)
{
   ASSERT_EQ(ramsortntuple(kUnsortedFile, kNameSortFile, true), 0);
   auto reader    = ROOT::RNTupleReader::Open("RAM", kNameSortFile);
   ASSERT_NE(reader, nullptr);
   auto viewQname = reader->GetView<std::string>("record.qname");
   std::string prev = "";
   for (uint64_t i = 0; i < reader->GetNEntries(); ++i) {
      std::string qname = viewQname(i);
      EXPECT_GE(qname, prev) << "qname out of order at entry " << i;
      prev = qname;
   }
}

TEST_F(RAMSortTest, IdempotentSort)
{
   ASSERT_EQ(ramsortntuple(kUnsortedFile, kSortedFile), 0);
   const char *doubleSorted = "sort_test_double.root";
   ASSERT_EQ(ramsortntuple(kSortedFile, doubleSorted), 0);
   auto r1 = ROOT::RNTupleReader::Open("RAM", kSortedFile);
   auto r2 = ROOT::RNTupleReader::Open("RAM", doubleSorted);
   ASSERT_NE(r1, nullptr);
   ASSERT_NE(r2, nullptr);
   EXPECT_EQ(r1->GetNEntries(), r2->GetNEntries());
   auto v1refid = r1->GetView<int32_t>("record.refid");
   auto v2refid = r2->GetView<int32_t>("record.refid");
   auto v1pos   = r1->GetView<int32_t>("record.pos");
   auto v2pos   = r2->GetView<int32_t>("record.pos");
   for (uint64_t i = 0; i < r1->GetNEntries(); ++i) {
      EXPECT_EQ(v1refid(i), v2refid(i));
      EXPECT_EQ(v1pos(i),   v2pos(i));
   }
   std::remove(doubleSorted);
}

TEST(RAMSortEdgeCases, MissingInputFileReturnsError)
{
   int ret = ramsortntuple("nonexistent.root", "/tmp/out.root");
   EXPECT_NE(ret, 0);
}

} // namespace
