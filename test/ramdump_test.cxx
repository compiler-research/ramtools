// Runs the built samtoramntuple, ramdump and ramstats on a small SAM and checks
// that the dump reproduces the input and that its counts follow samtools' rules.
#include <gtest/gtest.h>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>

namespace {

const char *const kSamFile = "ramdump_test.sam";
const char *const kRamFile = "ramdump_test.root";
const char *const kOutFile = "ramdump_test.out";

std::string Header()
{
   std::string h{};
   h += "@HD\tVN:1.6\tSO:coordinate\n";
   h += "@SQ\tSN:chr1\tLN:1000\n";
   h += "@SQ\tSN:chr2\tLN:1000\n";
   h += "@PG\tID:test\tPN:test\n";
   return h;
}

std::string Records()
{
   std::string r{};
   r += "r1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tXA:Z:a,b\n";
   r += "r2\t16\tchr1\t150\t30\t2S50M\t=\t100\t-50\t";
   r += std::string(52, 'A');
   r += "\t";
   r += std::string(52, 'F');
   r += "\n";
   r += "r3\t256\tchr1\t300\t0\t4M\t*\t0\t0\t*\t*\n";
   r += "r4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
   r += "r5\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
   return r;
}

std::string Command(const std::string &options, const std::string &file, const std::string &region)
{
   std::string cmd{};
   cmd += RAMDUMP_BIN;
   cmd += " ";
   cmd += options;
   cmd += " ";
   cmd += file;
   cmd += " ";
   cmd += region;
   return cmd;
}

// What a command printed on standard output.
std::string Output(const std::string &command)
{
   const std::unique_ptr<FILE, int (*)(FILE *)> pipe(popen(command.c_str(), "r"), pclose);
   std::string out{};
   std::array<char, 4096> buf{};
   while (pipe && fgets(buf.data(), buf.size(), pipe.get()))
      out += buf.data();
   return out;
}

// ramdump <options> <file> <region>, returning what it printed.
std::string Dump(const std::string &options, const std::string &region = "")
{
   return Output(Command(options, kRamFile, region));
}

// Runs a command with its output discarded; 0 on success.
int Quiet(const std::string &command)
{
   return std::system((command + " >/dev/null 2>&1").c_str()); // NOLINT(bugprone-command-processor)
}

// ramdump <options> <file> <region> with its output discarded; 0 on success.
int Status(const std::string &options, const std::string &file, const std::string &region = "")
{
   return Quiet(Command(options, file, region));
}

std::string ReadFile(const char *path)
{
   std::ifstream in(path);
   return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

class RamdumpTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void SetUp() override
   {
      std::ofstream(kSamFile) << Header() << Records();
      std::string cmd{};
      cmd += SAMTORAMNTUPLE_BIN;
      cmd += " ";
      cmd += kSamFile;
      cmd += " ";
      cmd += kRamFile;
      cmd += " >/dev/null 2>&1";
      ASSERT_EQ(std::system(cmd.c_str()), 0);
   }

   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override
   {
      std::remove(kSamFile);
      std::remove(kRamFile);
      std::remove(kOutFile);
   }
};

TEST_F(RamdumpTest, DumpWithHeaderReproducesTheInput)
{
   EXPECT_EQ(Dump("-h"), Header() + Records());
}

// Standard tags are stored in columns of their own; the dump has to give back
// every tag as it was, in its place, including those that stay text.
TEST_F(RamdumpTest, TagsComeBackAsTheyWere)
{
   const char *sam = "ramdump_test_tags.sam";
   const char *ram = "ramdump_test_tags.root";
   std::string records{};
   records += "t1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tMC:Z:4M\tPG:Z:MarkDuplicates\tMQ:i:60\tAS:i:4\t"
              "XS:i:0\tMD:Z:4\tNM:i:0\tRG:Z:lane1\n";
   records += "t2\t0\tchr1\t110\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:lane2\tSA:Z:chr2,500,+,2M2S,60,0;\tAS:i:-3\t"
              "pa:f:0.5\tXA:Z:chr1,+300,4M,0;\n";
   records += "t3\t0\tchr1\t120\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
   records += "t4\t0\tchr1\t130\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:007\tXS:i:4294967295\tNM:Z:x\tRG:i:5\tMD:i:4\n";
   records += "t5\t0\tchr1\t140\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:1\tNM:i:2\tRG:Z:\tRG:Z:lane1\tMD:Z:4\tMD:Z:3A0\n";
   std::ofstream(sam) << Header() << records;
   ASSERT_EQ(Quiet(std::string(SAMTORAMNTUPLE_BIN) + " " + sam + " " + ram + " -threads 2"), 0);
   EXPECT_EQ(Output(Command("-h", ram, "")), Header() + records);
   std::remove(sam);
   std::remove(ram);
}

// PNEXT and TLEN are stored against pos and the CIGAR; every kind of mate
// comes back as it was, a TLEN no aligner would write among them.
TEST_F(RamdumpTest, MateFieldsComeBackAsTheyWere)
{
   const char *sam = "ramdump_test_mates.sam";
   const char *ram = "ramdump_test_mates.root";
   std::string records{};
   records += "m1\t99\tchr1\t100\t60\t4M\t=\t300\t204\tACGT\tIIII\n";
   records += "m1\t147\tchr1\t300\t60\t4M\t=\t100\t-204\tACGT\tIIII\n";
   records += "m2\t99\tchr1\t400\t60\t2S2M\t=\t500\t103\tACGT\tIIII\n";
   records += "m3\t99\tchr1\t600\t60\t4M\tchr1\t700\t104\tACGT\tIIII\n";
   records += "m4\t65\tchr1\t800\t60\t4M\tchr2\t10\t0\tACGT\tIIII\n";
   records += "m5\t73\tchr1\t900\t60\t4M\t=\t900\t0\tACGT\tIIII\n";
   records += "m5\t133\tchr1\t900\t0\t*\t=\t900\t0\tACGT\tIIII\n";
   records += "m6\t0\tchr2\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
   records += "m7\t99\tchr2\t100\t60\t4M\t=\t200\t-5000\tACGT\tIIII\n";
   records += "m8\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
   std::ofstream(sam) << Header() << records;
   ASSERT_EQ(Quiet(std::string(SAMTORAMNTUPLE_BIN) + " " + sam + " " + ram + " -threads 2"), 0);
   EXPECT_EQ(Output(Command("-h", ram, "")), Header() + records);
   std::remove(sam);
   std::remove(ram);
}

TEST_F(RamdumpTest, HeaderOnly)
{
   EXPECT_EQ(Dump("-H"), Header());
}

TEST_F(RamdumpTest, CountsFollowSamtools)
{
   EXPECT_EQ(Dump("-c"), "5\n");
   EXPECT_EQ(Dump("-c", "chr1"), "3\n") << "secondary reads count, unmapped do not";
   EXPECT_EQ(Dump("-c", "chr1:150-200"), "1\n");
   EXPECT_EQ(Dump("-c", "chr1:103-149"), "1\n") << "r1 ends on 103";
   EXPECT_EQ(Dump("-c -F 0x904", "chr1"), "2\n") << "-F drops the secondary read";
   EXPECT_EQ(Dump("-c -f 0x10"), "1\n") << "-f keeps only the reverse read";
}

TEST_F(RamdumpTest, RegionAndOutputFile)
{
   ASSERT_EQ(Status("-o ramdump_test.out", kRamFile, "chr2"), 0);
   EXPECT_EQ(ReadFile(kOutFile), "r5\t0\tchr2\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
}

// -threads only changes how pages are compressed; the records must not change.
TEST_F(RamdumpTest, ThreadedConversionWritesTheSameRecords)
{
   const char *threaded = "ramdump_test_threads.root";
   std::string cmd{};
   cmd += SAMTORAMNTUPLE_BIN;
   cmd += " ";
   cmd += kSamFile;
   cmd += " ";
   cmd += threaded;
   cmd += " -threads 2 >/dev/null 2>&1";
   ASSERT_EQ(std::system(cmd.c_str()), 0);
   ASSERT_EQ(Status("-h -o ramdump_test.out", threaded), 0);
   EXPECT_EQ(ReadFile(kOutFile), Header() + Records());
   std::remove(threaded);
}

// The expected text is what samtools 1.18 idxstats prints for the same SAM.
TEST_F(RamdumpTest, IdxstatsFollowsSamtools)
{
   const std::string ramstats = RAMSTATS_BIN;
   EXPECT_EQ(Output(ramstats + " idxstats " + kRamFile), "chr1\t1000\t3\t0\nchr2\t1000\t1\t0\n*\t0\t0\t1\n");
   EXPECT_NE(Quiet(ramstats), 0);
   EXPECT_NE(Quiet(ramstats + " depth " + kRamFile), 0);
   EXPECT_NE(Quiet(ramstats + " idxstats missing.root"), 0);
}

// Every case bam_stat.c tells apart: proper and improper pairs, a singleton,
// mates on another reference with MAPQ above and below 5, an RNEXT that names
// the read's own reference, duplicates, secondary and supplementary records,
// and QC-failed records without pairs, for which samtools prints N/A. The
// expected texts are what samtools 1.18 flagstat prints for this SAM; bam_stat.c
// is the same up to 1.24.
TEST_F(RamdumpTest, FlagstatFollowsSamtools)
{
   const char *sam = "ramdump_test_flagstat.sam";
   const char *ram = "ramdump_test_flagstat.root";
   std::ofstream(sam) << "@HD\tVN:1.6\tSO:unsorted\n"
                         "@SQ\tSN:chr1\tLN:1000\n"
                         "@SQ\tSN:chr2\tLN:1000\n"
                         "p1\t99\tchr1\t100\t60\t4M\t=\t300\t204\tACGT\tIIII\n"
                         "p1\t147\tchr1\t300\t60\t4M\t=\t100\t-204\tACGT\tIIII\n"
                         "p2\t65\tchr1\t500\t60\t4M\tchr2\t10\t0\tACGT\tIIII\n"
                         "p2\t129\tchr2\t10\t3\t4M\tchr1\t500\t0\tACGT\tIIII\n"
                         "p3\t73\tchr1\t600\t60\t4M\t=\t600\t0\tACGT\tIIII\n"
                         "p3\t133\tchr1\t600\t0\t*\t=\t600\t0\tACGT\tIIII\n"
                         "p4\t67\tchr1\t700\t60\t4M\tchr1\t800\t104\tACGT\tIIII\n"
                         "p4\t131\tchr1\t800\t60\t4M\tchr1\t700\t-104\tACGT\tIIII\n"
                         "d1\t1024\tchr2\t50\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                         "s1\t256\tchr1\t100\t0\t4M\t*\t0\t0\t*\t*\n"
                         "u1\t2048\tchr2\t60\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                         "u2\t3072\tchr2\t70\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
                         "x1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
                         "q1\t516\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
                         "q2\t512\tchr2\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
   ASSERT_EQ(Quiet(std::string(SAMTORAMNTUPLE_BIN) + " " + sam + " " + ram), 0);
   const std::string flagstat = std::string(RAMSTATS_BIN) + " flagstat ";

   EXPECT_EQ(Output(flagstat + ram), "13 + 2 in total (QC-passed reads + QC-failed reads)\n"
                                     "10 + 2 primary\n"
                                     "1 + 0 secondary\n"
                                     "2 + 0 supplementary\n"
                                     "2 + 0 duplicates\n"
                                     "1 + 0 primary duplicates\n"
                                     "11 + 1 mapped (84.62% : 50.00%)\n"
                                     "8 + 1 primary mapped (80.00% : 50.00%)\n"
                                     "8 + 0 paired in sequencing\n"
                                     "4 + 0 read1\n"
                                     "4 + 0 read2\n"
                                     "4 + 0 properly paired (50.00% : N/A)\n"
                                     "6 + 0 with itself and mate mapped\n"
                                     "1 + 0 singletons (12.50% : N/A)\n"
                                     "2 + 0 with mate mapped to a different chr\n"
                                     "1 + 0 with mate mapped to a different chr (mapQ>=5)\n");

   EXPECT_EQ(Output(flagstat + "-O tsv " + ram), "13\t2\ttotal (QC-passed reads + QC-failed reads)\n"
                                                 "10\t2\tprimary\n"
                                                 "1\t0\tsecondary\n"
                                                 "2\t0\tsupplementary\n"
                                                 "2\t0\tduplicates\n"
                                                 "1\t0\tprimary duplicates\n"
                                                 "11\t1\tmapped\n"
                                                 "84.62%\t50.00%\tmapped %\n"
                                                 "8\t1\tprimary mapped\n"
                                                 "80.00%\t50.00%\tprimary mapped %\n"
                                                 "8\t0\tpaired in sequencing\n"
                                                 "4\t0\tread1\n"
                                                 "4\t0\tread2\n"
                                                 "4\t0\tproperly paired\n"
                                                 "50.00%\tN/A\tproperly paired %\n"
                                                 "6\t0\twith itself and mate mapped\n"
                                                 "1\t0\tsingletons\n"
                                                 "12.50%\tN/A\tsingletons %\n"
                                                 "2\t0\twith mate mapped to a different chr\n"
                                                 "1\t0\twith mate mapped to a different chr (mapQ>=5)\n");

   EXPECT_EQ(Output(flagstat + "-O json " + ram), "{\n"
                                                  " \"QC-passed reads\": { \n"
                                                  "  \"total\": 13, \n"
                                                  "  \"primary\": 10, \n"
                                                  "  \"secondary\": 1, \n"
                                                  "  \"supplementary\": 2, \n"
                                                  "  \"duplicates\": 2, \n"
                                                  "  \"primary duplicates\": 1, \n"
                                                  "  \"mapped\": 11, \n"
                                                  "  \"mapped %\": 84.62, \n"
                                                  "  \"primary mapped\": 8, \n"
                                                  "  \"primary mapped %\": 80.00, \n"
                                                  "  \"paired in sequencing\": 8, \n"
                                                  "  \"read1\": 4, \n"
                                                  "  \"read2\": 4, \n"
                                                  "  \"properly paired\": 4, \n"
                                                  "  \"properly paired %\": 50.00, \n"
                                                  "  \"with itself and mate mapped\": 6, \n"
                                                  "  \"singletons\": 1, \n"
                                                  "  \"singletons %\": 12.50, \n"
                                                  "  \"with mate mapped to a different chr\": 2, \n"
                                                  "  \"with mate mapped to a different chr (mapQ >= 5)\": 1 \n"
                                                  " },\n"
                                                  " \"QC-failed reads\": { \n"
                                                  "  \"total\": 2, \n"
                                                  "  \"primary\": 2, \n"
                                                  "  \"secondary\": 0, \n"
                                                  "  \"supplementary\": 0, \n"
                                                  "  \"duplicates\": 0, \n"
                                                  "  \"primary duplicates\": 0, \n"
                                                  "  \"mapped\": 1, \n"
                                                  "  \"mapped %\": 50.00, \n"
                                                  "  \"primary mapped\": 1, \n"
                                                  "  \"primary mapped %\": 50.00, \n"
                                                  "  \"paired in sequencing\": 0, \n"
                                                  "  \"read1\": 0, \n"
                                                  "  \"read2\": 0, \n"
                                                  "  \"properly paired\": 0, \n"
                                                  "  \"properly paired %\": null, \n"
                                                  "  \"with itself and mate mapped\": 0, \n"
                                                  "  \"singletons\": 0, \n"
                                                  "  \"singletons %\": null, \n"
                                                  "  \"with mate mapped to a different chr\": 0, \n"
                                                  "  \"with mate mapped to a different chr (mapQ >= 5)\": 0 \n"
                                                  " }\n"
                                                  "}\n");

   EXPECT_NE(Quiet(flagstat + "-O xml " + ram), 0);
   EXPECT_NE(Quiet(flagstat), 0);
   EXPECT_NE(Quiet(flagstat + "missing.root"), 0);
   std::remove(sam);
   std::remove(ram);
}

TEST_F(RamdumpTest, RejectsBadArguments)
{
   EXPECT_NE(Status("-f abc", kRamFile), 0);
   EXPECT_NE(Status("-x", kRamFile), 0);
   EXPECT_NE(Status("", "missing.root"), 0);
   EXPECT_EQ(Status("", kRamFile, "chr9"), 0) << "an unknown reference is empty, not an error";
}

} // namespace
