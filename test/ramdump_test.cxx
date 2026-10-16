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

// ramdump <options> <file> <region> with its output discarded; 0 on success.
int Status(const std::string &options, const std::string &file, const std::string &region = "")
{
   return std::system((Command(options, file, region) + " >/dev/null 2>&1").c_str());
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
   EXPECT_NE(std::system((ramstats + " >/dev/null 2>&1").c_str()), 0);
   EXPECT_NE(std::system((ramstats + " flagstat " + kRamFile + " >/dev/null 2>&1").c_str()), 0);
   EXPECT_NE(std::system((ramstats + " idxstats missing.root >/dev/null 2>&1").c_str()), 0);
}

TEST_F(RamdumpTest, RejectsBadArguments)
{
   EXPECT_NE(Status("-f abc", kRamFile), 0);
   EXPECT_NE(Status("-x", kRamFile), 0);
   EXPECT_NE(Status("", "missing.root"), 0);
   EXPECT_EQ(Status("", kRamFile, "chr9"), 0) << "an unknown reference is empty, not an error";
}

} // namespace
