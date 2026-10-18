#include "ramcore/BamtoNTuple.h"

#include "rntuple/RAMNTupleRecord.h"

#include <TROOT.h>
#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

// ROOT compression code: algorithm*100+level, with algorithm 1 ZLIB, 2 LZMA,
// 4 LZ4 or 5 ZSTD and level 1 to 9; 0 means no compression.
bool ParseCompression(const std::string &text, int &code)
{
   char *end = nullptr;
   const long value = std::strtol(text.c_str(), &end, 10);
   if (text.empty() || *end != '\0' || value < 0)
      return false;
   const long algorithm = value / 100;
   const long level = value % 100;
   const bool known = algorithm == 1 || algorithm == 2 || algorithm == 4 || algorithm == 5;
   if (value != 0 && (!known || level < 1 || level > 9))
      return false;
   code = static_cast<int>(value);
   return true;
}

// Thread count: a positive decimal integer and nothing else, so "3x" is an
// error rather than 3.
bool ParseThreads(const std::string &text, int &threads)
{
   char *end = nullptr;
   const long value = std::strtol(text.c_str(), &end, 10);
   if (text.empty() || !std::isdigit(static_cast<unsigned char>(text[0])) || *end != '\0' || value < 1 ||
       value > INT_MAX)
      return false;
   threads = static_cast<int>(value);
   return true;
}

} // namespace

int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cout << "Usage: " << argv[0] << " <input.bam> [output]\n"
                << "Options:\n"
                << "  -illumina    Use Illumina quality binning\n"
                << "  -dropqual    Drop quality scores\n"
                << "  -fqzcomp     Compress quality scores with fqzcomp instead of rANS: may be\n"
                << "               smaller for unbinned qualities, three times the CPU to read\n"
                << "  -compression N  ROOT compression code, algorithm*100+level (default 505, ZSTD level 5)\n"
                << "  -threads N   compress pages on N threads (default 1)\n";
      return 1;
   }

   const char *input = argv[1];
   const char *output = nullptr;

   uint32_t quality_mode = RAMNTupleRecord::kPhred33;
   bool fqzcomp = false;
   int compression = 505;
   bool want_compression = false;
   int threads = 1;
   bool want_threads = false;

   for (int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if (want_compression) {
         if (!ParseCompression(arg, compression)) {
            std::cerr << "invalid -compression value '" << arg << "'\n";
            return 1;
         }
         want_compression = false;
      } else if (want_threads) {
         if (!ParseThreads(arg, threads)) {
            std::cerr << "invalid -threads value '" << arg << "'\n";
            return 1;
         }
         want_threads = false;
      } else if (arg == "-threads")
         want_threads = true;
      else if (arg == "-illumina" || arg == "-dropqual")
         quality_mode = (arg == "-illumina") ? RAMNTupleRecord::kIlluminaBinning : RAMNTupleRecord::kDrop;
      else if (arg == "-fqzcomp")
         fqzcomp = true;
      else if (arg == "-compression")
         want_compression = true;
      else if (arg[0] != '-')
         output = argv[i];
   }
   if (want_compression) {
      std::cerr << "-compression needs a value\n";
      return 1;
   }
   if (want_threads) {
      std::cerr << "-threads needs a value\n";
      return 1;
   }
   if (fqzcomp) {
      if (quality_mode != RAMNTupleRecord::kPhred33) {
         std::cerr << "-fqzcomp compresses lossless qualities; drop -illumina or -dropqual\n";
         return 1;
      }
      quality_mode |= RAMNTupleRecord::kFqzcomp;
   }
   // The writer compresses pages on ROOT's thread pool once implicit
   // multithreading is on; the BAM reader already uses htslib threads.
   if (threads > 1)
      ROOT::EnableImplicitMT(threads);

   std::string outfile;
   if (output == nullptr) {
      outfile = input;
      const auto pos = outfile.rfind(".bam");
      if (pos != std::string::npos)
         outfile.erase(pos);
      output = outfile.c_str();
   }

   std::string ramfile = output;
   if (ramfile.find(".root") == std::string::npos && ramfile.find(".ram") == std::string::npos)
      ramfile += ".ram";

   bamtoramntuple(input, ramfile.c_str(),
                  /*split=*/false, /*cache=*/true,
                  /*compression_algorithm=*/compression, /*quality_policy=*/quality_mode);

   return 0;
}