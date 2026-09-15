#include "ramcore/BamtoNTuple.h"

#include "rntuple/RAMNTupleRecord.h"

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

} // namespace

int main(int argc, char *argv[])
{
   if (argc < 2) {
      std::cout << "Usage: " << argv[0] << " <input.bam> [output]\n"
                << "Options:\n"
                << "  -illumina    Use Illumina quality binning\n"
                << "  -dropqual    Drop quality scores\n"
                << "  -compression N  ROOT compression code, algorithm*100+level (default 505, ZSTD level 5)\n";
      return 1;
   }

   const char *input = argv[1];
   const char *output = nullptr;

   uint32_t quality_mode = RAMNTupleRecord::kPhred33;
   int compression = 505;
   bool want_compression = false;

   for (int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if (want_compression) {
         if (!ParseCompression(arg, compression)) {
            std::cerr << "invalid -compression value '" << arg << "'\n";
            return 1;
         }
         want_compression = false;
      } else if (arg == "-illumina" || arg == "-dropqual")
         quality_mode = (arg == "-illumina") ? RAMNTupleRecord::kIlluminaBinning : RAMNTupleRecord::kDrop;
      else if (arg == "-compression")
         want_compression = true;
      else if (arg[0] != '-')
         output = argv[i];
   }
   if (want_compression) {
      std::cerr << "-compression needs a value\n";
      return 1;
   }

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