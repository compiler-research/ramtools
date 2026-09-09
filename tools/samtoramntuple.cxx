#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

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

int main(int argc, char* argv[]) {
    if (argc < 2) {
       std::cout << "Usage: " << argv[0] << " <input.sam> [output]\n";
       std::cout << "Options:\n";
       std::cout << "  -split       Split by chromosome\n";
       std::cout << "  -noindex     Disable indexing\n";
       std::cout << "  -illumina    Use Illumina quality binning\n";
       std::cout << "  -dropqual    Drop quality scores\n";
       std::cout << "  -compression N  ROOT compression code, algorithm*100+level (default 505, ZSTD level 5)\n";
       return 1;
    }
    
    const char* input = argv[1];
    const char* output = nullptr;

    bool do_split = false;
    bool do_index = true;
    uint32_t quality_mode = RAMNTupleRecord::kPhred33;
    int compression = 505;
    bool want_compression = false;

    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (want_compression) {
           if (!ParseCompression(arg, compression)) {
              std::cerr << "invalid -compression value '" << arg << "'\n";
              return 1;
           }
           want_compression = false;
        } else if (arg == "-split") {
           do_split = true;
        } else if (arg == "-noindex") {
           do_index = false;
        } else if (arg == "-illumina") {
           quality_mode = RAMNTupleRecord::kIlluminaBinning;
        } else if (arg == "-dropqual") {
           quality_mode = RAMNTupleRecord::kDrop;
        } else if (arg == "-compression") {
           want_compression = true;
        } else if (arg[0] != '-') {
           output = argv[i];
        }
    }
    if (want_compression) {
       std::cerr << "-compression needs a value\n";
       return 1;
    }

    std::string outfile;
    if (!output) {
       outfile = std::string(input);
       size_t pos = outfile.rfind(".sam");
       if (pos != std::string::npos) {
          outfile.erase(pos);
       }
       output = outfile.c_str();
    }

    try {
       if (do_split) {
          samtoramntuple_split_by_chromosome(input, output, compression, quality_mode);
       } else {
          std::string ramfile = std::string(output);
          if (ramfile.find(".root") == std::string::npos && ramfile.find(".ram") == std::string::npos) {
             ramfile += ".ram";
          }
          samtoramntuple(input, ramfile.c_str(), do_index, true, true, compression, quality_mode);
       }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

