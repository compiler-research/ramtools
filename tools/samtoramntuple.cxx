#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <iostream>
#include <string>
#include <cstring>

int main(int argc, char* argv[]) {
    if (argc < 2) {
       std::cout << "Usage: " << argv[0] << " <input.sam> [output] [zlib|lzma|lz4|zstd]\n";
       std::cout << "Options:\n";
       std::cout << "  -split                Split by chromosome\n";
       std::cout << "  -noindex              Disable indexing\n";
       std::cout << "  -illumina             Use Illumina quality binning\n";
       std::cout << "  -dropqual             Drop quality scores\n";
       std::cout << "  zlib|lzma|lz4|zstd    RNTuple compression codec (default: zstd, level 5)\n";
       return 1;
    }

    const char* input = argv[1];
    const char* output = nullptr;

    bool do_split = false;
    bool do_index = true;
    uint32_t quality_mode = RAMNTupleRecord::kPhred33;
    // Default keeps existing behavior: ZSTD level 5 = algorithm 5 * 100 + 5.
    int compression = 505;
    const char* codec_name = "zstd";

    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-split") {
           do_split = true;
        } else if (arg == "-noindex") {
           do_index = false;
        } else if (arg == "-illumina") {
           quality_mode = RAMNTupleRecord::kIlluminaBinning;
        } else if (arg == "-dropqual") {
           quality_mode = RAMNTupleRecord::kDrop;
        } else if (arg == "zlib") {
           compression = 101; codec_name = "zlib";   // ZLIB level 1 (fast default)
        } else if (arg == "lzma") {
           compression = 207; codec_name = "lzma";   // LZMA level 7 (ROOT default)
        } else if (arg == "lz4") {
           compression = 404; codec_name = "lz4";    // LZ4 level 4
        } else if (arg == "zstd") {
           compression = 505; codec_name = "zstd";   // ZSTD level 5 (existing default)
        } else if (arg[0] != '-') {
           output = argv[i];
        }
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

    std::cout << "samtoramntuple: " << input << " -> " << output
              << " (compression=" << codec_name << ", code=" << compression << ")\n";

    try {
       if (do_split) {
          samtoramntuple_split_by_chromosome(input, output, compression, quality_mode, 4);
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

