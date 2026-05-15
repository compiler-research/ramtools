#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <iostream>
#include <string>
#include <cstring>
#include <Compression.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
       std::cout << "Usage: " << argv[0] << " <input.sam> [output]\n";
       std::cout << "Options:\n";
       std::cout << "  -split            Split by chromosome\n";
       std::cout << "  -noindex          Disable indexing\n";
       std::cout << "  -illumina         Use Illumina quality binning\n";
       std::cout << "  -dropqual         Drop quality scores\n";
       std::cout << "  --compression A   Compression algorithm: zstd (default), lz4, lzma, zlib\n";
       std::cout << "  --level N         Compression level 1-19 (default 5)\n";
       return 1;
    }

    const char* input = argv[1];
    const char* output = nullptr;

    bool do_split = false;
    bool do_index = true;
    uint32_t quality_mode = RAMNTupleRecord::kPhred33;
    int compression_algo  = ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD;
    int compression_level = 5;

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
        } else if (arg == "--compression" && i + 1 < argc) {
           std::string val = argv[++i];
           if      (val == "zstd" || val == "ZSTD") compression_algo = ROOT::RCompressionSetting::EAlgorithm::EValues::kZSTD;
           else if (val == "lz4"  || val == "LZ4")  compression_algo = ROOT::RCompressionSetting::EAlgorithm::EValues::kLZ4;
           else if (val == "lzma" || val == "LZMA") compression_algo = ROOT::RCompressionSetting::EAlgorithm::EValues::kLZMA;
           else if (val == "zlib" || val == "ZLIB") compression_algo = ROOT::RCompressionSetting::EAlgorithm::EValues::kZLIB;
           else {
              std::cerr << "Unknown compression algorithm: " << val
                        << ". Valid: zstd, lz4, lzma, zlib\n";
              return 1;
           }
        } else if (arg == "--level" && i + 1 < argc) {
           compression_level = std::atoi(argv[++i]);
           if (compression_level < 1 || compression_level > 19) {
              std::cerr << "Compression level must be 1-19\n";
              return 1;
           }
        } else if (arg[0] != '-') {
           output = argv[i];
        }
    }

    // ROOT compression setting format: level * 100 + algorithm_id
    int compression_setting = compression_level * 100 + compression_algo;

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
          samtoramntuple_split_by_chromosome(input, output, compression_setting, quality_mode, 4);
       } else {
          std::string ramfile = std::string(output);
          if (ramfile.find(".root") == std::string::npos && ramfile.find(".ram") == std::string::npos) {
             ramfile += ".ram";
          }
          samtoramntuple(input, ramfile.c_str(), do_index, true, true, compression_setting, quality_mode);
       }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}