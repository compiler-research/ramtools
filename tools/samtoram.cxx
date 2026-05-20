#include "ramcore/SamToTTree.h"
#include "ttree/RAMRecord.h"
#include <Compression.h>
#include <cstdio>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <input.sam> [output.root] [zlib|lzma|lz4|zstd]\n", argv[0]);
        printf("  compression defaults to lzma (unchanged from previous behavior)\n");
        return 1;
    }

    const char* input = argv[1];
    const char* output = (argc > 2 && argv[2][0] != '\0') ? argv[2] : "ramexample.root";

    // Optional 3rd arg: compression algorithm. Default stays LZMA so existing
    // callers/tests are unaffected.
    int algo = ROOT::RCompressionSetting::EAlgorithm::kLZMA;
    const char* algo_name = "lzma";
    if (argc > 3) {
        std::string a = argv[3];
        if (a == "zlib") {
            algo = ROOT::RCompressionSetting::EAlgorithm::kZLIB;
            algo_name = "zlib";
        } else if (a == "lzma") {
            algo = ROOT::RCompressionSetting::EAlgorithm::kLZMA;
            algo_name = "lzma";
        } else if (a == "lz4") {
            algo = ROOT::RCompressionSetting::EAlgorithm::kLZ4;
            algo_name = "lz4";
        } else if (a == "zstd") {
            algo = ROOT::RCompressionSetting::EAlgorithm::kZSTD;
            algo_name = "zstd";
        } else {
            fprintf(stderr, "Unknown compression '%s' (use zlib|lzma|lz4|zstd)\n", argv[3]);
            return 1;
        }
    }

    printf("samtoram: %s -> %s (compression=%s)\n", input, output, algo_name);
    samtoram(input, output, true, true, true, algo, RAMRecord::kPhred33);

    return 0;
}
