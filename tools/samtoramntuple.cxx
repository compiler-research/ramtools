#include "RAMCore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"
#include <iostream>
#include <string>
#include <cstring>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.sam> [output.ram]\n";
        std::cout << "Options:\n";
        std::cout << "  -noindex     Disable indexing\n";
        std::cout << "  -illumina    Use Illumina quality binning\n";
        std::cout << "  -dropqual    Drop quality scores\n";
        return 1;
    }
    
    const char* input = argv[1];
    const char* output = nullptr;
    
    std::string outfile;
    if (argc > 2 && argv[2][0] != '-') {
        output = argv[2];
    } else {
        outfile = std::string(input);
        size_t pos = outfile.rfind(".sam");
        if (pos != std::string::npos) {
            outfile.replace(pos, 4, ".ram");
        } else {
            outfile += ".ram";
        }
        output = outfile.c_str();
    }
    
   
    bool do_index = true;
    uint32_t quality_mode = RAMNTupleRecord::kPhred33;
    
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-noindex") {
            do_index = false;
        } else if (arg == "-illumina") {
            quality_mode = RAMNTupleRecord::kIlluminaBinning;
        } else if (arg == "-dropqual") {
            quality_mode = RAMNTupleRecord::kDrop;
        }
    }
    
    try {
        samtoramntuple(input, output, do_index, true, true, 505, quality_mode);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

