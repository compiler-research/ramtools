#include "RAMCore/SamToTTree.h"
#include "ttree/RAMRecord.h"
#include <Compression.h>
#include <cstdio>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s <input.sam> [output.root]\n", argv[0]);
        return 1;
    }
    
    const char* input = argv[1];
    const char* output = argc > 2 ? argv[2] : "ramexample.root";
    
    samtoram(input, output, true, true, true, 
             ROOT::RCompressionSetting::EAlgorithm::kLZMA, 
             RAMRecord::kPhred33);
    
    return 0;
}

