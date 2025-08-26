#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <utility>

int main(int argc, char* argv[]) {
    std::string filename = argc > 1 ? argv[1] : "samexample.sam";
    int num_reads = argc > 2 ? std::stoi(argv[2]) : 100;
    
    std::mt19937 rng(42);
    
    std::vector<std::pair<std::string, int>> chromosomes = {
        {"chr1", 249250621}, {"chr2", 243199373}, {"chr3", 198022430},
        {"chr4", 191154276}, {"chr5", 180915260}, {"chr6", 171115067},
        {"chr7", 159138663}, {"chr8", 146364022}, {"chr9", 141213431},
        {"chr10", 135534747}, {"chr11", 135006516}, {"chr12", 133851895},
        {"chr13", 115169878}, {"chr14", 107349540}, {"chr15", 102531392},
        {"chr16", 90354753}, {"chr17", 81195210}, {"chr18", 78077248},
        {"chr19", 59128983}, {"chr20", 63025520}, {"chr21", 48129895},
        {"chr22", 51304566}, {"chrM", 16571}, {"chrX", 155270560}
    };
    
    const char bases[4] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<> base_dist(0, 3);
    std::uniform_int_distribution<> flag_dist(0, 1);
    std::uniform_int_distribution<> chrom_dist(0, chromosomes.size() - 1);
    std::uniform_int_distribution<> qual_dist(33, 73);
    std::uniform_int_distribution<> fc_dist(10000, 99999);
    std::uniform_int_distribution<> lane_dist(1, 8);
    std::uniform_int_distribution<> tile_dist(1, 300);
    std::uniform_int_distribution<> coord_dist(1, 1000);
    std::uniform_int_distribution<> mismatch_dist(0, 9);
    
    std::ofstream out(filename);
    
    // Write header
    for (const auto& [chrom, length] : chromosomes) {
        out << "@SQ\tSN:" << chrom << "\tLN:" << length << "\n";
    }
    
    // Generate reads
    for (int i = 0; i < num_reads; ++i) {
        auto [chrom, chrom_length] = chromosomes[chrom_dist(rng)];
        int read_length = 36;
        std::uniform_int_distribution<> pos_dist(1, std::max(1, chrom_length - read_length));
        int position = pos_dist(rng);
        
        // Generate read name
        out << "SOLEXA-1GA-2_2_FC" << fc_dist(rng) 
            << ":" << lane_dist(rng) 
            << ":" << tile_dist(rng)
            << ":" << coord_dist(rng)
            << ":" << coord_dist(rng) << "\t";
        
        // Flag, chrom, position, mapq
        out << (flag_dist(rng) * 16) << "\t"
            << chrom << "\t"
            << position << "\t"
            << 25 << "\t";
        
        // CIGAR
        out << read_length << "M\t";
        
        // Mate info
        out << "*\t0\t0\t";
        
        // Sequence
        std::string sequence;
        for (int j = 0; j < read_length; ++j) {
            sequence += bases[base_dist(rng)];
        }
        out << sequence << "\t";
        
        // Quality
        for (int j = 0; j < read_length; ++j) {
            out << static_cast<char>(qual_dist(rng));
        }
        
        // Optional fields
        int nm = mismatch_dist(rng) < 7 ? 0 : 1;
        out << "\tNM:i:" << nm;
        out << "\tX" << nm << ":i:1";
        out << "\tMD:Z:" << read_length;
        
        out << "\n";
    }
    
    out.close();
    return 0;
}

