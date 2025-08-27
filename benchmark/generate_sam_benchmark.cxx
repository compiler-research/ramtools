#include <benchmark/benchmark.h>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <iostream>

void GenerateSAMFile(const std::string& filename, int num_reads) {
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
    
    std::ofstream out(filename);
    
    for (const auto& [chrom, length] : chromosomes) {
        out << "@SQ\tSN:" << chrom << "\tLN:" << length << "\n";
    }
    
    std::uniform_int_distribution<> base_dist(0, 3);
    std::uniform_int_distribution<> flag_dist(0, 1);
    std::uniform_int_distribution<> chrom_dist(0, chromosomes.size() - 1);
    std::uniform_int_distribution<> qual_dist(33, 73);
    std::uniform_int_distribution<> fc_dist(10000, 99999);
    std::uniform_int_distribution<> lane_dist(1, 8);
    std::uniform_int_distribution<> tile_dist(1, 300);
    std::uniform_int_distribution<> coord_dist(1, 1000);
    std::uniform_int_distribution<> mismatch_dist(0, 9);
    
    for (int i = 0; i < num_reads; ++i) {
        auto [chrom, chrom_length] = chromosomes[chrom_dist(rng)];
        int read_length = 36;
        std::uniform_int_distribution<> pos_dist(1, std::max(1, chrom_length - read_length));
        int position = pos_dist(rng);
        
        out << "SOLEXA-1GA-2_2_FC" << fc_dist(rng) 
            << ":" << lane_dist(rng) 
            << ":" << tile_dist(rng)
            << ":" << coord_dist(rng)
            << ":" << coord_dist(rng) << "\t";
        
        out << (flag_dist(rng) * 16) << "\t"
            << chrom << "\t"
            << position << "\t"
            << 25 << "\t";
        
        out << read_length << "M\t";
        
        out << "*\t0\t0\t";
        
        std::string sequence;
        sequence.reserve(read_length);
        for (int j = 0; j < read_length; ++j) {
            sequence += bases[base_dist(rng)];
        }
        out << sequence << "\t";
        
        for (int j = 0; j < read_length; ++j) {
            out << static_cast<char>(qual_dist(rng));
        }
        
        int nm = mismatch_dist(rng) < 7 ? 0 : 1;
        out << "\tNM:i:" << nm;
        out << "\tX" << nm << ":i:1";
        out << "\tMD:Z:" << read_length;
        
        out << "\n";
    }
    
    out.close();
}

static void BM_GenerateSAM(benchmark::State& state) {
    int num_reads = state.range(0);
    for (auto _ : state) {
        GenerateSAMFile("benchmark_temp.sam", num_reads);
        std::remove("benchmark_temp.sam");
    }
    
    state.counters["reads_per_second"] = benchmark::Counter(
        num_reads, benchmark::Counter::kIsRate);
    state.counters["bytes_per_second"] = benchmark::Counter(
        num_reads * 200, benchmark::Counter::kIsRate);
}

int main(int argc, char** argv) {
    // MODE 1: File generation (used during build process and manual testing)
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--generate" && i + 1 < argc) {
            std::string filename = argv[i + 1];
            int num_reads = (i + 2 < argc) ? std::stoi(argv[i + 2]) : 100;
            GenerateSAMFile(filename, num_reads);
            std::cout << "Generated SAM file: " << filename << " with " << num_reads << " reads" << std::endl;
            return 0;
        }
    }
    
    // MODE 2: Benchmark with configurable range
    // Parse the external parameter for maximum number of entries to test
    int min_reads = 100;
    int max_reads = 100000;  // Default maximum if not specified
    
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--max_reads" && i + 1 < argc) {
            max_reads = std::stoi(argv[i + 1]);
            for (int j = i; j < argc - 2; ++j) {
                argv[j] = argv[j + 2];
            }
            argc -= 2;
            break;
        }
    }
    
    std::cout << "Benchmarking SAM generation from " << min_reads << " to " << max_reads << " reads" << std::endl;
    std::cout << "Scaling factor: 10x per step" << std::endl;
    
    // Register benchmark with dynamic range based on external parameter
    ::benchmark::RegisterBenchmark("BM_GenerateSAM", BM_GenerateSAM)
        ->RangeMultiplier(10)  // Test 100, 1000, 10000, etc.
        ->Range(min_reads, max_reads)
        ->Unit(benchmark::kMicrosecond);
    
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    return 0;
}

