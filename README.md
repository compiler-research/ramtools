# RAMTools - ROOT Alignment/Map Format Tools

RAMTools provides efficient tools for converting SAM files to ROOT's modern, columnar RNTuple format (RAM - ROOT Alignment/Map) and working with genomic alignment data.

## Features

- High-performance SAM to RAM conversion
- Chromosome-based splitting for parallel processing
- Region-based querying capabilities

## Requirements

- ROOT 6.38+
- C++17 compatible compiler
- CMake 3.16+

## Quick Start
```bash
# 1. Build the tools
mkdir build && cd build
cmake ..
make -j$(nproc)

# 2. Convert a SAM file to the RAM format
./tools/samtoramntuple input.sam output.root

# 3. Query a specific region from the command line
./tools/ramntupleview output.root "chr1:15700-15800"
```

## Command-Line Tools

The primary way to interact with RAMTools is through these command-line executables.

### SAM to RAM Conversion

Convert a standard SAM file into the optimized RNTuple-based RAM format.
```bash
# Basic conversion
./tools/samtoramntuple input.sam output.root

# Split by chromosome for parallel processing
# (Creates output_chr1.root, output_chr2.root, etc.)
./tools/samtoramntuple input.sam output -split
```

Options: `-noindex` skips the region index, `-illumina` stores 8-level binned
quality scores, `-dropqual` stores none, `-compression N` sets the ROOT
compression code (algorithm*100+level; the default 505 is ZSTD level 5).

The index needs the input in coordinate order. An unsorted input converts
fine but gets no index, and region queries on it read the whole file.

### Region Querying

Count the records overlapping a genomic region, as `samtools view -c` does.
```bash
# Usage: ./tools/ramntupleview [input.root] "[chromosome]:[start]-[end]"
./tools/ramntupleview output.root "chr1:10150-10300"
```

### Dumping back to SAM

`ramdump` writes a RAM file out as SAM and takes the `samtools view` options
`-h`, `-H`, `-c`, `-f`, `-F` and `-o`, so its output can be checked against
samtools directly.
```bash
# The whole file with its header; should reproduce the input SAM
./tools/ramdump -h output.root > roundtrip.sam

# Count primary alignments in a region
./tools/ramdump -c -F 0x900 output.root "chr1:10150-10300"
```

## Benchmark Results

Tested with the HG00154 sample from the 1000 Genomes Project (196M reads, 72.1 GB SAM). RAM is RNTuple-only. The same sample compresses to **11.4 GB**.

### Region Query Performance (LZMA)

| Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|----------|---------|-----------|-------------|
| Small (100bp) | 4.10 | 2.55 | 2.7 | 7 |
| Gene (BRCA2) | 1.20 | 1.15 | 37,308 | 42,961 |
| 10Mb | 7.40 | 6.67 | 446,331 | 2,977,922 |
| 100Mb | 6.46 | 6.36 | 448,688 | 2,852,438 |

### Region Query Performance (LZ4)

| Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|----------|---------|-----------|-------------|
| Small (100bp) | 2.98 | 1.67 | 4.2 | 7 |
| Gene (BRCA2) | 1.01 | 0.948 | 45,321 | 42,961 |
| 10Mb | 7.43 | 6.68 | 445,709 | 2,977,922 |
| 100Mb | 6.47 | 6.29 | 453,736 | 2,852,438 |

### Region Query Performance (ZLIB)

| Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|----------|---------|-----------|-------------|
| Small (100bp) | 2.85 | 1.73 | 4.0 | 7 |
| Gene (BRCA2) | 1.19 | 1.14 | 37,529 | 42,961 |
| 10Mb | 7.40 | 6.62 | 449,599 | 2,977,922 |
| 100Mb | 6.49 | 6.41 | 445,148 | 2,852,438 |

**Key findings:**
- Large-region throughput is about **445,000–454,000 reads/sec** across LZMA, LZ4, and ZLIB
- LZ4 is the fastest of the three on the 100Mb query (**453,736 reads/sec**)
- Small 100bp queries are dominated by setup overhead (a few seconds for 7 reads)

