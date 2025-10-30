# RAMTools - ROOT Alignment/Map Format Tools

RAMTools provides efficient tools for converting SAM files to ROOT's modern, columnar RNTuple format (RAM - ROOT Alignment/Map) and working with genomic alignment data.

## Features

- High-performance SAM to RAM conversion with an RNTuple backend
- Chromosome-based splitting for parallel processing
- Region-based querying capabilities

## Requirements

- ROOT 6.26+
- C++17 compatible compiler
- CMake 3.16+

## Quick Start
```bash
# 1. Build the tools
mkdir build && cd build
cmake ..
make -j$(nproc)

# 2. Convert a SAM file to the RAM format
./tools/samtoramntuple ../test/samexample.sam output.root

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
# (Creates output-chr1.root, output-chr2.root, etc.)
./tools/samtoramntuple input.sam output -split
```

### Region Querying (RNTuple)

Query a specific genomic region from a RAM file, similar to samtools view.
```bash
# Usage: ./tools/ramntupleview [input.root] "[chromosome]:[start]-[end]"
./tools/ramntupleview output.root "chr1:10150-10300"
```

## Benchmark Results

Tested with HG00154 sample from the 1000 Genomes Project (196M reads, 72GB SAM file):

![File Size Comparison](./img/file-size-comparison.png)

### Region Query Performance(LZMA Compression)

| Format | Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|--------|----------|---------|-----------|-------------|
| **TTree** | Small (100bp) | 4.18 | 1.50 | 4.7 | 7 |
| | Gene (BRCA2) | 1.87 | 1.79 | 24,010 | 42,961 |
| | 10Mb | 41.2 | 39.7 | 75,105 | 2,977,922 |
| | 100Mb | 36.3 | 35.4 | 80,492 | 2,852,438 |
| **RNTuple** | Small (100bp) | 4.10 | 2.55 | 2.7 | 7 |
| | Gene (BRCA2) | 1.20 | 1.15 | 37,308 | 42,961 |
| | 10Mb | 7.40 | 6.67 | 446,331 | 2,977,922 |
| | 100Mb | 6.46 | 6.36 | 448,688 | 2,852,438 |

### Region Query Performance (LZ4 Compression)

| Format | Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|--------|----------|---------|-----------|-------------|
| **TTree** | Small (100bp) | 3.23 | 1.86 | 3.8 | 7 |
| | Gene (BRCA2) | 0.941 | 0.842 | 51,030 | 42,961 |
| | 10Mb | 14.5 | 11.0 | 269,797 | 2,977,922 |
| | 100Mb | 9.05 | 9.05 | 315,088 | 2,852,438 |
| **RNTuple** | Small (100bp) | 2.98 | 1.67 | 4.2 | 7 |
| | Gene (BRCA2) | 1.01 | 0.948 | 45,321 | 42,961 |
| | 10Mb | 7.43 | 6.68 | 445,709 | 2,977,922 |
| | 100Mb | 6.47 | 6.29 | 453,736 | 2,852,438 |

### Region Query Performance (ZLIB Compression)

| Format | Region | Time (s) | CPU (s) | Reads/sec | Total Reads |
|--------|--------|----------|---------|-----------|-------------|
| **TTree** | Small (100bp) | 4.19 | 2.31 | 3.0 | 7 |
| | Gene (BRCA2) | 1.22 | 1.11 | 38,661 | 42,961 |
| | 10Mb | 18.2 | 16.3 | 183,021 | 2,977,922 |
| | 100Mb | 14.4 | 14.4 | 197,815 | 2,852,438 |
| **RNTuple** | Small (100bp) | 2.85 | 1.73 | 4.0 | 7 |
| | Gene (BRCA2) | 1.19 | 1.14 | 37,529 | 42,961 |
| | 10Mb | 7.40 | 6.62 | 449,599 | 2,977,922 |
| | 100Mb | 6.49 | 6.41 | 445,148 | 2,852,438 |

**Key Findings**: 
- RNTuple demonstrates **1.4-2.5x faster** query performance for large regions compared to TTree
- LZ4 compression provides the best query performance among all compression algorithms
- For a 100Mb region query: RNTuple+LZ4 processes **453,736 reads/sec** vs TTree+ZLIB's **197,815 reads/sec**

## TTree Implementation (Legacy)

ROOT scripts to convert a SAM file to a RAM (ROOT Alignment/Map) file using the older TTree format and to work with those files.

### Convert SAM to RAM with TTree
```bash
$ root
root [0] .x samtoram.C
root [1] .q
```

### Read a RAM file (TTree)
```bash
$ root
root [0] .x ramreader.C
root [1] .q
```

### View a specific region (TTree)

To view a region, the equivalent of `samtools view bamexample.bam chr1:10150-10300`:
```bash
$ root
root [0] .x ramview.C("ramexample.root","chr1:10150-10300")
root [1] .q
```

