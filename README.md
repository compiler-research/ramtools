# RAMTools

RAM (ROOT Alignment/Map) is a file format for aligned sequencing reads, the same
data as a BAM, stored column by column in ROOT's RNTuple format. A RAM file is
about a quarter smaller than the BAM of the same reads, needs no reference to be
read, and answers region queries such as `chr13:32889611-32973805`. RAMTools
converts SAM and BAM files to RAM, queries them, and writes them back out as SAM
so every result can be checked against samtools.

The tools:

| Tool | What it does |
|------|--------------|
| `samtoramntuple` | converts a SAM file to RAM |
| `bamtoramntuple` | converts a BAM file to RAM |
| `ramdump` | writes a RAM file, or a region of it, out as SAM; takes the `samtools view` options |
| `ramntupleview` | counts the records in a region and reports the query time |

## Installation

RAMTools builds on Linux. The steps below are for Ubuntu; other distributions
need the same three things: ROOT, htslib and a C++17 compiler with CMake 3.16 or
newer.

### 1. System packages

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git pkg-config libhts-dev libtbb-dev libvdt-dev
```

`libhts-dev` is htslib, which reads the BAM input. `libtbb-dev` and `libvdt-dev`
are needed by ROOT. Quality scores are compressed with fqzcomp from htscodecs,
the codec CRAM 3.1 uses; if `libhtscodecs-dev` is installed it is used,
otherwise the build downloads htscodecs 1.6.7 and compiles the two files it
needs.

### 2. ROOT

RAMTools needs ROOT 6.38 or newer. The quickest way is the prebuilt archive from
the ROOT project; pick the one matching your Ubuntu version from
https://root.cern/install/all_releases/. For Ubuntu 24.04:

```bash
wget https://github.com/root-project/root/releases/download/v6-38-06/root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
sudo tar -xzf root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz -C /opt/
```

ROOT has to be put on your path in every terminal you use it from:

```bash
source /opt/root/bin/thisroot.sh
```

Add that line to your `~/.bashrc` to make it permanent. `root --version` should
then print 6.38. If you use conda, `conda install -c conda-forge root` works as
well and needs no `thisroot.sh`.

### 3. Build RAMTools

```bash
git clone https://github.com/compiler-research/ramtools.git
cd ramtools
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The first configure downloads GoogleTest and Google Benchmark for the tests and
benchmarks, and htscodecs unless `libhtscodecs-dev` is installed. Without
network access, install `libhtscodecs-dev` and turn the tests and benchmarks
off with `-DRAMTOOLS_BUILD_TESTS=OFF -DRAMTOOLS_BUILD_BENCHMARKS=OFF`.

The tools are in `build/tools/`. Running one without arguments prints its
usage:

```bash
./build/tools/samtoramntuple
```

The tests take under a minute:

```bash
ctest --test-dir build
```

## Converting to RAM

From a SAM file:

```bash
./build/tools/samtoramntuple reads.sam reads.ram
```

From a BAM file:

```bash
./build/tools/bamtoramntuple reads.bam reads.ram
```

Without an output name the file is written next to the input, `reads.sam`
becoming `reads.ram`. Options for either tool:

| Option | Effect |
|--------|--------|
| `-compression N` | ROOT compression code, algorithm times 100 plus level; the default 505 is ZSTD level 5 (see below) |
| `-threads N` | converts on N threads, see below |
| `-illumina` | stores quality scores in Illumina's 8 bins, which makes the file smaller |
| `-dropqual` | stores no quality scores |
| `-split` | `samtoramntuple` only: one file per reference, `reads_chr1.root`, `reads_chr2.root`, and so on |

With `-threads N`, `samtoramntuple` reads the input in 64 MB blocks and every
thread parses, encodes and compresses blocks of its own. The blocks are written
in input order, so the records come out in the order they went in, whatever
the thread count. It needs about 2 × N × 64 MB of memory on top of ROOT's own.
`bamtoramntuple`, and `samtoramntuple` with `-split`, parse on one thread and
use the others to compress.

Sort the input by coordinate first, as `samtools sort` does. Region queries
find a region by binary search over the positions, which only works on a sorted
file. An unsorted file converts fine, but every query on it reads the whole
file; the converter warns when that happens.

## Querying a region

Regions are 1-based and inclusive, as in samtools:

| Region | Records |
|--------|---------|
| `chr1` | every record on `chr1` |
| `chr1:1000-2000` | records overlapping bases 1000 to 2000 |
| `chr1:1000` | records overlapping base 1000 only; samtools reads this as 1000 to the end of `chr1` |
| none, or `*` | every record in the file |

Reference names are whatever the input used, so a GRCh37 file is queried as
`13:32889611-32973805`, not `chr13`.

Count the records overlapping a region, as `samtools view -c` does:

```bash
./build/tools/ramdump -c reads.ram chr1:10150-10300
```

Print them as SAM lines:

```bash
./build/tools/ramdump reads.ram chr1:10150-10300
```

`ramdump` takes the `samtools view` options `-h` (include the header), `-H`
(header only), `-c` (count only), `-f` and `-F` (keep records with all, or none,
of these FLAG bits) and `-o FILE`. For example, primary alignments only:

```bash
./build/tools/ramdump -c -F 0x900 reads.ram chr1:10150-10300
```

`ramntupleview` counts as well and prints how long the query took:

```bash
./build/tools/ramntupleview reads.ram chr1:10150-10300
```

## Checking a RAM file against samtools

A RAM file holds everything the input did, so it can be written back out and
compared:

```bash
./build/tools/ramdump -h reads.ram > roundtrip.sam
samtools view -c reads.bam chr1:10150-10300
./build/tools/ramdump -c reads.ram chr1:10150-10300
```

The two counts agree, and for a SAM input `roundtrip.sam` is the input byte for
byte.

## Benchmarks

HG00154 from the 1000 Genomes Project: 196,040,370 records, 72.06 GB as SAM,
coordinate sorted, GRCh37 reference names (`1`, not `chr1`). Run on an Intel
i3-7020U laptop (2 cores, 4 threads, 12 GB of memory) with input and output on
the same hard disk, samtools 1.13. Every output holds the same 196,040,370
records.

### File size

![File size of HG00154 in each format](assets/benchmark_sizes.svg)

BAM is `samtools view -b`, CRAM is `samtools view -C -T ref.fasta`, and each RAM
file is `samtoramntuple -compression N` with the setting shown. A CRAM cannot
be read without its reference, so the stacked bar is the size of a
self-contained copy. Uncompressed RAM (`-compression 0`) is 83.7 GB, larger
than the SAM, because every string and vector column carries an offset column
that the codecs then remove. ZSTD 9 is the only setting that comes in under
CRAM plus its reference, at ten times the conversion time of ZSTD 5 for 9% less
space.

### Conversion

From the SAM file. CPU is user plus system time; wall time is bounded by the
hard disk, which reads and writes about 77 MB/s.

| Output | Command | Wall | CPU (s) | Size (GB) |
|--------|---------|------|---------|-----------|
| BAM | `samtools view -@ 4 -b` | 24:08 | 3,589 | 15.22 |
| CRAM | `samtools view -@ 4 -C -T ref.fasta` | 12:34 | 1,817 | 7.77 |
| RAM, ZSTD 5 | `samtoramntuple` | 1:30:06 | 5,258 | 11.43 |
| RAM, ZSTD 5 | `samtoramntuple -threads 4` | 32:59 | 7,664 | 11.43 |
| RAM, ZSTD 3 | `samtoramntuple -compression 503 -threads 4` | 17:09 | 3,513 | 11.74 |
| RAM, ZSTD 1 | `samtoramntuple -compression 501 -threads 4` | 14:32 | 1,273 | 12.27 |

Almost all of the conversion time is ZSTD; parsing the SAM is about 4% of it.
Lowering the level is what makes it faster: ZSTD 3 needs half the CPU of ZSTD 5
for a 3% larger file, and ZSTD 1 converts faster than samtools writes a BAM and
still gives a file a fifth smaller.

### Region queries

Records overlapping a region, counted with `samtools view -c` for BAM and CRAM
and `ramdump -c` for RAM (ZSTD 5). Best of three runs, warm cache. All three
formats return the same records for every region.

| Region | Records | BAM (s) | CRAM (s) | RAM (s) |
|--------|---------|---------|----------|---------|
| 100 bp, `1:1000000-1000100` | 2 | 0.15 | 0.10 | 0.48 |
| BRCA2, `13:32889611-32973805` | 6,057 | 0.10 | 0.09 | 0.39 |
| 10 Mb, `1:10000000-20000000` | 582,985 | 0.80 | 1.75 | 0.40 |
| 100 Mb, `2:1-100000000` | 7,048,385 | 8.59 | 21.43 | 0.85 |

A RAM query spends about 0.4 s starting ROOT and opening the file, then finds
the region by binary search and counts about 15 million records per second,
since a count reads only the position and CIGAR columns. It is slower than
samtools on small regions and faster from about 10 Mb up. The codec makes
little difference: the 10 Mb query takes 0.46 to 0.55 s across all of them.

### Compression settings

`-compression N` takes a ROOT compression code, algorithm times 100 plus level:
1 is ZLIB, 2 LZMA, 4 LZ4, 5 ZSTD, with levels 1 to 9, and 0 means no
compression. The default is 505, ZSTD level 5. For a quick conversion use 501
or 503 with `-threads`; for the smallest file, 509 at a much higher conversion
cost.

## Contributing

Bug reports, test cases and pull requests are welcome. [CONTRIBUTING.md](CONTRIBUTING.md)
describes the development build, the test suites, the formatting and clang-tidy
checks that run on every pull request, and what a reviewer looks for. The short
version: branch from `develop`, keep each PR to one topic, and put the samtools
command that confirms your result in the description.

## Reporting a problem

Open an issue at https://github.com/compiler-research/ramtools/issues. Include
the command you ran, its output, your ROOT version, and for a wrong result the
samtools command that gives the answer you expected. A few SAM lines that
reproduce the problem are the most useful thing you can attach.

## License

RAMTools is released under the Apache License 2.0; see [LICENSE](LICENSE).
