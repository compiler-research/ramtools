# RAMtools Benchmark Suite

Performance benchmarks for the RAMtools SAM → ROOT (TTree / RNTuple) pipeline, built on
[Google Benchmark](https://github.com/google/benchmark). Every benchmark runs on
synthetic data out of the box and can be pointed at a real dataset from the command line.

> These benchmarks measure **performance** (time, file size, throughput). Correctness is
> covered separately by the gtest suite in [`../test/`](../test). Memory profiling is not
> yet wired up (planned follow-up).

## Benchmarks

| Binary | What it measures |
|--------|------------------|
| `sam_to_ram_benchmark` | SAM→TTree vs SAM→RNTuple conversion, output sizes, compression ratio |
| `conversion_time_benchmark` | SAM→TTree and SAM→RNTuple conversion time separately, output size |
| `region_query_benchmark` | Region-query latency, TTree vs RNTuple, across a scale ladder of regions |
| `chromosome_split_benchmark` | RNTuple parallel chromosome split vs a samtools split pipeline (needs `samtools`) |
| `columnar_read_benchmark` | Single-column (`flag`/`mapq`) scan vs full-record read — RNTuple's columnar advantage |

## Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # RAMTOOLS_BUILD_BENCHMARKS=ON by default
cmake --build build
```

Binaries land in `build/benchmark/`.

## Run

### Everything, with JSON output (one command)

```bash
cmake --build build --target benchmark
```

Runs every benchmark and writes one `*.json` result file per binary into
`build/benchmark/` (the build tree — nothing is written to the source tree).

### A single benchmark

```bash
./build/benchmark/region_query_benchmark
./build/benchmark/columnar_read_benchmark --benchmark_out=cols.json --benchmark_out_format=json
```

All standard Google Benchmark flags work (`--benchmark_filter=`, `--benchmark_repetitions=`,
`--benchmark_out=`, `--benchmark_out_format=json|csv`, …).

### On a real dataset

By default the benchmarks generate a synthetic SAM (`GenerateSAMFile`). Pass `--sam=` (and
optionally pre-converted `.root` files) to use your own data:

```bash
# Convert + query a real SAM (region_query / columnar convert once, then read):
./build/benchmark/region_query_benchmark   --sam=/data/HG002.sam --regions=all
./build/benchmark/conversion_time_benchmark --sam=/data/HG002.sam
./build/benchmark/columnar_read_benchmark   --sam=/data/HG002.sam

# Reuse already-converted files (skip conversion):
./build/benchmark/region_query_benchmark \
    --ttree-root=/data/HG002_ttree.root --rntuple-root=/data/HG002_rntuple.root
```

To forward dataset flags to *all* binaries at once and collect JSON:

```bash
scripts/run_benchmarks.sh build -- --sam=/data/HG002.sam --regions=all
```

## Command-line flags (shared)

Parsed by `benchmark_config.h`; stripped before Google Benchmark sees argv.

| Flag | Default | Meaning |
|------|---------|---------|
| `--sam=PATH` | *(generate)* | SAM file to benchmark |
| `--ttree-root=PATH` | *(derive from `--sam`)* | Pre-converted TTree `.root` |
| `--rntuple-root=PATH` | *(derive from `--sam`)* | Pre-converted RNTuple `.root` |
| `--compression=N` | `505` (ZSTD-5) | ROOT compression setting for RNTuple |
| `--quality=N` | `0` | Quality policy bitmask |
| `--threads=N` | `4` | Worker threads (split benchmarks) |
| `--reads=N` | `100000` | Synthetic read count when no `--sam` given |
| `--regions=all\|i,j,k` | `0,3,6,9` | Region indices for `region_query` |

When no `--sam`/`.root` is supplied, temp synthetic files are generated and cleaned up
automatically.

## Reading the results

Each row reports wall time, CPU time, iterations, and benchmark-specific counters:
`file_size_mb`, `ttree_size_mb` / `rntuple_size_mb` / `compression_ratio`,
`reads_per_second`, `size_MB`, `threads`, `region_idx`. The label column shows the read
count for the run (e.g. region matches, dataset rows).

The JSON files are consumed by the analysis notebooks in
[`io_performance/`](io_performance) (`ramview_perf.ipynb`, `samtoram_perf.ipynb`).

## Reproducibility

When reporting numbers, record the environment:

```
Date:
CPU model / cores:
RAM:
Disk (SSD/NVMe/HDD):
OS / kernel:
ROOT version:
samtools version (for chromosome_split):
Build type: Release
Dataset (name, #reads, source):
Command line:
```

## Notes

- `chromosome_split_benchmark` shells out to `samtools`; without it installed those rows
  produce no output.
