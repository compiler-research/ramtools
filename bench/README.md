# Benchmarks on real data

`bench/ramtools-bench.py` reproduces the sizes and timings quoted in
`CHANGELOG.md` on any Linux or macOS machine. It downloads the data, builds two
revisions of ramtools the same way, measures them in alternating order, checks
every output against samtools, and writes a report.

## What you need

- Linux or macOS, Python 3.8 or newer, git, CMake 3.16 or newer and a C++17
  compiler.
- ROOT 6.38 or newer, set up in the shell (`source <ROOT>/bin/thisroot.sh`).
- htslib with its development files, as for building ramtools.
- samtools 1.16 or newer, built with libcurl, which reads the data from EBI.
- About 15 GB of disk for the data, and 3 GB more while all of chr20 is
  measured.

`bench/ramtools-bench.py check` lists what is there and what is missing.

## Running it

```bash
git clone https://github.com/compiler-research/ramtools.git
cd ramtools
source <ROOT>/bin/thisroot.sh
bench/ramtools-bench.py check
bench/ramtools-bench.py fetch                                  # once: 20 to 60 minutes
bench/ramtools-bench.py run --base v0.1.0 --head develop       # 2 million records
bench/ramtools-bench.py run --base develop~1 --head develop --input chr20
```

`--base` and `--head` take anything git does: a branch, a tag, a commit. The
default head, `WORKTREE`, is the checkout with its uncommitted changes. To
measure a pull request, fetch it first:
`git fetch origin pull/115/head:pr-115`, then `--head pr-115`.

Everything goes into `--work`, by default `../ramtools-bench-work` next to the
clone, or `$RAMTOOLS_BENCH_WORK`: the data, a build per revision (built once
and kept), and the results.

On a machine with several processor sockets, `--pin 0` runs every measured
command on NUMA node 0 with `numactl`, so that threads and memory stay on one
socket. On Linux, the `performance` CPU frequency governor gives steadier
timings; `check` says which governor is set.

## The data

| | |
|---|---|
| Sample | NA12878, 1000 Genomes 30x high-coverage set (NYGC), GRCh38, NovaSeq, qualities binned to four values |
| Source | `https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram`, chr20 only |
| Reference | chr20 of `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa`; md5 of the sequence `b18e6c531b0bd70e949a7fc20859cb01`, the M5 of chr20 in the CRAM |
| Records | 17,705,654; md5 of `samtools view` of the BAM, without the header, `89386f84300ce2a3b5100ec5760a197f` |

`fetch` checks the reference and the records against these and stops when
they differ; `--accept-different-data` goes on, and the report then shows
another md5 for the SAM. It writes, in `--work`:

- `GRCh38.chr20.fa`: `samtools faidx <reference> chr20`.
- `NA12878.chr20.bam`: `samtools view --no-PG -b -T GRCh38.chr20.fa -X <CRAM> <CRAI> chr20`, indexed.
- `NA12878.chr20.sam`: `samtools view --no-PG -h` of the BAM, about 9 GB.
- `NA12878.chr20.cram`: CRAM 3.1 with samtools' default profile, indexed, for
  the comparison with samtools.
- The first 2,000,000 and 50,000 records of the SAM, each also as BAM and CRAM.

`--no-PG` keeps the commands out of the header, so every machine gets the same
files: the report gives the md5 of the SAM it measured.

## What is measured

`run --input` picks the records: `small` (50,000, a minute or two, to try the
script out), `subset` (2 million, the default, for comparing changes) or
`chr20` (all 17.7 million, once per change). The steps:

| Step | Command | small | subset | chr20 |
|---|---|---:|---:|---:|
| convert SAM, 1 thread | `samtoramntuple in.sam out.ram` | 3 | 5 | 0 |
| convert SAM, threads | `samtoramntuple in.sam out.ram -threads 4` | 3 | 5 | 2 |
| convert BAM | `bamtoramntuple in.bam out.ram` | 3 | 5 | 1 |
| whole file | `ramdump -h` | 3 | 5 | 3 |
| large region | `ramdump file region`, about half the records (10 Mb on chr20) | 3 | 5 | 3 |
| three small regions | `ramdump file region`, 500 to 2,000 bp | 5 | 5 | 3 |
| count | `ramdump -c` of the large region | 5 | 5 | 3 |
| statistics | `ramstats flagstat`, `ramstats idxstats` | 5 | 5 | 3 |

The numbers are the timed runs per revision; `--repeat` and `--repeat-convert`
change them. The two revisions take turns, base and head, then head and base,
so a machine that slows down part of the way slows both alike. Reads have one
untimed run first, so the file is in the page cache for every timed run.

For each run the script records the wall time, the CPU time (user plus
system, of all the process's threads), the peak memory, and the load average
before it. Peak memory is the peak physical footprint on macOS and the peak
RSS on Linux; the two are not comparable with each other. On Linux a command
runs under a small helper that the script compiles with `cc`, because a child
forked from Python would count Python's own 20 to 30 MB in its peak.

samtools runs the same reads on the BAM and on the CRAM 3.1, one thread each,
as the yardstick (`--no-samtools` leaves it out).

After the timings, every result is checked: `ramdump -h` of the file converted
from SAM gives the SAM byte for byte; `ramdump` of the file converted from BAM
gives what `samtools view` prints; the count, the small regions, `flagstat`
and `idxstats` equal samtools' on the BAM. `run` exits with 1 when a check
fails.

## Reading the results

`run` prints the report and saves it as Markdown and JSON in `WORK/results`,
named after the time, the machine, the input and the two builds. The JSON has
every run; `bench/ramtools-bench.py report <file.json>` prints the report
again.

Each time is the median of the runs, with the fastest and the slowest in
brackets. Compare CPU times when the machine is busy: wall times then vary
with the load, CPU times much less. When the brackets of base and head
overlap, the difference is not shown by this run; run more repeats, or run
again when the machine is quiet. A size, a correctness check and an output are
exact.

To compare two machines, run `fetch` on both and check that the reports give
the same md5 for the SAM.
