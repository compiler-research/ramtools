# Changelog

All notable changes to RAMTools are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html); see "Versioning" in
the README for what that means before 1.0.

## [Unreleased]

### Added

- `bench/ramtools-bench.py` reproduces the benchmarks on any Linux or macOS
  machine: it downloads chr20 of NA12878 from EBI and checks it against known
  md5 sums, builds two revisions the same way, measures conversion, size, the
  round trip, dumps, regions and statistics in alternating order, checks every
  output against samtools, and writes a report. See `bench/README.md`
  ([#116]). On Linux it runs each command under a small helper it compiles,
  because a command started from Python counted Python's 20 to 30 MB in its
  peak memory: `samtools idxstats` showed 22 to 30 MB instead of 10 ([#117]).
- `-fqzcomp` on both converters compresses lossless quality scores with
  fqzcomp, the codec of CRAM 3.1, in blocks of 10,000 records. On 2 million
  records of chromosome 11 the file is 8.1% smaller, conversion takes 20 to 30%
  longer and a full dump about 30% longer; region counts are unchanged. It is
  off by default, and files written with it need this release or later to read
  their qualities ([#90], [#99]).
- `ramstats idxstats` prints the records per reference as `samtools idxstats`
  does. The converters count them while writing and store the counts in
  METADATA, so it reads no record: 0.2 s and 5 MB for any file size. The format
  version becomes 3, and files written before must be converted again ([#102]).
- `ramstats flagstat` prints what `samtools flagstat` prints, in its default,
  TSV (`-O tsv`) and JSON (`-O json`) formats. It reads only the FLAG,
  reference, mate reference and MAPQ columns. On chr20 of NA12878 at 30x
  (17.7 million records) it takes 1.3 s of CPU and 20 MB on one thread, against
  9.4 s for `samtools flagstat` on the BAM and 3.9 s on the CRAM 3.1. The file
  format is unchanged ([#108]).

### Changed

- `samtoramntuple -threads N` holds up to N + 1 blocks of 64 MB in memory
  instead of 2N + 1. On 2 million records of chromosome 11 the peak goes from
  695 to 426 MB on 4 threads and from 1052 to 725 MB on 8, in the same time
  ([#98]).
- RAM files record a format version, and a reader opens only files of its own
  version: files written by 0.1.0 must be converted again. Until 1.0 the format
  may change between minor releases; see "Versioning" in the README ([#101]).
- The SAM header is stored in METADATA instead of a ROOT `TList`, so
  `ramdump -H` and `-h` no longer start ROOT's interpreter. On 2 million
  records of chromosome 11, `-H` goes from 0.86 s and 50 MB to 0.35 s and 6 MB,
  a 1 kb region with `-h` from 0.86 s and 90 MB to 0.40 s and 47 MB, and the
  whole file with `-h` from 6.2 s and 96 MB to 5.3 s and 55 MB ([#101]).
- Library: `RAMFile::Open()` opens a RAM file and keeps its reference names
  and METADATA in the object it returns instead of in static members, so
  several files can be open in one process. `ramntuplescan()` and
  `QualityBlockReader` take that object ([#105]).
- Library: the converters keep the state of the file they write, the
  reference tables and the METADATA values, in their own objects instead of in
  static members of `RAMNTupleRecord`, so several conversions can run in one
  process at once. Before, two at once mixed their reference names ([#107]).
- The tags NM, AS, XS, MQ, RG, PG, MD, MC, SA and XA are stored in columns of
  their own instead of as text, RG and PG as ids into a dictionary in METADATA.
  The other tags stay text, and `ramdump -h` still reproduces the input byte for
  byte. On chr20 of NA12878 at 30x the file is 9.5% smaller (622.5 to
  563.1 MB), conversion takes 5 to 13% less CPU, and `ramdump -h` 8% less. The
  format version becomes 5, and files written before must be converted again
  ([#111], [#112]).
- PNEXT and TLEN are stored as the difference from what POS, the strands and
  the CIGAR predict, in the columns `pnext_delta` and `tlen_delta`, and still
  come back exactly. On chr20 of NA12878 at 30x the file goes from 563.1 to
  531.0 MB (-5.7%), and conversion and dumps take the same time. A program
  that reads the columns directly gets the SAM values with
  `RAMNTupleRecord::GetPNEXT()` and `GetTLEN()`; see CONTRIBUTING. The format
  version becomes 6 ([#113]).
- Read names are compressed with tok3, the name codec of CRAM 3.1, in blocks
  of 10,000 records, instead of stored in a `qname` column. On chr20 of NA12878
  at 30x the file goes from 531.0 to 478.2 MB (-9.9%), and `ramdump -h` still
  reproduces the input byte for byte. Small region queries, `ramdump -c`,
  `flagstat` and `idxstats` take the same time, and so does conversion. Reading
  every name costs CPU: on macOS a dump of the whole file takes about 6% longer
  (28 s) and 28% more CPU, with the names decoded ahead on a second thread. The
  format version becomes 7 ([#114]).
- Library: `QualityBlocks.h` becomes `RecordBlocks.h` and `QualityBlockWriter`
  becomes `BlockWriter`, which also writes the name blocks. `NameBlockReader`
  gives the name of any row; a record read from a file has an empty
  `GetQNAME()`. METADATA's `qual_block_ends` becomes `block_ends` ([#114]).
- Lossless quality scores are stored in the blocks of 10,000 records with
  rANS order 1, the codec CRAM 3.1 uses for them by default, instead of in the
  `qual` column. On chr20 of NA12878 at 30x the file goes from 478.2 to
  435.7 MB (-8.9%), and `ramdump -h` still reproduces the input byte for byte.
  Conversion takes the same CPU or up to 10% less. Reading every quality costs
  CPU: on macOS a dump of the whole file takes about 20% more CPU and 8% more
  time; small region queries, `ramdump -c`, `flagstat` and `idxstats` take the
  same time. The format version becomes 8 ([#115]).
- `-fqzcomp` compresses the same blocks with fqzcomp instead of rANS. On chr20
  of NA12878 at 30x, whose qualities are binned, it saves about 1% of the
  qualities and takes three times the CPU to read them ([#115]).
- Library: lossless qualities go into blocks whatever the quality policy, so
  `kQualInBlock` only says that a record's quality is in its block, and
  `kFqzcomp` in the policy asks for fqzcomp. A record read from a file has an
  empty `GetQUAL()` when its quality is in a block; `QualityBlockReader` gives
  the quality of any record ([#115]).

### Removed

- Reading files written by 0.1.0 and earlier, and with them SEQ in the 4-bit
  packing: `kSeqRaw`, `RAMNTupleUtils::EncodeSequence` and
  `RAMNTupleUtils::DecodeSequence` ([#101]).
- `RAMNTupleRecord::OpenRAMFile()` and `ReadAllRefs()`; open files with
  `RAMFile::Open()` ([#105]).
- `RAMNTupleRecord`'s static writer state and the calls on it:
  `InitializeRefs()`, `WriteAllRefs()`, `GetRnameRefs()`, `GetRnextRefs()` and
  the static getters and setters of METADATA values. A converter now fills a
  `RAMFileMetadata` and calls `WriteMetadata()`. Also removed are the record's
  `SetRNAME()`, `SetRNEXT()`, `SetREFID(name)`, `SetREFNEXT(name)`,
  `GetRNAME()`, `GetRNEXT()` and `Print()`. A record holds ids, and
  `RAMFile::GetRname()` and `GetRnext()` give their names ([#107]).
- `RAMNTupleRecord`'s `AddTag()`, `SetOPT()`, `ClearTags()`, `ResetNOPT()`,
  `GetTags()`, `GetNOPT()` and `GetOPT()`. `tags` now holds only the tags
  without a column of their own. Converters split the tags with `RAMTagWriter`,
  and `RAMFile::AppendTags()` writes all of them back ([#111]).
- `RAMNTupleRecord`'s `SetPNEXT()` and `SetTLEN()`, and `GetPNEXT()` and
  `GetTLEN()` without arguments. Converters set both with `SetMate()`, and
  readers pass `RAMFile::MateOnSameRef()` to `GetPNEXT()` and `GetTLEN()`
  ([#113]).

## [0.1.0] - 2026-10-14

The first tagged release. RAMTools began in 2017 at CERN openlab as a format
stored in ROOT's TTree; this release is its rewrite on RNTuple. Files written by
it stay readable by every later release.

### Added

- The RAM format on ROOT's RNTuple and `samtoramntuple`, which converts SAM to
  it ([#1], [#2]).
- `bamtoramntuple`, which converts BAM through htslib ([#20]).
- `ramdump`, which writes a file or a region of it back out as SAM and takes the
  `samtools view` options `-h`, `-H`, `-c`, `-f`, `-F` and `-o` ([#75]).
- `ramntupleview`, which counts the records in a region and reports the query
  time ([#10], [#69]).
- `-split`, one output file per reference, written as a stream and indexed like
  any other file ([#9], [#79]).
- `-compression N` on both converters, taking a ROOT compression code ([#73]).
- `-threads N` on both converters; `samtoramntuple` parses, encodes and
  compresses 64 MB blocks on every thread and writes them in input order
  ([#86], [#87], [#88]).
- Region queries find the region by binary search over the positions of a
  coordinate-sorted file; the converters record whether the input was sorted,
  and an unsorted file is read in full ([#79], [#83]).
- Benchmarks: a SAM generator and Google Benchmark suites for conversion and
  region queries ([#4], [#6], [#7], [#8]), and scripts that run them on real data
  ([#54]).
- Continuous integration: build and tests, clang-tidy review, formatting checks
  and Codecov gating ([#3], [#5], [#11], [#17], [#18], [#19], [#21], [#32], [#35],
  [#47], [#57], [#77]).
- Tests for the index, the quality encodings and the record getters ([#39],
  [#41], [#43]).
- The Apache License 2.0 ([#60]); `CONTRIBUTING.md`, issue and pull request
  templates, and a README that walks through installation, conversion and
  queries ([#12], [#89]).

### Changed

- ROOT 6.38 is the minimum version ([#76]).
- SEQ is stored as uppercase text instead of 4-bit packed bases, which
  compresses better: 3.9% smaller files on chromosome 11. Files with packed
  sequences still read ([#85]).
- Opening a RAM file no longer starts ROOT's interpreter. On 2 million records
  of chromosome 11, counting a 10 Mb region went from 0.59 s and 150 MB to
  0.26 s and 19 MB ([#92]), and dumping a 1 kb region from 0.61 s and 204 MB to
  0.30 s and 78 MB ([#93]). `-h` and `-H` still start it, to read the header.
- Region queries decide overlap by each read's CIGAR span and accept `rname`,
  `rname:pos`, `rname:start-end` and `*` ([#15]).

### Removed

- The TTree backend ([#71]).
- `RAMNTupleConverter`, replaced by `ramdump` ([#75]).
- The sparse index of mapped reads ([#16]), replaced by the binary search
  ([#83]).

### Fixed

- Sequences decoded wrongly by `DecodeSequence` ([#45]).
- The SAM parser accepted non-numeric and out-of-range integer fields, and a
  query for an unknown reference added it to the reference list ([#59]).
- Float tags lost precision ([#64]).
- A CIGAR of `*` came back as an empty field ([#66]).
- Illumina quality binning was indexed by the ASCII character instead of the
  Phred score ([#68]).
- Fields lost or shifted on the way in or out: lines over 10 kB, consecutive
  tabs, unknown and lowercase bases, unknown CIGAR operators ([#70]).
- The SAM header was written as one key per line and could not be read back
  ([#72]).
- A mapped read with a malformed CIGAR such as `10Q` was stored without one; it
  is now reported with its line number and skipped ([#74]).
- A region query past the last index anchor of a reference started from the
  first row and read the whole file ([#82]); the index was later replaced by the
  binary search ([#83]).
- README examples that did not run ([#78]).
- `bamtoramntuple` did not link when htslib is not in a default library
  directory, as on macOS ([#95]).
- clang-tidy findings and suppressions ([#14], [#38]).

[Unreleased]: https://github.com/compiler-research/ramtools/compare/v0.1.0...develop
[0.1.0]: https://github.com/compiler-research/ramtools/releases/tag/v0.1.0
[#1]: https://github.com/compiler-research/ramtools/pull/1
[#2]: https://github.com/compiler-research/ramtools/pull/2
[#3]: https://github.com/compiler-research/ramtools/pull/3
[#4]: https://github.com/compiler-research/ramtools/pull/4
[#5]: https://github.com/compiler-research/ramtools/pull/5
[#6]: https://github.com/compiler-research/ramtools/pull/6
[#7]: https://github.com/compiler-research/ramtools/pull/7
[#8]: https://github.com/compiler-research/ramtools/pull/8
[#9]: https://github.com/compiler-research/ramtools/pull/9
[#10]: https://github.com/compiler-research/ramtools/pull/10
[#11]: https://github.com/compiler-research/ramtools/pull/11
[#12]: https://github.com/compiler-research/ramtools/pull/12
[#14]: https://github.com/compiler-research/ramtools/pull/14
[#15]: https://github.com/compiler-research/ramtools/pull/15
[#16]: https://github.com/compiler-research/ramtools/pull/16
[#17]: https://github.com/compiler-research/ramtools/pull/17
[#18]: https://github.com/compiler-research/ramtools/pull/18
[#19]: https://github.com/compiler-research/ramtools/pull/19
[#20]: https://github.com/compiler-research/ramtools/pull/20
[#21]: https://github.com/compiler-research/ramtools/pull/21
[#32]: https://github.com/compiler-research/ramtools/pull/32
[#35]: https://github.com/compiler-research/ramtools/pull/35
[#38]: https://github.com/compiler-research/ramtools/pull/38
[#39]: https://github.com/compiler-research/ramtools/pull/39
[#41]: https://github.com/compiler-research/ramtools/pull/41
[#43]: https://github.com/compiler-research/ramtools/pull/43
[#45]: https://github.com/compiler-research/ramtools/pull/45
[#47]: https://github.com/compiler-research/ramtools/pull/47
[#54]: https://github.com/compiler-research/ramtools/pull/54
[#57]: https://github.com/compiler-research/ramtools/pull/57
[#59]: https://github.com/compiler-research/ramtools/pull/59
[#60]: https://github.com/compiler-research/ramtools/pull/60
[#64]: https://github.com/compiler-research/ramtools/pull/64
[#66]: https://github.com/compiler-research/ramtools/pull/66
[#68]: https://github.com/compiler-research/ramtools/pull/68
[#69]: https://github.com/compiler-research/ramtools/pull/69
[#70]: https://github.com/compiler-research/ramtools/pull/70
[#71]: https://github.com/compiler-research/ramtools/pull/71
[#72]: https://github.com/compiler-research/ramtools/pull/72
[#73]: https://github.com/compiler-research/ramtools/pull/73
[#74]: https://github.com/compiler-research/ramtools/pull/74
[#75]: https://github.com/compiler-research/ramtools/pull/75
[#76]: https://github.com/compiler-research/ramtools/pull/76
[#77]: https://github.com/compiler-research/ramtools/pull/77
[#78]: https://github.com/compiler-research/ramtools/pull/78
[#79]: https://github.com/compiler-research/ramtools/pull/79
[#82]: https://github.com/compiler-research/ramtools/pull/82
[#83]: https://github.com/compiler-research/ramtools/pull/83
[#85]: https://github.com/compiler-research/ramtools/pull/85
[#86]: https://github.com/compiler-research/ramtools/pull/86
[#87]: https://github.com/compiler-research/ramtools/pull/87
[#88]: https://github.com/compiler-research/ramtools/pull/88
[#89]: https://github.com/compiler-research/ramtools/pull/89
[#90]: https://github.com/compiler-research/ramtools/pull/90
[#92]: https://github.com/compiler-research/ramtools/pull/92
[#93]: https://github.com/compiler-research/ramtools/pull/93
[#95]: https://github.com/compiler-research/ramtools/pull/95
[#98]: https://github.com/compiler-research/ramtools/pull/98
[#99]: https://github.com/compiler-research/ramtools/pull/99
[#101]: https://github.com/compiler-research/ramtools/pull/101
[#102]: https://github.com/compiler-research/ramtools/pull/102
[#105]: https://github.com/compiler-research/ramtools/pull/105
[#107]: https://github.com/compiler-research/ramtools/pull/107
[#108]: https://github.com/compiler-research/ramtools/pull/108
[#111]: https://github.com/compiler-research/ramtools/pull/111
[#112]: https://github.com/compiler-research/ramtools/pull/112
[#113]: https://github.com/compiler-research/ramtools/pull/113
[#114]: https://github.com/compiler-research/ramtools/pull/114
[#115]: https://github.com/compiler-research/ramtools/pull/115
[#116]: https://github.com/compiler-research/ramtools/pull/116
[#117]: https://github.com/compiler-research/ramtools/pull/117
