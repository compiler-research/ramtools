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
- `ramdump -@ N` reads and formats the whole file on N threads, each taking
  every N-th block of 10,000 records with a reader of its own; the records come
  out in order. On chromosome 20 of NA12878 a full read takes 4.6 s on 8
  threads and 7.9 s on 4 instead of 16.5 s on one, with 546 MB instead of
  107 MB on 8; `samtools view -@ 8` takes 3.6 s from BAM and 7.7 s from
  CRAM 3.1. With a region, `-@` is ignored ([#124]).
- `ramdump` takes several regions and reads them in turn from the file kept
  open, each region's records in order, as `samtools view` does without `-M`;
  `-c` counts them all. On chromosome 20 of NA12878, 100 regions of 1 kb take
  0.73 s in one process instead of 6.55 s in a process each ([#125]).

### Changed

- `ramdump` with a region searches for it with only the reference and position
  columns open, opens the others at the first record found, and reads without
  RNTuple's cluster cache, which loaded every column of each cluster the search
  visited and of the next one. On chromosome 20 of NA12878 a region of 1 kb
  reads 1.0 MB of the file instead of 53.4 MB and takes 68 ms and 38 MB
  instead of 103 ms and 50 MB, or 80 ms instead of 246 ms from a cold cache.
  A count with `-f` or `-F` reads only the flags: 0.31 s instead of 1.76 s in
  a 10 Mb region. Library: `RAMFile::Open()` takes RNTuple's read options;
  without the cluster cache, a region of 1 kb in a file kept open takes 6.3 ms
  instead of 34.4 ms ([#119]).
- `ramstats flagstat` reads its four columns 65,536 rows at a time with
  RNTuple's bulk reads instead of one view call per value, which took four
  fifths of its time. On chromosome 20 of NA12878 it takes 0.44 s instead of
  1.38 s, against 5.9 s for `samtools flagstat` on the BAM and 2.9 s on the
  CRAM 3.1 ([#120]).
- `bamtoramntuple -threads N` fills and compresses the records on N threads, as
  `samtoramntuple` does: the BAM is read in blocks of 64 MB, each filled on a
  thread of its own and committed in input order, and htslib inflates the
  input on N threads instead of on every core of the machine. On chromosome 20
  of NA12878 a conversion on 8 threads takes 22.8 s instead of 81.2 s and
  203 s of CPU instead of 357 s, with 1136 MB instead of 419 MB. The file has a
  cluster per block, 115 of about 3.8 MB instead of 6 of about 70 MB, as from
  SAM. A read error in the BAM now stops the conversion with a message ([#122]).
- `ramdump` reads the fifteen numeric fields of a record in bulk, a run of rows
  at a time, instead of one view call per value. On chromosome 20 of NA12878 a
  full read takes 16.8 s instead of 18.1 s and a 10 Mb region 3.2 s instead of
  3.4 s; a small region reads as many pages as before ([#123]).
- `ramstats flagstat` prints counts the converters store in METADATA while
  they write the records, as `idxstats` does, instead of reading four columns of
  every record: on chromosome 20 of NA12878 it takes 0.05 s instead of 0.41 s,
  and as long on any file. `--scan` counts them again from the columns. The
  format version becomes 9, and files written before must be converted again
  ([#126]).
- The converters cap pages at 128 kB before compression instead of 64 kB.
  Opening a file reads the page list of every cluster, and a file now has half
  the pages. On all of NA12878 converted from BAM, 1.99 million pages instead
  of 3.99 million: `ramstats idxstats` takes 0.37 s and 157 MB instead of
  0.59 s and 268 MB, a region of 1 kb 0.59 s and 209 MB instead of 0.90 s and
  345 MB, and the file is 18.16 GB instead of 18.42 GB. On chromosome 20 a
  region of 1 kb in a file kept open takes 8.3 ms instead of 6.4 ms, because
  its pages are larger, and as long as before in a process of its own. Files
  written before read as before ([#127]).
- The bases go into blocks of 10,000 records, in a field `seqblock`, compressed
  with libdeflate at level 6, instead of into a column that RNTuple compresses.
  ROOT runs zstd at twice the level it is given, so the converters' 505 is
  zstd's level 10, which on text of four letters is both larger and twice as
  slow; the converters' help and the README said level 5. On chromosome 20 of
  NA12878 the file is 417.9 MB instead of 429.6 MB, and converting from SAM on
  one thread takes 112.5 s instead of 128.8 s. A full read takes 15.1 s instead
  of 16.3 s on one thread, and 25.7 s of CPU instead of 26.4 s: the block
  readers now decode into the memory of a block they are done with, and rANS
  decodes the qualities in place. Reading only the bases takes 2.9 s instead of
  3.9 s, but 4.3 s of CPU instead of 3.9 s.
  `record.seq` now holds only `*`: a program that reads the columns directly
  has to decode `seqblock`, with `SeqBlockReader` in the library. The format
  version becomes 10, and files written before must be converted again.
  libdeflate is the system's, or version 1.26 built from source ([#128]).
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
- `ramdump` writes numbers, CIGARs and qualities into its output without
  temporary strings. On chr20 of NA12878 at 30x a dump of the whole file takes
  27.9 s of CPU and 17.6 s instead of 29.8 s and 19.3 s (-6% and -9%, one
  socket of a Xeon E5-2690 v4), and a 10 Mb region 6% less CPU. The output is
  unchanged ([#<N>]).
- Library: `QualityBlockReader::Get` returns a `std::string_view` into the
  decoded block, valid until the next call, instead of a copy.
  `RAMNTupleRecord::AppendCIGAR` and `RAMNTupleUtils::AppendInt` append to a
  string without making a temporary one ([#<N>]).

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
[#119]: https://github.com/compiler-research/ramtools/pull/119
[#120]: https://github.com/compiler-research/ramtools/pull/120
[#122]: https://github.com/compiler-research/ramtools/pull/122
[#123]: https://github.com/compiler-research/ramtools/pull/123
[#124]: https://github.com/compiler-research/ramtools/pull/124
[#125]: https://github.com/compiler-research/ramtools/pull/125
[#126]: https://github.com/compiler-research/ramtools/pull/126
[#127]: https://github.com/compiler-research/ramtools/pull/127
[#128]: https://github.com/compiler-research/ramtools/pull/128
[#<N>]: https://github.com/compiler-research/ramtools/pull/<N>
