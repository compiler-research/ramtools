# Changelog

All notable changes to RAMTools are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-07-13

First tagged release. RAMTools is experimental; formats and APIs may change without notice.

### Added

- SAM to RNTuple conversion producing the RAM (ROOT Alignment/Map) format on top of ROOT RNTuple ([#1](https://github.com/compiler-research/ramtools/pull/1), [#2](https://github.com/compiler-research/ramtools/pull/2)).
- BAM input with direct BAM to RNTuple conversion ([#20](https://github.com/compiler-research/ramtools/pull/20)).
- Region-based genomic querying through RAMNTupleView, comparable to samtools view ([#10](https://github.com/compiler-research/ramtools/pull/10), [#15](https://github.com/compiler-research/ramtools/pull/15)).
- Chromosome-based file splitting ([#9](https://github.com/compiler-research/ramtools/pull/9)).
- Command-line tools: `samtoram`, `samtoramntuple`, `bamtoramntuple`, `ramntupleview`, `rammerge`, `ramreader`, `ramrandom`, `checkindex`, `parsetreestats`.
- Benchmarking infrastructure covering conversion time and region-query throughput ([#4](https://github.com/compiler-research/ramtools/pull/4), [#6](https://github.com/compiler-research/ramtools/pull/6), [#7](https://github.com/compiler-research/ramtools/pull/7), [#8](https://github.com/compiler-research/ramtools/pull/8), [#54](https://github.com/compiler-research/ramtools/pull/54)).
- Continuous integration with build checks, clang-tidy, and code-coverage gating via Codecov ([#3](https://github.com/compiler-research/ramtools/pull/3), [#5](https://github.com/compiler-research/ramtools/pull/5), [#11](https://github.com/compiler-research/ramtools/pull/11), [#17](https://github.com/compiler-research/ramtools/pull/17)).

### Changed

- Set the project version to 0.1.0 to reflect experimental status.
- Refactored RAMNTupleView ([#10](https://github.com/compiler-research/ramtools/pull/10)).
- Upgraded ROOT in CI from v36 to v38 and aligned the clang-tidy ROOT version ([#47](https://github.com/compiler-research/ramtools/pull/47), [#57](https://github.com/compiler-research/ramtools/pull/57)).
- Upgraded the Codecov action to v5 ([#32](https://github.com/compiler-research/ramtools/pull/32)).

### Fixed

- Validate integer fields in SamParser and stop fRefVec pollution ([#59](https://github.com/compiler-research/ramtools/pull/59)).
- Correct DecodeSequence behavior ([#45](https://github.com/compiler-research/ramtools/pull/45)).
- Improve region-query logic in RAMNTupleView ([#15](https://github.com/compiler-research/ramtools/pull/15)).
- Restore casting-through-void and reinterpret-cast suppression ([#38](https://github.com/compiler-research/ramtools/pull/38)).

### Performance

- On the 1000 Genomes HG00154 sample (196M reads, 72GB SAM), RNTuple region queries run 1.4 to 2.5x faster than TTree for large regions. A 100Mb query reaches about 453k reads/sec with LZ4, versus about 198k reads/sec for TTree with ZLIB.
