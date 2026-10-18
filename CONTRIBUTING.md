# Contributing to RAMTools

Thanks for helping. This page is about getting a change in: how to build for
development, what the tests expect, what CI checks, and what a reviewer looks
for. If you only want to use the tools, the README is enough.

## Reporting a problem

Open an issue at https://github.com/compiler-research/ramtools/issues with:

- the command you ran and what it printed, including any warning on standard
  error;
- the ROOT version (`root --version`) and the operating system;
- what the input looks like: SAM or BAM, coordinate sorted or not, roughly how
  many records, and the reference naming (`chr1` or `1`);
- for a wrong result, the samtools command that gives the answer you expected
  and its output. The tools are meant to agree with samtools record for record,
  so a mismatch is a bug and the samtools output is the test case.

A small input that reproduces the problem is worth more than any description.
Ten SAM lines are usually enough.

## Building for development

Follow the installation steps in the README, then configure a debug build with
the tests on:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DRAMTOOLS_BUILD_TESTS=ON
cmake --build build -j$(nproc)
ctest --test-dir build
```

The first configure downloads GoogleTest and Google Benchmark, so it needs
network access once.

## Tests

The suites live in `test/`:

| Suite | Covers |
|-------|--------|
| `ramcoretests` | the SAM parser, the encoders, coordinate order and region queries, with small hand-written SAM inputs whose expected counts are stated in the test |
| `bam_conversion` | BAM input through htslib, compared with the same reads converted from SAM |
| `chromosome_split_test` | split output against the unsplit file |
| `parallel_conversion_test` | `-threads` conversion with tiny blocks: the same records as one block, the sort check across block boundaries, long headers, malformed lines and read errors |
| `record_blocks_test` | name and quality blocks, with rANS and with `-fqzcomp`: which names and qualities go into a block, reading them back in any order, the codec each block records, that names are not a column, and damaged blocks and block ends |
| `ramdump_test` | the built tools on a fixture: the dump reproduces the input and counts follow samtools' rules |

Run one suite on its own while you work, it is faster than `ctest`:

```bash
./build/test/ramcoretests --gtest_filter='*Cigar*'
```

When you fix a bug, add the input that exposed it as a test. When you change
what a query returns, state the expected count in the test and say where it
comes from.

## Style and CI

Every pull request runs four checks:

- **precheckin** runs `git clang-format` with clang-format 18 over the lines
  the PR touches and fails on any difference. Format only what you changed:

  ```bash
  git clang-format --binary clang-format-18 develop
  ```

  Formatting whole files makes the diff unreadable and is not what the check
  asks for.
- **Build and Test** builds on Ubuntu with the ROOT release named in
  `.github/workflows/ci.yml` and runs `ctest`.
- **Code Coverage** posts to Codecov, which fails a PR when fewer than 85% of
  the lines it adds or changes are covered, or when total coverage drops by
  more than 2%.
- **review** runs clang-tidy with the repository's `.clang-tidy` and posts
  findings as PR comments. New code should be clean. Findings on lines you did
  not change are not yours to fix in the same PR.

## Pull requests

- Branch from `develop` and open the PR against it.
- One topic per PR, and one change per commit. Write the commit message as
  plain sentences that say what was wrong and what changed; the first line is
  the summary.
- Put the evidence in the description: what you ran, what it printed, and what
  you compared it against. Reviewers check claims against samtools, so give
  them the command that shows the agreement.
- Add a line under `[Unreleased]` in `CHANGELOG.md`, in Added, Changed, Removed
  or Fixed, with a link to the PR. Say what a user notices and, for speed or
  size, by how much. Changes users do not see, such as tests, CI or a
  refactoring, need no line.
- Keep the code minimal. A reviewer should be able to read the whole diff in
  one sitting; if a change is getting long, split it.
- Public functions and members in `inc/` carry a short `///` comment saying
  what they do and what the caller must ensure. That is what the documentation
  is generated from, so keep it accurate when you change behaviour.

## Per-file state

A file being read keeps what its METADATA holds, the reference names, the
longest span, the sort flag, the block ends, the SAM header and the records
per reference, in the `RAMFile` that `RAMFile::Open()` returns. Pass that
object to `ramntuplescan`, `NameBlockReader`, `QualityBlockReader` and the
reference name lookups;
nothing about an open file is static, so several can be open in one process.
The test `TwoOpenFilesKeepTheirOwnState` checks that.

A converter keeps the state of the file it writes in its own objects: two
`RAMNTupleRefs` tables, for the ids it stores in `refid` and `refnext`, and a
`RAMFileMetadata` it fills while it writes the records and passes to
`RAMNTupleRecord::WriteMetadata()` at the end. A record holds ids only, and
`RAMNTupleRecord` has no static members, so conversions can run side by side.
The test `TwoConversionsAtOnceKeepTheirOwnState` checks that.

A new METADATA field goes into `RAMFileMetadata`, is written in
`RAMNTupleRecord::WriteMetadata()` and read in `RAMFile::Open()`.

## Changing the file format

Any change to what a RAM file holds, in the records or in METADATA, bumps
`RAMNTupleRecord::kFormatVersion`. Until 1.0 a reader opens only files of its
own version, so there is no code for reading older layouts: a file written
before the change is refused with a message to convert it again.

## Tags

The standard tags in `kRAMTypedTags` (`inc/rntuple/RAMTags.h`) have columns of
their own in the record; every other tag stays SAM text in `tags`. The
converters split a record's tags with `RAMTagWriter`, and readers write them
back with `RAMFile::AppendTags()`, which keeps the input's order, so
`ramdump -h` still reproduces the input byte for byte. A typed tag goes into
its column only when writing the column back gives the same text; anything
else, `NM:i:007` or `NM:Z:x`, stays text.

A new typed tag needs a member in `RAMNTupleRecord`, an entry in
`kRAMTypedTags`, its handling in `RAMTagWriter::Split()` and
`RAMTagReader::Append()`, and a view in `RecordReader` (`tools/ramdump.cxx`),
which reads the record field by field and does not see new members by itself.

## PNEXT and TLEN

A record does not store PNEXT and TLEN as SAM has them, because other columns
already say most of what they hold. When the mate is on the record's own
reference, `pnext_delta` is PNEXT minus POS; otherwise it is PNEXT.
`tlen_delta` is TLEN minus `PredictTLEN()`, the distance between the 5' ends
of the pair as BWA and Picard compute it, with the mate's span taken to be the
record's own. Both are differences modulo 2^32, so every value comes back
exactly; for most records `tlen_delta` is 0. Write them with `SetMate()` and
read them with `GetPNEXT()` and `GetTLEN()`, which take
`RAMFile::MateOnSameRef()`. A program that reads the columns directly, such as
an RDataFrame, has to do the same.

## Read names

There is no QNAME column. The converters collect the names of every 10,000
records, the same blocks as the `-fqzcomp` qualities, and compress them with
tok3, htscodecs' name tokeniser at level 3 without arithmetic coding, as CRAM
3.1 does. The last row of a block carries the block in `nameblock`, and
METADATA's `block_ends` says where each block ends. `qname` stays a member of
`RAMNTupleRecord` for writing, but it is transient: a record read from a file
has an empty `GetQNAME()`. Read the names with `NameBlockReader`, which decodes
one block at a time and, while the rows come in order, decodes the next blocks
ahead on other threads. A region query decodes every block it touches, about
3 ms each.

## Quality scores

Lossless qualities go into the same blocks, in `qualblock`, and the record's
`qual` stays empty, with `kQualInBlock` set. The first byte of a block says how
it is compressed: by default with rANS Nx16 at order 1 with RLE and PACK, as
CRAM 3.1 compresses its quality series, the lengths of the qualities first and
then the qualities; with `-fqzcomp`, which passes `kFqzcomp` in the quality
policy, with fqzcomp. Qualities binned with `-illumina` or dropped with
`-dropqual` stay in `qual`, as does `*`. Read them with `QualityBlockReader`,
which returns the quality of any record, from its block or from `qual`.
