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
| `quality_blocks_test` | `-fqzcomp` quality blocks: which qualities go into a block, reading them back in any order, and damaged block ends |
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

## One rule about the shared state

`RAMNTupleRecord` keeps the reference-name tables, the longest span, the sort
flag, the quality block ends and the SAM header as static members that
describe the file currently open or being written. `InitializeRefs()` resets
them, and only the converters and `OpenRAMFile()` may call it. RNTuple
constructs a record whenever a view or a writer model is created, so the
constructor must never touch that state; a reset there once made an unsorted
file look sorted and dropped records from queries. If you add per-file state,
reset it in `InitializeRefs()` and extend the test
`ConstructingARecordKeepsTheOpenFileState`.

## Changing the file format

Any change to what a RAM file holds, in the records or in METADATA, bumps
`RAMNTupleRecord::kFormatVersion`. Until 1.0 a reader opens only files of its
own version, so there is no code for reading older layouts: a file written
before the change is refused with a message to convert it again.
