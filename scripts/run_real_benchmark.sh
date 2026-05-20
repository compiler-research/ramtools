#!/usr/bin/env bash
# Real benchmark driver for the README-style "Benchmark Results" section.
#
# Dataset (default): $SAM = complete chromosome-11 low-coverage alignment
# (1000 Genomes phase 3, GRCh37/hs37d5, coordinate-sorted, 101 bp). This is
# NOT the README's HG00154/72GB sample - results must be labeled honestly.
#
# Full "best-shot per format" matrix, on the SAME storage medium:
#   - RNTuple .root x4: LZMA, LZ4, ZLIB, ZSTD (samtoramntuple <codec>)
#   - TTree   .root x4: LZMA, LZ4, ZLIB, ZSTD (samtoram <codec>)
#   - BAM   (BGZF, the source download — single canonical config)
#   - CRAM  x2: 3.1 archive (smallest), 3.0 default (faster encode/decode)
#
# Each format is measured by region_query_benchmark (one binary, one run,
# one JSON) so the comparison is symmetric. Sizes are storage-independent;
# timings are warm-cache on the medium documented at the call site.
#
# ROOT must already be sourced (root-config on PATH), or set THISROOT to the
# path of a thisroot.sh and it will be sourced automatically.
#
# Overridable via environment:
#   SAM       input SAM           (default: <repo>/data/input.sam)
#   OUT       output dir for .root (default: dir of $SAM)
#   RESULTS   results dir         (default: <repo>/results_real)
#   REPS      benchmark reps      (default: 7)
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$SCRIPT_DIR/.." && pwd)
cd "$REPO"

if ! command -v root-config >/dev/null 2>&1; then
  if [ -n "${THISROOT:-}" ] && [ -f "$THISROOT" ]; then
    # shellcheck disable=SC1090
    source "$THISROOT"
  else
    echo "ERROR: ROOT not found. Source thisroot.sh or set THISROOT=/path/to/thisroot.sh" >&2
    exit 1
  fi
fi

SAM=${SAM:-"$REPO/data/input.sam"}
OUT=${OUT:-"$(dirname "$SAM")"}
RESULTS=${RESULTS:-"$REPO/results_real"}
REPS=${REPS:-7}
BUILD="$REPO/build"
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
mkdir -p "$RESULTS"

bytes() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }
log()   { echo "[$(date '+%H:%M:%S')] $*"; }

[ -f "$SAM" ] || { echo "ERROR: SAM not found: $SAM" >&2; exit 1; }
log "ROOT $(root-config --version); SAM=$SAM ($(bytes "$SAM") bytes)"

# --- TTree at all 4 codecs (idempotent) ---
for codec in lzma lz4 zlib zstd; do
  f="$OUT/ttree_${codec}.root"
  [ -f "$f" ] || { log "convert TTree $codec -> $f"; "$BUILD/tools/samtoram" "$SAM" "$f" "$codec"; }
done

# --- RNTuple at all 4 codecs (idempotent). Migrate legacy data/rn.root if present
#     (it was produced by the old samtoramntuple defaulting to ZSTD-5). ---
if [ -f "$OUT/rn.root" ] && [ ! -f "$OUT/rn_zstd.root" ]; then
  log "migrate $OUT/rn.root -> $OUT/rn_zstd.root (legacy default = ZSTD-5)"
  mv "$OUT/rn.root" "$OUT/rn_zstd.root"
fi
for codec in lzma lz4 zlib zstd; do
  f="$OUT/rn_${codec}.root"
  [ -f "$f" ] || { log "convert RNTuple $codec -> $f"; "$BUILD/tools/samtoramntuple" "$SAM" "$f" "$codec"; }
done

# --- BAM (source) + two CRAM presets (3.1 archive, 3.0 default) ---
BAM=${BAM:-"$OUT/input.bam"}
REF=${REF:-"$OUT/reference.fa"}
CRAM_ARCHIVE=${CRAM_ARCHIVE:-"$OUT/input.archive.cram"}
CRAM_V30=${CRAM_V30:-"$OUT/input.v30.cram"}
[ -f "$BAM" ]     || { echo "ERROR: BAM not found: $BAM" >&2; exit 1; }
[ -f "$BAM.bai" ] || { log "index BAM"; samtools index "$BAM"; }

# Migrate legacy data/input.cram (which was the 3.1 archive preset).
if [ -f "$OUT/input.cram" ] && [ ! -f "$CRAM_ARCHIVE" ]; then
  log "migrate $OUT/input.cram -> $CRAM_ARCHIVE (legacy = 3.1 archive)"
  mv "$OUT/input.cram" "$CRAM_ARCHIVE"
  [ -f "$OUT/input.cram.crai" ] && mv "$OUT/input.cram.crai" "$CRAM_ARCHIVE.crai"
fi

if [ ! -f "$CRAM_ARCHIVE" ]; then
  [ -f "$REF" ] || { echo "ERROR: reference FASTA missing: $REF" >&2; exit 1; }
  log "convert BAM -> CRAM 3.1 archive (small, slow)"
  samtools view -@ "$NCPU" -C -T "$REF" \
    --output-fmt-option version=3.1 \
    --output-fmt-option archive=1 \
    --output-fmt-option use_fqz=1 \
    --output-fmt-option use_tok=1 \
    --output-fmt-option use_arith=1 \
    -o "$CRAM_ARCHIVE" "$BAM"
fi
[ -f "$CRAM_ARCHIVE.crai" ] || { log "index CRAM archive"; samtools index "$CRAM_ARCHIVE"; }

if [ ! -f "$CRAM_V30" ]; then
  [ -f "$REF" ] || { echo "ERROR: reference FASTA missing: $REF" >&2; exit 1; }
  log "convert BAM -> CRAM 3.0 default (larger, faster)"
  samtools view -@ "$NCPU" -C -T "$REF" -o "$CRAM_V30" "$BAM"
fi
[ -f "$CRAM_V30.crai" ] || { log "index CRAM v30"; samtools index "$CRAM_V30"; }

# --- File sizes (storage-independent, fully trustworthy) ---
{
  echo "sam        $(bytes "$SAM")"
  echo "bam        $(bytes "$BAM")"
  echo "bam_bai    $(bytes "$BAM.bai")"
  echo "cram_archive $(bytes "$CRAM_ARCHIVE")"
  echo "cram_archive_crai $(bytes "$CRAM_ARCHIVE.crai")"
  echo "cram_v30   $(bytes "$CRAM_V30")"
  echo "cram_v30_crai $(bytes "$CRAM_V30.crai")"
  for codec in lzma lz4 zlib zstd; do
    echo "ttree_${codec} $(bytes "$OUT/ttree_${codec}.root")"
    echo "rntuple_${codec} $(bytes "$OUT/rn_${codec}.root")"
  done
} | tee "$RESULTS/file_sizes.txt"

# --- One unified gbench run, all 11 fixtures, REPS repetitions ---
log "region_query_benchmark (REPS=$REPS, all variants in one run)"
RAMBENCH_TTREE_LZMA="$OUT/ttree_lzma.root" \
RAMBENCH_TTREE_LZ4="$OUT/ttree_lz4.root" \
RAMBENCH_TTREE_ZLIB="$OUT/ttree_zlib.root" \
RAMBENCH_TTREE_ZSTD="$OUT/ttree_zstd.root" \
RAMBENCH_RNTUPLE_LZMA="$OUT/rn_lzma.root" \
RAMBENCH_RNTUPLE_LZ4="$OUT/rn_lz4.root" \
RAMBENCH_RNTUPLE_ZLIB="$OUT/rn_zlib.root" \
RAMBENCH_RNTUPLE_ZSTD="$OUT/rn_zstd.root" \
RAMBENCH_BAM="$BAM" \
RAMBENCH_CRAM_ARCHIVE="$CRAM_ARCHIVE" \
RAMBENCH_CRAM_V30="$CRAM_V30" \
RAMBENCH_REFERENCE="$REF" \
"$BUILD/benchmark/region_query_benchmark" \
    --benchmark_repetitions="$REPS" \
    --benchmark_report_aggregates_only=true \
    --benchmark_format=json \
    --benchmark_out="$RESULTS/query_all.json" \
    2>"$RESULTS/query_all.err" | tee "$RESULTS/query_all.console" >/dev/null

log "DONE - artifacts in $RESULTS/"
