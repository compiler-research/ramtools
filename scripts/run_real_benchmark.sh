#!/usr/bin/env bash
# Real benchmark driver: converts a SAM to every supported format and runs
# region_query_benchmark against all of them. Thin wrapper around
# scripts/convert_all_formats.sh (single source of truth for conversions) +
# the unified gbench binary.
#
# Dataset (default): chromosome-11 low-coverage 1000G phase-3 sample at
# data/input.sam. Override SAM= for any other input.
#
# ROOT must already be sourced (root-config on PATH), or set THISROOT to a
# thisroot.sh and it will be sourced automatically.
#
# Env overrides:
#   SAM       input SAM                  (default: <repo>/data/input.sam)
#   OUT       conversion output dir      (default: dir of SAM)
#   BAM       BAM path                   (default: $OUT/input.bam)
#   REF       reference FASTA            (default: $OUT/reference.fa)
#   RESULTS   results dir for JSON/txt   (default: <repo>/results_real)
#   REPS      gbench repetitions         (default: 7)
#   CODECS    ROOT codecs to produce     (default: "lzma lz4 zlib zstd")
#   CRAM_PRESETS  CRAM presets           (default: "archive v30")
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
BAM=${BAM:-"$OUT/input.bam"}
REF=${REF:-"$OUT/reference.fa"}
RESULTS=${RESULTS:-"$REPO/results_real"}
REPS=${REPS:-7}
CODECS=${CODECS:-"lzma lz4 zlib zstd"}
CRAM_PRESETS=${CRAM_PRESETS:-"archive v30"}
BUILD="$REPO/build"
mkdir -p "$RESULTS"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

[ -f "$SAM" ] || { echo "ERROR: SAM not found: $SAM" >&2; exit 1; }
log "ROOT $(root-config --version); SAM=$SAM; OUT=$OUT; RESULTS=$RESULTS; REPS=$REPS"

# --- One source of truth for conversions: convert_all_formats.sh ---
SAM="$SAM" OUT="$OUT" BAM="$BAM" REF="$REF" \
  CODECS="$CODECS" CRAM_PRESETS="$CRAM_PRESETS" \
  SIZES="$RESULTS/file_sizes.txt" \
  BUILD="$BUILD" \
  bash "$SCRIPT_DIR/convert_all_formats.sh" "$SAM"

# --- Unified gbench run: all 11 fixtures, REPS repetitions, single JSON ---
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
RAMBENCH_CRAM_ARCHIVE="$OUT/input.archive.cram" \
RAMBENCH_CRAM_V30="$OUT/input.v30.cram" \
RAMBENCH_REFERENCE="$REF" \
"$BUILD/benchmark/region_query_benchmark" \
    --benchmark_repetitions="$REPS" \
    --benchmark_report_aggregates_only=true \
    --benchmark_format=json \
    --benchmark_out="$RESULTS/query_all.json" \
    2>"$RESULTS/query_all.err" | tee "$RESULTS/query_all.console" >/dev/null

log "DONE - artifacts in $RESULTS/"
