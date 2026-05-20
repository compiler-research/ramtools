#!/usr/bin/env bash
# Convert one SAM file to every format ramtools compares against:
#   RNTuple .root x4 (LZMA, LZ4, ZLIB, ZSTD)
#   TTree   .root x4 (LZMA, LZ4, ZLIB, ZSTD)
#   BAM            (BGZF default; generated from SAM if not present)
#   CRAM           (3.1 archive preset + 3.0 default preset; needs $REF)
#
# Idempotent: any output that already exists is left alone. Writes a
# sizes summary to $SIZES (default: $OUT/file_sizes.txt).
#
# Does NOT run benchmarks. Pair with scripts/run_real_benchmark.sh (or
# call the gbench binary directly) for the timing measurement.
#
# Usage:
#   scripts/convert_all_formats.sh <input.sam>
#   SAM=path OUT=dir [BAM=path] [REF=path] scripts/convert_all_formats.sh
#
# ROOT must be sourced (root-config on PATH), or set THISROOT=/path/to/thisroot.sh.
#
# Env overrides:
#   SAM            input SAM                  (positional arg or env)
#   OUT            output directory           (default: dir of SAM)
#   BAM            BAM path                   (default: $OUT/input.bam)
#   REF            reference FASTA            (default: $OUT/reference.fa; CRAM skipped if absent)
#   SIZES          sizes summary file         (default: $OUT/file_sizes.txt)
#   CODECS         ROOT codecs to produce     (default: "lzma lz4 zlib zstd")
#   CRAM_PRESETS   CRAM presets to produce    (default: "archive v30")
#   BUILD          ramtools build directory   (default: <repo>/build)
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$SCRIPT_DIR/.." && pwd)

if ! command -v root-config >/dev/null 2>&1; then
  if [ -n "${THISROOT:-}" ] && [ -f "$THISROOT" ]; then
    # shellcheck disable=SC1090
    source "$THISROOT"
  else
    echo "ERROR: ROOT not found. Source thisroot.sh or set THISROOT=/path/to/thisroot.sh" >&2
    exit 1
  fi
fi

SAM=${1:-${SAM:-}}
[ -n "$SAM" ] || { echo "Usage: $0 <input.sam>  (or set SAM= and other env vars; see header)" >&2; exit 1; }
[ -f "$SAM" ] || { echo "ERROR: SAM not found: $SAM" >&2; exit 1; }
SAM=$(cd "$(dirname "$SAM")" && pwd)/$(basename "$SAM")

OUT=${OUT:-"$(dirname "$SAM")"}
BAM=${BAM:-"$OUT/input.bam"}
REF=${REF:-"$OUT/reference.fa"}
SIZES=${SIZES:-"$OUT/file_sizes.txt"}
CODECS=${CODECS:-"lzma lz4 zlib zstd"}
CRAM_PRESETS=${CRAM_PRESETS:-"archive v30"}
BUILD=${BUILD:-"$REPO/build"}
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)

bytes() { stat -f%z "$1" 2>/dev/null || stat -c%s "$1"; }
log()   { echo "[$(date '+%H:%M:%S')] $*"; }
have()  { command -v "$1" >/dev/null 2>&1; }

for t in "$BUILD/tools/samtoram" "$BUILD/tools/samtoramntuple"; do
  [ -x "$t" ] || { echo "ERROR: $t not built (cmake --build $BUILD --target samtoram samtoramntuple)" >&2; exit 1; }
done
have samtools || echo "WARN: samtools not on PATH - BAM/CRAM steps will be skipped" >&2

mkdir -p "$OUT"
log "SAM     : $SAM ($(bytes "$SAM") bytes)"
log "OUT     : $OUT"
log "CODECS  : $CODECS"
log "CRAM    : $CRAM_PRESETS"

# --- TTree (one .root per ROOT compression codec) ---
for codec in $CODECS; do
  f="$OUT/ttree_${codec}.root"
  if [ -f "$f" ]; then log "TTree $codec: cached ($f)"; continue; fi
  log "TTree $codec: convert -> $f"
  "$BUILD/tools/samtoram" "$SAM" "$f" "$codec"
done

# --- RNTuple (one .root per ROOT compression codec) ---
for codec in $CODECS; do
  f="$OUT/rn_${codec}.root"
  if [ -f "$f" ]; then log "RNTuple $codec: cached ($f)"; continue; fi
  log "RNTuple $codec: convert -> $f"
  "$BUILD/tools/samtoramntuple" "$SAM" "$f" "$codec"
done

# --- BAM (generate from SAM if missing) + .bai ---
if have samtools; then
  if [ ! -f "$BAM" ]; then
    log "BAM: convert from SAM -> $BAM"
    samtools view -@ "$NCPU" -b -o "$BAM" "$SAM"
  else
    log "BAM: cached ($BAM)"
  fi
  [ -f "$BAM.bai" ] || { log "BAM: index"; samtools index "$BAM"; }
fi

# --- CRAM (skip if no reference or no samtools) ---
if ! have samtools; then
  :
elif [ ! -f "$REF" ]; then
  echo "WARN: reference FASTA missing: $REF -- skipping CRAM (set REF= to enable)" >&2
elif [ ! -f "$BAM" ]; then
  echo "WARN: BAM missing -- skipping CRAM" >&2
else
  for preset in $CRAM_PRESETS; do
    cf="$OUT/input.${preset}.cram"
    if [ -f "$cf" ]; then log "CRAM $preset: cached ($cf)"; else
      case "$preset" in
        archive)
          log "CRAM 3.1 archive (smallest, slow): convert -> $cf"
          samtools view -@ "$NCPU" -C -T "$REF" \
            --output-fmt-option version=3.1 \
            --output-fmt-option archive=1 \
            --output-fmt-option use_fqz=1 \
            --output-fmt-option use_tok=1 \
            --output-fmt-option use_arith=1 \
            -o "$cf" "$BAM"
          ;;
        v30)
          log "CRAM 3.0 default (faster, larger): convert -> $cf"
          samtools view -@ "$NCPU" -C -T "$REF" -o "$cf" "$BAM"
          ;;
        *)
          echo "WARN: unknown CRAM preset '$preset'; skipping" >&2
          continue
          ;;
      esac
    fi
    [ -f "$cf.crai" ] || samtools index "$cf"
  done
fi

# --- Sizes summary ---
{
  echo "sam        $(bytes "$SAM")"
  [ -f "$BAM" ]      && echo "bam        $(bytes "$BAM")"
  [ -f "$BAM.bai" ]  && echo "bam_bai    $(bytes "$BAM.bai")"
  for codec in $CODECS; do
    [ -f "$OUT/ttree_${codec}.root" ] && echo "ttree_${codec} $(bytes "$OUT/ttree_${codec}.root")"
    [ -f "$OUT/rn_${codec}.root" ]    && echo "rntuple_${codec} $(bytes "$OUT/rn_${codec}.root")"
  done
  for preset in $CRAM_PRESETS; do
    cf="$OUT/input.${preset}.cram"
    [ -f "$cf" ]      && echo "cram_${preset}      $(bytes "$cf")"
    [ -f "$cf.crai" ] && echo "cram_${preset}_crai $(bytes "$cf.crai")"
  done
} | tee "$SIZES"

log "DONE -> $SIZES"
