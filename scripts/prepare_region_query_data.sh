#!/usr/bin/env bash
# Prepare the data the region_query_benchmark needs: a SAM from the 1000 Genomes
# HG00097 chromosome-11 sample (same as .github/workflows/benchmark.yml),
# converted to both the TTree and RNTuple RAM formats.
#
# Idempotent: any artifact that already exists is reused. Delete the files in
# DATA_DIR to force a rebuild.
#
# Requires `curl` and `samtools` on PATH, and ROOT already sourced
# (the converters are ROOT executables and need ROOT libs at runtime).
#
# Env overrides:
#   BUILD_DIR       build directory containing tools/ (default: build)
#   DATA_DIR        where to place data files          (default: data)
#   ALIGNMENT_URL   BAM/SAM/SAM.gz URL to benchmark on
set -euo pipefail

BUILD_DIR="${BUILD_DIR:-build}"
DATA_DIR="${DATA_DIR:-data}"
ALIGNMENT_URL="${ALIGNMENT_URL:-http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00097/alignment/HG00097.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam}"

BAM="${DATA_DIR}/input.bam"
SAM="${DATA_DIR}/input.sam"
TTREE="${DATA_DIR}/input_ttree.root"
RNTUPLE="${DATA_DIR}/input_rntuple.root"

mkdir -p "${DATA_DIR}"

if [[ -f "${TTREE}" && -f "${RNTUPLE}" ]]; then
  echo "[prepare] ${TTREE} and ${RNTUPLE} already exist - nothing to do."
  exit 0
fi

for bin in curl samtools; do
  command -v "${bin}" >/dev/null 2>&1 || { echo "[prepare] ERROR: '${bin}' not found on PATH" >&2; exit 1; }
done
for tool in samtoram samtoramntuple; do
  [[ -x "${BUILD_DIR}/tools/${tool}" ]] || { echo "[prepare] ERROR: ${BUILD_DIR}/tools/${tool} not built" >&2; exit 1; }
done

if [[ ! -f "${SAM}" ]]; then
  if [[ ! -f "${BAM}" ]]; then
    echo "[prepare] downloading ${ALIGNMENT_URL}"
    curl -fL "${ALIGNMENT_URL}" -o "${BAM}"
  fi
  echo "[prepare] samtools view -> ${SAM}"
  samtools view -h -o "${SAM}" "${BAM}"
fi

if [[ ! -f "${TTREE}" ]]; then
  echo "[prepare] converting SAM -> TTree (${TTREE})"
  "${BUILD_DIR}/tools/samtoram" "${SAM}" "${TTREE}"
fi

if [[ ! -f "${RNTUPLE}" ]]; then
  echo "[prepare] converting SAM -> RNTuple (${RNTUPLE})"
  "${BUILD_DIR}/tools/samtoramntuple" "${SAM}" "${RNTUPLE}"
fi

echo "[prepare] done: ${TTREE} ${RNTUPLE}"
