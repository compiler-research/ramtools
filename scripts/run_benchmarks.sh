#!/usr/bin/env bash
#
# Run the full RAMtools benchmark suite, forwarding any dataset flags to every binary
# and writing one JSON result file per benchmark into the build tree.
#
# Usage:
#   scripts/run_benchmarks.sh [BUILD_DIR] [-- <extra benchmark flags...>]
#
# Examples:
#   scripts/run_benchmarks.sh                       # build/, synthetic data
#   scripts/run_benchmarks.sh build
#   scripts/run_benchmarks.sh build -- --sam=/data/HG002.sam --regions=all
#
set -euo pipefail

BUILD_DIR="${1:-build}"
shift || true
if [[ "${1:-}" == "--" ]]; then
   shift
fi
EXTRA_ARGS=("$@")

BENCH_DIR="${BUILD_DIR}/benchmark"
if [[ ! -d "${BENCH_DIR}" ]]; then
   echo "error: ${BENCH_DIR} not found - build first: cmake --build ${BUILD_DIR}" >&2
   exit 1
fi

# Absolute path so JSON lands in the build tree and temp files stay out of the source tree.
ABS_BENCH_DIR="$(cd "${BENCH_DIR}" && pwd)"
cd "${ABS_BENCH_DIR}"

BINARIES=(
   sam_to_ram_benchmark
   conversion_time_benchmark
   region_query_benchmark
   chromosome_split_benchmark
   columnar_read_benchmark
)

for b in "${BINARIES[@]}"; do
   if [[ ! -x "./${b}" ]]; then
      echo "skip: ./${b} not built" >&2
      continue
   fi
   echo "=== ${b} ==="
   "./${b}" "${EXTRA_ARGS[@]}" \
      --benchmark_out="${ABS_BENCH_DIR}/${b}.json" \
      --benchmark_out_format=json
done

echo "JSON results written to ${ABS_BENCH_DIR}"
