"""Grouped bar chart of REAL region-query times: RNTuple vs TTree(LZMA/LZ4/ZLIB).

Reads results_real/query_{lzma,lz4,zlib}.json (Google Benchmark output from
scripts/run_real_benchmark.sh). Log y-axis because times span 0.04 s .. 87 s.

Usage: python3 scripts/plot_query_comparison.py [-o OUT.png] [--results DIR]
"""
import argparse
import json
import pathlib
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REGIONS = ["Small\n~100 bp", "1 Mb", "10 Mb", "full chr11\n~135 Mb"]
COLORS = {"RNTuple": "#2a9d8f", "BAM (samtools)": "#264653",
          "TTree LZMA": "#e76f51", "TTree LZ4": "#f4a261", "TTree ZLIB": "#e9c46a"}


def mean_times(results, codec, fmt):
    d = json.loads((results / f"query_{codec}.json").read_text())
    out = {}
    for b in d["benchmarks"]:
        if b.get("aggregate_name") != "mean":
            continue
        nm = b["name"]  # RegionQueryFixture/<fmt>/<idx>_mean
        if f"/{fmt}/" not in nm:
            continue
        idx = int(nm.split("/")[2].split("_")[0])
        out[idx] = b["real_time"]
    return [out[i] for i in range(4)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="img/query-time-comparison.png")
    ap.add_argument("--results", default="results_real")
    args = ap.parse_args()
    R = pathlib.Path(args.results)

    bam = json.loads((R / "query_bam.json").read_text())["by_idx"]
    series = {
        "TTree LZMA": mean_times(R, "lzma", "TTree"),
        "TTree LZ4": mean_times(R, "lz4", "TTree"),
        "TTree ZLIB": mean_times(R, "zlib", "TTree"),
        "BAM (samtools)": [bam[str(i)]["mean"] for i in range(4)],
        # RNTuple file is identical across runs; average the 3 as the estimate.
        "RNTuple": [statistics.mean(v) for v in zip(
            mean_times(R, "lzma", "RNTuple"),
            mean_times(R, "lz4", "RNTuple"),
            mean_times(R, "zlib", "RNTuple"))],
    }
    order = ["RNTuple", "BAM (samtools)", "TTree LZ4", "TTree ZLIB", "TTree LZMA"]

    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=100)
    n = len(order)
    width = 0.8 / n
    x = range(len(REGIONS))
    for k, name in enumerate(order):
        offs = [i + (k - (n - 1) / 2) * width for i in x]
        bars = ax.bar(offs, series[name], width, label=name, color=COLORS[name])
        for b, v in zip(bars, series[name]):
            ax.text(b.get_x() + b.get_width() / 2, v * 1.05,
                    f"{v:.3g}", ha="center", va="bottom", fontsize=7.5, rotation=90)

    ax.set_yscale("log")
    ax.set_ylabel("Query time (s, log scale) — lower is better")
    ax.set_xlabel("Region size")
    ax.set_xticks(list(x))
    ax.set_xticklabels(REGIONS)
    ax.set_ylim(0.02, 200)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#e6e6e6")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Region query time — RNTuple vs BAM(samtools) vs TTree codecs\n"
                 "REAL: HG00097 chr11 low-coverage, 12,876,510 reads "
                 "(ROOT viewers 3 reps, BAM 7 reps)",
                 fontsize=9)
    ax.legend(frameon=False, ncol=4, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
