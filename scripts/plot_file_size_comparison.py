"""Generate the "File Size Comparison" chart for the Benchmark Results section.

These are REAL measured sizes from results_real/file_sizes.txt, produced by
scripts/run_real_benchmark.sh on the complete HG00097 chromosome-11
low-coverage 1000 Genomes sample (12,876,510 reads, 5.36 GiB SAM).
NOT the README's HG00154/72GB figures.

Usage:
    python3 scripts/plot_file_size_comparison.py [-o OUTPUT.png] [--sizes file_sizes.txt]
"""
import argparse
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Map file_sizes.txt keys -> chart labels, in display order.
ORDER = [
    ("sam", "SAM"),
    ("bam", "BAM"),
    ("rntuple", "RNTuple ROOT"),
    ("ttree_lz4", "TTree ROOT(LZ4)"),
    ("ttree_lzma", "TTree ROOT(LZMA)"),
    ("ttree_zlib", "TTree ROOT(ZLIB)"),
]
BAR_COLOR = "#69bfba"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="img/file-size-comparison.png")
    ap.add_argument("--sizes", default="results_real/file_sizes.txt")
    args = ap.parse_args()

    raw = dict(
        line.split() for line in pathlib.Path(args.sizes).read_text().splitlines() if line.strip()
    )
    labels = [lbl for key, lbl in ORDER]
    values = [int(raw[key]) / 1e9 for key, _ in ORDER]  # bytes -> GB (decimal, like README)

    fig, ax = plt.subplots(figsize=(9.5, 5.0), dpi=100)
    bars = ax.bar(labels, values, color=BAR_COLOR, width=0.55)

    ax.set_xlabel("File format")
    ax.set_ylabel("File Sizes(GB)")
    ax.set_ylim(0, max(values) * 1.18)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#e6e6e6")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "HG00097 chr11 low-coverage (1000G phase3) — 12,876,510 reads — REAL measured sizes",
        fontsize=9,
    )

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(values) * 0.012,
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"wrote {args.output} from {args.sizes}")


if __name__ == "__main__":
    main()
