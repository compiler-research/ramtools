import argparse
import json
import pathlib
import statistics


def fmt_size(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.2f} {unit}"
        n /= 1024


def fmt_dist(samples):
    if not samples:
        return "n/a"
    mean = statistics.mean(samples)
    stddev = statistics.stdev(samples) if len(samples) > 1 else 0.0
    if mean < 1e-3:
        unit, scale = "μs", 1e6
    elif mean < 1:
        unit, scale = "ms", 1e3
    else:
        unit, scale = "s", 1.0
    if stddev > 0:
        return f"{mean * scale:.3g} ± {stddev * scale:.2g} {unit}"
    return f"{mean * scale:.3g} {unit}"


def ratio(base, pr):
    if not base or not pr:
        return "-"
    pm = statistics.mean(pr)
    return f"{statistics.mean(base) / pm:.2f}" if pm else "-"


def region_label(q):
    if q.get("records") is not None:
        return f"`{q['region']}` ({q['records']:,} reads)"
    return f"`{q['region']}`"


def records_line(pr):
    sam = pr.get("sam_records")
    ram = pr["ram"].get("records")
    cram = pr.get("cram", {}).get("records")
    parts = [f"SAM `{sam:,}`" if sam is not None else "SAM `?`"]
    if ram is not None:
        tag = "match" if sam is not None and ram == sam else "MISMATCH"
        parts.append(f"RAM `{ram:,}` ({tag})")
    if cram is not None:
        tag = "match" if sam is not None and cram == sam else "MISMATCH"
        parts.append(f"CRAM `{cram:,}` ({tag})")
    return f"**Records:** {' / '.join(parts)}"


def render_comparison(pr, base, base_sha, pr_sha):
    bl = f"base ({base_sha[:7]})" if base_sha else "base"
    pl = f"PR ({pr_sha[:7]})" if pr_sha else "PR"
    lines = [
        f"Comparing **{bl}** vs **{pl}** on `{fmt_size(pr['sam_size_bytes'])}` SAM input",
        "",
        f"| Benchmark | {bl} | {pl} | base / PR |",
        "|-----------|------|------|----------:|",
        f"| RAM size | {fmt_size(base['ram']['size_bytes'])} | {fmt_size(pr['ram']['size_bytes'])} "
        f"| {base['ram']['size_bytes'] / pr['ram']['size_bytes']:.2f} |",
        f"| RAM conversion | {fmt_dist(base['ram']['conversion_seconds'])} "
        f"| {fmt_dist(pr['ram']['conversion_seconds'])} "
        f"| {ratio(base['ram']['conversion_seconds'], pr['ram']['conversion_seconds'])} |",
    ]
    for b, p in zip(base["ram"]["queries"], pr["ram"]["queries"]):
        lines.append(
            f"| query {region_label(p)} | {fmt_dist(b['seconds'])} | {fmt_dist(p['seconds'])} "
            f"| {ratio(b['seconds'], p['seconds'])} |"
        )
    return lines + [""]


def render_absolute(pr, pr_sha):
    pl = f"PR ({pr_sha[:7]})" if pr_sha else "PR"
    lines = [
        f"**{pl}** on `{fmt_size(pr['sam_size_bytes'])}` SAM input",
        "",
        "| Benchmark | Value |",
        "|-----------|-------|",
        f"| RAM size | {fmt_size(pr['ram']['size_bytes'])} |",
        f"| RAM conversion | {fmt_dist(pr['ram']['conversion_seconds'])} |",
    ]
    for q in pr["ram"]["queries"]:
        lines.append(f"| query {region_label(q)} | {fmt_dist(q['seconds'])} |")
    return lines + [""]


def render_cram(pr):
    cram, sam = pr["cram"], pr["sam_size_bytes"]
    lines = [
        "### CRAM 3.1 reference (archive preset: fqzcomp + name tokenizer + rANS)",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| CRAM size | {fmt_size(cram['size_bytes'])} ({cram['size_bytes'] / sam:.1%} of SAM) |",
        f"| CRAM conversion | {fmt_dist(cram['conversion_seconds'])} |",
        f"| CRAM index | {fmt_dist(cram['index_seconds'])} |",
    ]
    for q in cram["queries"]:
        lines.append(f"| query {region_label(q)} (CRAM) | {fmt_dist(q['seconds'])} |")
    return lines + [""]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pr_json")
    ap.add_argument("--base")
    ap.add_argument("--base-sha", default="")
    ap.add_argument("--pr-sha", default="")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    pr = json.loads(pathlib.Path(args.pr_json).read_text())
    base = json.loads(pathlib.Path(args.base).read_text()) if args.base else None

    out = ["## Benchmark results", "", records_line(pr), ""]
    out += render_comparison(pr, base, args.base_sha, args.pr_sha) if base else render_absolute(pr, args.pr_sha)
    if "cram" in pr:
        out += render_cram(pr)

    text = "\n".join(out)
    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(args.output).write_text(text)
    print(text)


if __name__ == "__main__":
    main()
