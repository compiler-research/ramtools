import argparse
import json
import os
import pathlib
import re
import subprocess
import time

SAM = "data/input.sam"
REF = "data/reference.fa"
CRAM = "results/input.cram"
QUERY_REPEAT = 3


def measure(cmd, n, cleanup=None):
    samples, last = [], ""
    for _ in range(n):
        if cleanup:
            cleanup()
        start = time.perf_counter()
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
        samples.append(time.perf_counter() - start)
        last = proc.stdout
    return samples, last


def parse_records(output):
    m = re.findall(r"Found\s+(\d+)\s+records", output)
    return int(m[-1]) if m else None


def count_int(cmd):
    return int(subprocess.run(cmd, text=True, capture_output=True, check=True).stdout.strip())


def count_ram(cmd):
    return parse_records(subprocess.run(cmd, text=True, capture_output=True, check=True).stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default="build")
    ap.add_argument("--ram", default="results/input.ram")
    ap.add_argument("--no-cram", action="store_true")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    threads = str(os.cpu_count() or 2)
    regions = os.environ.get("REGIONS", "").split()
    ram = pathlib.Path(args.ram)
    cram = pathlib.Path(CRAM)
    ram.parent.mkdir(parents=True, exist_ok=True)

    sam_records = count_int(["samtools", "view", "-c", SAM])

    convert, _ = measure(
        [f"{args.build_dir}/tools/samtoramntuple", SAM, str(ram)],
        1,
        cleanup=lambda: ram.unlink(missing_ok=True),
    )
    ram_total = count_ram([f"{args.build_dir}/tools/ramntupleview", str(ram)])
    ram_queries = []
    for r in regions:
        samples, out = measure([f"{args.build_dir}/tools/ramntupleview", str(ram), r], QUERY_REPEAT)
        ram_queries.append({"region": r, "seconds": samples, "records": parse_records(out)})

    result = {
        "sam_size_bytes": os.path.getsize(SAM),
        "sam_records": sam_records,
        "ram": {
            "size_bytes": ram.stat().st_size,
            "records": ram_total,
            "conversion_seconds": convert,
            "queries": ram_queries,
        },
    }

    if not args.no_cram:
        def rm_cram():
            cram.unlink(missing_ok=True)
            pathlib.Path(str(cram) + ".crai").unlink(missing_ok=True)

        cram_convert, _ = measure(
            [
                "samtools", "view", "-@", threads, "-C", "-T", REF,
                "--output-fmt-option", "version=3.1",
                "--output-fmt-option", "archive=1",
                "--output-fmt-option", "use_fqz=1",
                "--output-fmt-option", "use_tok=1",
                "--output-fmt-option", "use_arith=1",
                "-o", str(cram), SAM,
            ],
            1,
            cleanup=rm_cram,
        )
        cram_index, _ = measure(
            ["samtools", "index", "-@", threads, str(cram)],
            1,
            cleanup=lambda: pathlib.Path(str(cram) + ".crai").unlink(missing_ok=True),
        )
        cram_total = count_int(["samtools", "view", "-@", threads, "-c", str(cram)])
        cram_queries = []
        for r in regions:
            samples, out = measure(
                ["samtools", "view", "-@", threads, "-c", "-F", "2308", str(cram), r], QUERY_REPEAT
            )
            cram_queries.append({"region": r, "seconds": samples, "records": int(out.strip())})
        result["cram"] = {
            "size_bytes": cram.stat().st_size,
            "records": cram_total,
            "conversion_seconds": cram_convert,
            "index_seconds": cram_index,
            "queries": cram_queries,
        }

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
