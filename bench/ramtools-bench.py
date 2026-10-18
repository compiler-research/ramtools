#!/usr/bin/env python3
"""Reproduce RAMTools' benchmarks on any Linux or macOS machine.

The data is chr20 of NA12878 from the 1000 Genomes 30x high-coverage set (NYGC,
GRCh38, NovaSeq, 17,705,654 records), read from EBI. Two git revisions of
ramtools are built the same way and measured in alternating order: conversion
from SAM and BAM, file size by field, the round trip, full and region dumps,
counts and statistics, with samtools on the same files as the yardstick.

  bench/ramtools-bench.py check                    # what is installed and what is missing
  bench/ramtools-bench.py fetch                    # download and prepare the data, about 15 GB
  bench/ramtools-bench.py run --base develop --head HEAD
  bench/ramtools-bench.py run --base develop --head HEAD --input chr20
  bench/ramtools-bench.py report RESULT.json       # print a saved result again

bench/README.md describes the data, the protocol and how to read the results.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Bump when what is measured, or how, changes; results of different protocols
# are not comparable.
PROTOCOL = 1

CRAM_URL = "https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram"
REF_URL = ("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
           "GRCh38_full_analysis_set_plus_decoy_hla.fa")
CHROM = "chr20"
CHROM_M5 = "b18e6c531b0bd70e949a7fc20859cb01"  # M5 of chr20 in the CRAM's header
RECORDS = 17_705_654
RECORDS_MD5 = "89386f84300ce2a3b5100ec5760a197f"  # md5 of `samtools view` of the chr20 BAM, no header

STEM = "NA12878.chr20"
SUBSETS = {"subset": 2_000_000, "small": 50_000}  # the first records of chr20
INPUTS = ["small", "subset", "chr20"]

# Regions on all of chr20, as in the numbers quoted in CHANGELOG.md.
CHR20_REGIONS = {"large": "chr20:30000000-40000000",
                 "small": ["chr20:30000000-30001000", "chr20:45000000-45000500", "chr20:10000000-10002000"]}

TOOLS = ["samtoramntuple", "bamtoramntuple", "ramdump", "ramstats"]
# What a WORKTREE build is made from; a change anywhere else does not rebuild it.
SOURCES = ["CMakeLists.txt", "cmake", "inc", "src", "tools", "test", "benchmark"]

# Timed runs of each step per input; 0 leaves the step out.
REPEATS = {
    "small": {"convert_sam": 3, "convert_threads": 3, "convert_bam": 3, "read": 3, "fast": 5},
    "subset": {"convert_sam": 5, "convert_threads": 5, "convert_bam": 5, "read": 5, "fast": 5},
    "chr20": {"convert_sam": 0, "convert_threads": 2, "convert_bam": 1, "read": 3, "fast": 3},
}


# ---------------------------------------------------------------- helpers

def say(text):
    print(text, flush=True)


def sh(cmd, cwd=None):
    """stdout of a command, or "" when it cannot run."""
    try:
        return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, check=True,
                              cwd=cwd).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def check_call(cmd, cwd=None, log=None):
    """Runs a command and stops with its output when it fails."""
    cmd = [str(c) for c in cmd]
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if log:
        with open(log, "a") as f:
            f.write("$ " + " ".join(cmd) + "\n" + p.stdout)
    if p.returncode != 0:
        tail = "\n".join(p.stdout.splitlines()[-20:])
        sys.exit(f"failed ({p.returncode}): {' '.join(cmd)}\n{tail}")
    return p.stdout


def md5_of_output(cmd):
    """md5 of what a command writes to stdout; stops when it fails."""
    cmd = [str(c) for c in cmd]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    h = hashlib.md5()
    for chunk in iter(lambda: proc.stdout.read(1 << 20), b""):
        h.update(chunk)
    err = proc.stderr.read().decode(errors="replace")
    if proc.wait() != 0:
        sys.exit(f"failed ({proc.returncode}): {' '.join(cmd)}\n{err[-2000:]}")
    return h.hexdigest()


def md5_of_file(path):
    """md5 of a file, cached next to it while its size and time do not change."""
    path = Path(path)
    st = path.stat()
    cache = path.with_name(path.name + ".md5")
    key = f"{st.st_size} {int(st.st_mtime)}"
    if cache.exists():
        size_time, _, digest = cache.read_text().strip().rpartition(" ")
        if size_time == key:
            return digest
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    cache.write_text(f"{key} {h.hexdigest()}\n")
    return h.hexdigest()


def git(*args, check=True):
    p = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True)
    if p.returncode != 0:
        if check:
            sys.exit(f"git {' '.join(args)}: {p.stderr.strip()}")
        return None
    return p.stdout.strip()


def version_tuple(text):
    m = re.search(r"(\d+)\.(\d+)(?:[./](\d+))?", text or "")
    return tuple(int(x or 0) for x in m.groups()) if m else (0, 0, 0)


# ---------------------------------------------------------------- environment

def environment(pin):
    """What the numbers are measured on."""
    env = {"python": platform.python_version(), "root": sh(["root-config", "--version"]),
           "samtools": (sh(["samtools", "--version"]).splitlines() or [""])[0],
           "pinned": " ".join(pin) or None}
    if sys.platform == "darwin":
        env.update(cpu=sh(["sysctl", "-n", "machdep.cpu.brand_string"]), model=sh(["sysctl", "-n", "hw.model"]),
                   sockets=1, cores=int(sh(["sysctl", "-n", "hw.physicalcpu"]) or 0),
                   threads=int(sh(["sysctl", "-n", "hw.logicalcpu"]) or 0),
                   memory_gb=round(int(sh(["sysctl", "-n", "hw.memsize"]) or 0) / 1e9, 1),
                   os="macOS " + sh(["sw_vers", "-productVersion"]), memory_metric="peak physical footprint")
    else:
        lscpu = {}
        for line in sh(["lscpu"]).splitlines():
            key, _, value = line.partition(":")
            lscpu[key.strip()] = value.strip()
        sockets = int(lscpu.get("Socket(s)", "1") or 1)
        per_socket = int(lscpu.get("Core(s) per socket", "0") or 0)
        mem = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        governor = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
        env.update(cpu=lscpu.get("Model name", platform.processor()), model=platform.machine(), sockets=sockets,
                   cores=sockets * per_socket or os.cpu_count() or 0, threads=os.cpu_count() or 0,
                   memory_gb=round(mem / 1e9, 1), os=f"{platform.system()} {platform.release()}",
                   distribution=_distribution(), memory_metric="peak RSS",
                   governor=governor.read_text().strip() if governor.exists() else None)
    env["machine_id"] = hashlib.md5(f"{env['cpu']}|{env['cores']}|{env['memory_gb']}|{env['os']}".encode()
                                    ).hexdigest()[:8]
    return env


def _distribution():
    try:
        for line in Path("/etc/os-release").read_text().splitlines():
            if line.startswith("PRETTY_NAME="):
                return line.split("=", 1)[1].strip('"')
    except OSError:
        pass
    return None


def check(args):
    """Reports what is installed; True when everything needed is there."""
    ok = True

    def row(name, found, need, hint=""):
        nonlocal ok
        good = bool(found) and (need is None or version_tuple(found) >= need)
        ok = ok and good
        want = f" (need {'.'.join(map(str, need))} or newer)" if need else ""
        say(f"  {'ok     ' if good else 'MISSING'} {name:10} {found or 'not found'}{want}"
            f"{'' if good or not hint else '  ' + hint}")

    say("prerequisites:")
    row("python", platform.python_version(), (3, 8))
    row("git", sh(["git", "--version"]), None)
    row("cmake", (sh(["cmake", "--version"]).splitlines() or [""])[0], (3, 16))
    row("c++", (sh(["c++", "--version"]).splitlines() or [""])[0], None, "install a C++17 compiler")
    row("ROOT", sh(["root-config", "--version"]), (6, 38), "source <ROOT>/bin/thisroot.sh")
    row("htslib", sh(["pkg-config", "--modversion", "htslib"]), None, "install htslib's development package")
    samtools = sh(["samtools", "--version"])
    row("samtools", (samtools.splitlines() or [""])[0], (1, 16))
    if samtools and "libcurl=yes" not in samtools.replace(" ", ""):
        say("  WARNING samtools was built without libcurl: fetch cannot read the data from EBI")
    if git("rev-parse", "--git-dir", check=False) is None:
        say(f"  MISSING {REPO} is not a git clone: run needs one to build revisions")
        ok = False
    if args.pin is not None and not shutil.which("numactl"):
        say("  MISSING numactl, which --pin needs")
        ok = False
    env = environment([])
    say(f"\nmachine: {env['cpu']}, {env['sockets']} socket(s), {env['cores']} cores, {env['threads']} threads, "
        f"{env['memory_gb']} GB, {env.get('distribution') or env['os']}")
    if env.get("governor") and env["governor"] != "performance":
        say(f"  note: the CPU frequency governor is '{env['governor']}'; 'performance' gives steadier timings")
    return ok


# ---------------------------------------------------------------- data

def files(work, name):
    stem = STEM if name == "chr20" else f"{STEM}_{SUBSETS[name]}"
    return {"sam": work / f"{stem}.sam", "bam": work / f"{stem}.bam", "cram": work / f"{stem}.cram",
            "regions": work / f"{stem}.regions.json", "ref": work / "GRCh38.chr20.fa"}


def sequence_md5(fasta):
    h = hashlib.md5()
    with open(fasta, "rb") as f:
        for line in f:
            if not line.startswith(b">"):
                h.update(line.strip().upper())
    return h.hexdigest()


def fetch(args):
    """Downloads chr20 and its reference and writes every input the benchmarks read."""
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    t = str(args.threads)
    log = work / "fetch.log"
    whole = files(work, "chr20")
    ref, bam = whole["ref"], whole["bam"]

    if not ref.exists():
        say(f"reference: {CHROM} of GRCh38 from 1000 Genomes ...")
        part = work / "GRCh38.chr20.fa.part"
        check_call(["samtools", "faidx", "-o", part, REF_URL, CHROM], cwd=work, log=log)
        part.rename(ref)
    if sequence_md5(ref) != CHROM_M5:
        sys.exit(f"{ref}: the sequence's md5 is not {CHROM_M5}, the M5 of {CHROM} in the CRAM; delete it and fetch again")
    if not Path(str(ref) + ".fai").exists():
        check_call(["samtools", "faidx", ref], log=log)

    if not bam.exists():
        say(f"records of {args.region or CHROM} from the 30x CRAM at EBI (about 360 MB for all of {CHROM}) ...")
        part = work / f"{STEM}.part.bam"
        check_call(["samtools", "view", "--no-PG", "-b", "-@", t, "-T", ref, "-X", "-o", part,
                    CRAM_URL, CRAM_URL + ".crai", args.region or CHROM], cwd=work, log=log)
        part.rename(bam)
    if not Path(str(bam) + ".bai").exists():
        check_call(["samtools", "index", "-@", t, bam], log=log)
    verify_records(bam, args.region, args.accept_different_data)

    sam = whole["sam"]
    if not sam.exists():
        say(f"{sam.name} (about 9 GB) ...")
        part = work / f"{STEM}.part.sam"
        check_call(["samtools", "view", "--no-PG", "-h", "-@", t, "-o", part, bam], log=log)
        part.rename(sam)
    write_cram(whole, t, log)
    if not whole["regions"].exists():
        whole["regions"].write_text(json.dumps(CHR20_REGIONS) + "\n")

    for name, n in SUBSETS.items():
        f = files(work, name)
        if not (f["sam"].exists() and f["regions"].exists()):
            say(f"the first {n:,} records as {f['sam'].name} ...")
            write_subset(sam, f["sam"], f["regions"], n)
        if not f["bam"].exists():
            part = work / f"{f['bam'].stem}.part.bam"
            check_call(["samtools", "view", "--no-PG", "-b", "-@", t, "-o", part, f["sam"]], log=log)
            part.rename(f["bam"])
        if not Path(str(f["bam"]) + ".bai").exists():
            check_call(["samtools", "index", f["bam"]], log=log)
        write_cram(f, t, log)
    say(f"data ready in {work}")


def verify_records(bam, region, accept):
    """Stops when the BAM does not hold the records the benchmarks were measured on,
    unless told to accept different data."""
    stamp = Path(str(bam) + ".verified")
    if stamp.exists():
        return
    if region:
        say(f"  {region} only: the records are not checked against all of {CHROM}")
        return
    say("checking the records (count and md5) ...")
    count = sum(int(line.split("\t")[2]) + int(line.split("\t")[3])
                for line in sh(["samtools", "idxstats", bam]).splitlines() if line.strip())
    digest = md5_of_output(["samtools", "view", bam])
    if count != RECORDS or digest != RECORDS_MD5:
        problem = (f"{bam} holds {count:,} records with md5 {digest}, not {RECORDS:,} with md5 {RECORDS_MD5}. "
                   f"A download cut short gives fewer records: delete the BAM and fetch again. The same count with "
                   f"another md5 means this samtools prints the records differently; please report it with "
                   f"`samtools --version`.")
        if not accept:
            sys.exit(problem + " --accept-different-data goes on with these records.")
        say("  WARNING " + problem)
    stamp.write_text(f"{count} {digest}\n")
    say(f"  {count:,} records, md5 {digest}")


def write_cram(f, threads, log):
    """The same records as CRAM 3.1 with samtools' default profile, for the yardstick."""
    if f["cram"].exists():
        return
    part = f["cram"].with_suffix(".part.cram")
    check_call(["samtools", "view", "--no-PG", "-C", "--output-fmt-option", "version=3.1", "-T", f["ref"],
                "-@", threads, "-o", part, f["bam"]], log=log)
    part.rename(f["cram"])
    check_call(["samtools", "index", f["cram"]], log=log)


def write_subset(sam, out, regions, n):
    """The header and the first n records of sam, and regions at a quarter, half and
    three quarters of them."""
    marks = {n // 4: None, n // 2: None, (3 * n) // 4: None}
    kept = 0
    part = out.with_suffix(".part.sam")
    with open(sam, "rb") as src, open(part, "wb") as dst:
        for line in src:
            if line.startswith(b"@"):
                dst.write(line)
                continue
            if kept == n:
                break
            if kept in marks:
                fields = line.split(b"\t", 4)
                marks[kept] = (fields[2].decode(), int(fields[3]))
            dst.write(line)
            kept += 1
    if kept < n:
        sys.exit(f"{sam} has only {kept:,} records, fewer than {n:,}")
    part.rename(out)
    (r1, p1), (r2, p2), (r3, p3) = (marks[k] for k in sorted(marks))
    regions.write_text(json.dumps({"large": f"{r1}:{p1}-{p3}" if r1 == r3 else f"{r2}:{p2}-{p2 + 9_999_999}",
                                   "small": [f"{r2}:{p2}-{p2 + 999}", f"{r1}:{p1}-{p1 + 499}",
                                             f"{r3}:{p3}-{p3 + 1999}"]}) + "\n")


# ---------------------------------------------------------------- builds

def resolve(ref):
    """(id, sha, title) of a revision; WORKTREE is the checkout as it is, uncommitted changes included."""
    if ref == "WORKTREE":
        head = git("rev-parse", "HEAD")
        diff = git("diff", "HEAD", "--binary", "--", *SOURCES)
        untracked = git("ls-files", "--others", "--exclude-standard", "--", *SOURCES)
        h = hashlib.md5((head + diff).encode())
        for name in untracked.splitlines():
            h.update(name.encode() + (REPO / name).read_bytes())
        return f"wt-{h.hexdigest()[:8]}", head, f"working tree on {head[:7]}"
    sha = git("rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}", check=False)
    if not sha:
        sys.exit(f"{ref}: not a commit in {REPO}; fetch it first (git fetch upstream ...)")
    return sha[:12], sha, git("log", "-1", "--format=%s", sha)


def build(ref, args):
    """Builds a revision once into --work/builds; returns its description and tools directory."""
    bid, sha, title = resolve(ref)
    root = Path(args.work) / "builds" / bid
    tools = root / "build" / "tools"
    info = {"ref": ref, "id": bid, "sha": sha, "title": title}
    if (root / "DONE").exists():
        info.update(json.loads((root / "DONE").read_text()))
        return info, tools

    say(f"building {ref} ({bid}) ...")
    shutil.rmtree(root, ignore_errors=True)
    src = root / "src"
    src.mkdir(parents=True)
    if ref == "WORKTREE":
        names = git("ls-files", "--cached", "--others", "--exclude-standard", "--", *SOURCES).splitlines()
        for name in names:
            if (REPO / name).is_file():
                (src / name).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(REPO / name, src / name)
    else:
        archive = subprocess.run(["git", "-C", str(REPO), "archive", "--format=tar", sha],
                                 stdout=subprocess.PIPE, check=True).stdout
        subprocess.run(["tar", "xf", "-", "-C", str(src)], input=archive, check=True)

    log = root / "build.log"
    out = check_call(["cmake", "-S", src, "-B", root / "build", "-DCMAKE_BUILD_TYPE=Release",
                      "-DRAMTOOLS_BUILD_TESTS=OFF", "-DRAMTOOLS_BUILD_BENCHMARKS=OFF"], log=log)
    check_call(["cmake", "--build", root / "build", "-j", str(os.cpu_count() or 2)], log=log)
    m = re.search(r"htscodecs: (.+)", out)
    cache = (root / "build" / "CMakeCache.txt").read_text()
    cxx = re.search(r"^CMAKE_CXX_COMPILER:\w+=(.+)$", cache, re.M)
    meta = {"htscodecs": m.group(1).strip() if m else None,
            "compiler": (sh([cxx.group(1), "--version"]).splitlines() or [""])[0] if cxx else None,
            "tools": [t for t in TOOLS if (tools / t).exists()]}
    (root / "DONE").write_text(json.dumps(meta) + "\n")
    info.update(meta)
    return info, tools


# ---------------------------------------------------------------- measuring

# On macOS ru_maxrss also counts memory malloc has freed and the kernel has not
# taken back, so the peak there is the peak physical footprint, what Activity
# Monitor shows, read from the exited child before wait4 reaps it.
if sys.platform == "darwin":
    import ctypes

    _libc = ctypes.CDLL(None, use_errno=True)
    _libproc = ctypes.CDLL("/usr/lib/libproc.dylib")
    _P_PID, _WEXITED, _WNOWAIT, _RUSAGE_INFO_V4 = 1, 0x04, 0x20, 4
    _LIFETIME_MAX_PHYS_FOOTPRINT = 240  # byte offset in struct rusage_info_v4

    def peak_footprint(pid):
        siginfo = ctypes.create_string_buffer(256)
        if _libc.waitid(_P_PID, pid, siginfo, _WEXITED | _WNOWAIT) != 0:
            return None
        info = (ctypes.c_uint64 * 64)()
        if _libproc.proc_pid_rusage(pid, _RUSAGE_INFO_V4, info) != 0:
            return None
        return info[_LIFETIME_MAX_PHYS_FOOTPRINT // 8]
else:
    def peak_footprint(pid):
        return None


def exit_code(status):
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    return -os.WTERMSIG(status) if os.WIFSIGNALED(status) else 1


def timed(cmd, pin):
    """One run: wall and CPU seconds, peak memory in bytes, the load average before it."""
    cmd = [*pin, *[str(c) for c in cmd]]
    load = os.getloadavg()[0]
    with tempfile.TemporaryFile() as err:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err)
        footprint = peak_footprint(proc.pid)
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
        proc.returncode = exit_code(status)
        if proc.returncode != 0:
            err.seek(0)
            sys.exit(f"failed ({proc.returncode}): {' '.join(cmd)}\n{err.read().decode(errors='replace')[-2000:]}")
    peak = footprint or usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    return {"wall": round(wall, 4), "cpu": round(usage.ru_utime + usage.ru_stime, 4),
            "peak": int(peak), "load": round(load, 2)}


def series(name, repeat, sides, pin, warm, before=None):
    """repeat timed runs of every side in alternating order (A B, B A, ...), after one
    untimed run of each when warm; sides maps a side to its command line."""
    if repeat == 0 or not sides:
        return {}
    say(f"  {name}: {len(sides)} x {repeat} ...")
    if warm:
        for cmd in sides.values():
            if before:
                before()
            timed(cmd, pin)
    runs = {side: [] for side in sides}
    order = list(sides)
    for r in range(repeat):
        for side in (order if r % 2 == 0 else order[::-1]):
            if before:
                before()
            runs[side].append(timed(sides[side], pin))
    return runs


def storage(ram):
    """Bytes on storage per field of the RAM RNTuple, through ROOT."""
    if not shutil.which("root"):
        return {}
    macro = (f'auto r = ROOT::RNTupleReader::Open("RAM", "{ram}"); '
             f'r->PrintInfo(ROOT::ENTupleInfo::kStorageDetails);')
    out = subprocess.run(["root", "-l", "-b", "-q", "-e", macro], capture_output=True, text=True).stdout
    fields, column = {}, None
    for line in out.splitlines():
        m = re.match(r"\s+(\S+) \[#\d+\]\s+--", line)
        if m:
            column = re.sub(r"(\._\d+)+$", "", m.group(1))  # record.seq and record.seq._0 are one field
            continue
        m = re.match(r"\s+Size on storage:\s+(\d+) B", line)
        if m and column:
            fields[column] = fields.get(column, 0) + int(m.group(1))
            column = None
    return fields


def run(args):
    work = Path(args.work)
    f = files(work, args.input)
    if not all(f[k].exists() for k in ("sam", "bam", "cram", "regions", "ref")):
        say("the data is not all there; fetching it first")
        fetch(args)
    if not check(args):
        sys.exit("install what is missing first")
    pin = ["numactl", f"--cpunodebind={args.pin}", f"--membind={args.pin}"] if args.pin is not None else []
    if pin and sys.platform == "darwin":
        sys.exit("--pin needs numactl, which macOS does not have")

    base, base_tools = build(args.base, args)
    head, head_tools = build(args.head, args)
    if base["id"] == head["id"]:
        sys.exit("--base and --head are the same build")
    builds = {"base": (base, base_tools), "head": (head, head_tools)}
    regions = json.loads(f["regions"].read_text())
    sam, bam, cram, ref = f["sam"], f["bam"], f["cram"], f["ref"]
    reps = dict(REPEATS[args.input])
    for key in reps:
        if args.repeat is not None and reps[key]:
            reps[key] = args.repeat
    if args.repeat_convert is not None:
        for key in ("convert_sam", "convert_threads", "convert_bam"):
            reps[key] = args.repeat_convert

    out_dir = work / "runs" / args.input
    out_dir.mkdir(parents=True, exist_ok=True)
    out = {side: {kind: out_dir / f"{info['id']}.{kind}.ram" for kind in ("sam", "threads", "bam")}
           for side, (info, _) in builds.items()}

    def has(side, tool):
        return tool in builds[side][0]["tools"]

    def warm_input(path):
        return lambda: md5_of_file(path)  # reads the input once, cached afterwards

    steps = {}
    say(f"measuring {args.input}: base {base['ref']} ({base['id']}), head {head['ref']} ({head['id']})")
    md5_of_file(sam)  # the input is in the page cache before the first conversion
    for kind, step, src, tool, extra in (("sam", "convert_sam", sam, "samtoramntuple", []),
                                         ("threads", "convert_threads", sam, "samtoramntuple",
                                          ["-threads", str(args.threads)]),
                                         ("bam", "convert_bam", bam, "bamtoramntuple", [])):
        sides = {s: [t / tool, src, out[s][kind], *extra] for s, (_, t) in builds.items() if has(s, tool)}
        steps[step] = series(f"{step} {' '.join(extra)}".strip(), reps[step], sides, pin, warm=False)

    # Reads use the file converted from SAM on one thread, or with threads when that was left out.
    readable = {s: next((out[s][k] for k in ("sam", "threads") if out[s][k].exists()), None) for s in builds}
    if not all(readable.values()):
        sys.exit("no file converted from SAM to read; give --repeat-convert")
    dump = {s: [t / "ramdump", "-h", readable[s]] for s, (_, t) in builds.items()}
    steps["dump_h"] = series("ramdump -h", reps["read"], dump, pin, warm=True)
    steps["dump_large"] = series(f"ramdump {regions['large']}", reps["read"],
                                 {s: [t / "ramdump", readable[s], regions["large"]] for s, (_, t) in builds.items()},
                                 pin, warm=True)
    for i, region in enumerate(regions["small"]):
        steps[f"dump_small_{i + 1}"] = series(f"ramdump {region}", reps["fast"],
                                              {s: [t / "ramdump", readable[s], region]
                                               for s, (_, t) in builds.items()}, pin, warm=True)
    steps["count_large"] = series(f"ramdump -c {regions['large']}", reps["fast"],
                                  {s: [t / "ramdump", "-c", readable[s], regions["large"]]
                                   for s, (_, t) in builds.items()}, pin, warm=True)
    for stat in ("flagstat", "idxstats"):
        steps[stat] = series(f"ramstats {stat}", reps["fast"],
                             {s: [t / "ramstats", stat, readable[s]] for s, (_, t) in builds.items()
                              if has(s, "ramstats")}, pin, warm=True)

    yard = {}
    if not args.no_samtools:
        say("  samtools on the same records:")
        yard_cmds = {
            "dump_h": {"BAM": ["samtools", "view", "-h", bam], "CRAM 3.1": ["samtools", "view", "-h", "-T", ref, cram]},
            "dump_large": {"BAM": ["samtools", "view", bam, regions["large"]],
                           "CRAM 3.1": ["samtools", "view", "-T", ref, cram, regions["large"]]},
            "count_large": {"BAM": ["samtools", "view", "-c", bam, regions["large"]],
                            "CRAM 3.1": ["samtools", "view", "-c", "-T", ref, cram, regions["large"]]},
            "flagstat": {"BAM": ["samtools", "flagstat", bam], "CRAM 3.1": ["samtools", "flagstat", cram]},
            "idxstats": {"BAM": ["samtools", "idxstats", bam], "CRAM 3.1": ["samtools", "idxstats", cram]},
        }
        for step, sides in yard_cmds.items():
            repeat = reps["read"] if step.startswith("dump") else reps["fast"]
            yard[step] = series(f"samtools {step}", repeat, sides, pin, warm=True)

    say("  checking the results ...")
    correctness = check_results(builds, out, readable, f, regions)
    sizes = {s: {k: p.stat().st_size for k, p in out[s].items() if p.exists()} for s in builds}
    fields = {s: storage(readable[s]) for s in builds}
    doc = {"protocol": PROTOCOL, "measured_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "env": environment(pin), "input": {"name": args.input, "sam_md5": md5_of_file(sam), "regions": regions,
                                              "bam_size": bam.stat().st_size, "cram_size": cram.stat().st_size},
           "threads": args.threads, "repeats": reps, "builds": {s: info for s, (info, _) in builds.items()},
           "steps": steps, "samtools": yard, "sizes": sizes, "fields": fields, "correctness": correctness}
    results = Path(args.results or work / "results")
    results.mkdir(parents=True, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    path = results / f"{stamp}_{doc['env']['machine_id']}_{args.input}_{base['id']}_{head['id']}.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")
    text = render(doc)
    path.with_suffix(".md").write_text(text)
    print("\n" + text)
    say(f"saved {path} and {path.with_suffix('.md')}")
    if not args.keep:
        for side in out.values():
            for p in side.values():
                p.unlink(missing_ok=True)
    if not all(all(v is not False for v in c.values()) for c in correctness.values()):
        sys.exit(1)


def check_results(builds, out, readable, f, regions):
    """Every output that samtools or the input defines, compared; None when a build cannot say."""
    sam_md5 = md5_of_file(f["sam"])
    expect = {"bam_records": md5_of_output(["samtools", "view", f["bam"]]),
              "flagstat": sh(["samtools", "flagstat", f["bam"]]),
              "idxstats": sh(["samtools", "idxstats", f["bam"]]),
              "count": sh(["samtools", "view", "-c", f["bam"], regions["large"]]),
              "small": [md5_of_output(["samtools", "view", f["bam"], r]) for r in regions["small"]]}
    result = {}
    for side, (info, t) in builds.items():
        c = {"sam_round_trip": md5_of_output([t / "ramdump", "-h", readable[side]]) == sam_md5}
        if out[side]["threads"].exists() and readable[side] != out[side]["threads"]:
            c["threads_round_trip"] = md5_of_output([t / "ramdump", "-h", out[side]["threads"]]) == sam_md5
        if out[side]["bam"].exists():
            c["bam_records"] = md5_of_output([t / "ramdump", out[side]["bam"]]) == expect["bam_records"]
        c["count"] = sh([t / "ramdump", "-c", readable[side], regions["large"]]) == expect["count"]
        c["small_regions"] = all(md5_of_output([t / "ramdump", readable[side], r]) == m
                                 for r, m in zip(regions["small"], expect["small"]))
        if "ramstats" in info["tools"]:
            c["flagstat"] = sh([t / "ramstats", "flagstat", readable[side]]) == expect["flagstat"]
            c["idxstats"] = sh([t / "ramstats", "idxstats", readable[side]]) == expect["idxstats"]
        result[side] = c
    return result


# ---------------------------------------------------------------- reporting

LABELS = {"convert_sam": "convert SAM, 1 thread", "convert_threads": "convert SAM, -threads {threads}",
          "convert_bam": "convert BAM", "dump_h": "ramdump -h, whole file", "dump_large": "dump {large}",
          "dump_small_1": "dump {small1}", "dump_small_2": "dump {small2}", "dump_small_3": "dump {small3}",
          "count_large": "count {large}", "flagstat": "flagstat", "idxstats": "idxstats"}
CHECKS = {"sam_round_trip": "SAM round trip: ramdump -h equals the input",
          "threads_round_trip": "the same with -threads",
          "bam_records": "BAM round trip: ramdump equals samtools view",
          "count": "count of the large region equals samtools view -c",
          "small_regions": "small regions equal samtools view",
          "flagstat": "flagstat equals samtools flagstat", "idxstats": "idxstats equals samtools idxstats"}


def mb(n):
    return f"{n / 1e6:,.1f} MB"


def pct(old, new):
    return f"{(new - old) / old:+.1%}" if old else ""


def stat(runs, key):
    values = [r[key] for r in runs]
    return statistics.median(values), min(values), max(values)


def cell(runs, key, unit):
    med, lo, hi = stat(runs, key)
    if unit == "s":
        return f"{med:.2f} s" + (f" ({lo:.2f}-{hi:.2f})" if len(runs) > 1 else "")
    return f"{med / 1e6:.0f} MB"


def render(doc):
    env, inp, b = doc["env"], doc["input"], doc["builds"]
    regions = inp["regions"]
    names = {"threads": doc["threads"], "large": regions["large"],
             **{f"small{i + 1}": r for i, r in enumerate(regions["small"])}}
    lines = [f"# RAMTools benchmark: {inp['name']}, {doc['measured_at'][:10]}", "",
             f"- Machine: {env['cpu']}, {env['sockets']} socket(s), {env['cores']} cores, {env['memory_gb']} GB, "
             f"{env.get('distribution') or env['os']}" + (f"; pinned with `{env['pinned']}`" if env["pinned"] else ""),
             f"- Software: ROOT {env['root']}, {env['samtools']}, {b['head'].get('compiler')}, "
             f"htscodecs {b['head'].get('htscodecs')}",
             f"- base: `{b['base']['ref']}` {b['base']['sha'][:7]} {b['base']['title']}",
             f"- head: `{b['head']['ref']}` {b['head']['sha'][:7]} {b['head']['title']}",
             f"- Input: md5 of the SAM {inp['sam_md5']}", "",
             "## Correctness", "", "| check | base | head |", "|---|---|---|"]
    corr = doc["correctness"]
    for key, label in CHECKS.items():
        if key in corr["base"] or key in corr["head"]:
            vals = ["yes" if corr[s].get(key) else ("NO" if corr[s].get(key) is False else "n/a")
                    for s in ("base", "head")]
            lines.append(f"| {label} | {vals[0]} | {vals[1]} |")

    lines += ["", "## Size", "", "| file | base | head | change |", "|---|---:|---:|---:|"]
    for kind, label in (("sam", "converted from SAM"), ("threads", "converted from SAM with -threads"),
                        ("bam", "converted from BAM")):
        old, new = doc["sizes"]["base"].get(kind), doc["sizes"]["head"].get(kind)
        if old and new:
            lines.append(f"| {label} | {mb(old)} | {mb(new)} | {pct(old, new)} |")
    lines += [f"| samtools BAM | {mb(inp['bam_size'])} | | |", f"| samtools CRAM 3.1 | {mb(inp['cram_size'])} | | |"]
    fb, fh = doc["fields"]["base"], doc["fields"]["head"]
    if fb or fh:
        lines += ["", "By field, on storage:", "", "| field | base | head | change |", "|---|---:|---:|---:|"]
        for n in sorted(set(fb) | set(fh), key=lambda n: -max(fb.get(n, 0), fh.get(n, 0))):
            old, new = fb.get(n, 0), fh.get(n, 0)
            if max(old, new) < 100_000:
                continue
            # A field that held almost nothing on one side has no meaningful percentage.
            change = "new" if old < 10_000 else ("gone" if new < 10_000 else pct(old, new))
            lines.append(f"| {n} | {mb(old)} | {mb(new)} | {change} |")

    loads = [r["load"] for step in doc["steps"].values() for runs in step.values() for r in runs]
    lines += ["", "## Time and memory", "",
              f"Medians of the runs in alternating order, with the fastest and slowest in brackets; CPU is user plus "
              f"system time of the process and its threads; memory is the {env['memory_metric']}. Runs per step: "
              + ", ".join(f"{k} {v}" for k, v in doc["repeats"].items()) + "."
              + (f" Load average before the runs: {min(loads):.1f} to {max(loads):.1f}." if loads else ""), "",
              "| step | wall base | wall head | change | CPU base | CPU head | change | memory base | memory head |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for step, runs in doc["steps"].items():
        if "base" not in runs or "head" not in runs:
            continue
        label = LABELS.get(step, step).format(**names)
        wb, wh = stat(runs["base"], "wall")[0], stat(runs["head"], "wall")[0]
        cb, ch = stat(runs["base"], "cpu")[0], stat(runs["head"], "cpu")[0]
        lines.append(f"| {label} | {cell(runs['base'], 'wall', 's')} | {cell(runs['head'], 'wall', 's')} | "
                     f"{pct(wb, wh)} | {cell(runs['base'], 'cpu', 's')} | {cell(runs['head'], 'cpu', 's')} | "
                     f"{pct(cb, ch)} | {cell(runs['base'], 'peak', 'MB')} | {cell(runs['head'], 'peak', 'MB')} |")
    if doc["samtools"]:
        lines += ["", "## samtools on the same records", "",
                  "| step | wall BAM | CPU BAM | memory BAM | wall CRAM 3.1 | CPU CRAM 3.1 | memory CRAM 3.1 |",
                  "|---|---:|---:|---:|---:|---:|---:|"]
        for step, runs in doc["samtools"].items():
            label = LABELS.get(step, step).format(**names).replace("ramdump -h", "samtools view -h")
            row = [f"| {label} "]
            for side in ("BAM", "CRAM 3.1"):
                r = runs.get(side)
                row.append(f"| {cell(r, 'wall', 's')} | {cell(r, 'cpu', 's')} | {cell(r, 'peak', 'MB')} " if r
                           else "| | | ")
            lines.append("".join(row) + "|")
    return "\n".join(lines) + "\n"


def report(args):
    print(render(json.loads(Path(args.result).read_text())), end="")


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    default_work = os.environ.get("RAMTOOLS_BENCH_WORK", str(REPO.parent / "ramtools-bench-work"))

    def common(p):
        p.add_argument("--work", default=default_work,
                       help="data, builds and results (default %(default)s, or $RAMTOOLS_BENCH_WORK)")
        p.add_argument("--threads", type=int, default=4, help="threads for -threads and samtools (default 4)")
        p.add_argument("--pin", type=int, metavar="NODE",
                       help="Linux: run every measured command on one NUMA node with numactl")
        p.add_argument("--accept-different-data", action="store_true",
                       help="go on when the downloaded records differ from the ones the numbers were measured on")

    p = sub.add_parser("check", help="report what is installed")
    common(p)
    p = sub.add_parser("fetch", help="download and prepare the data")
    common(p)
    p.add_argument("--region", help=argparse.SUPPRESS)  # a slice of chr20, to try the script out quickly
    p = sub.add_parser("run", help="build two revisions and measure them")
    common(p)
    p.add_argument("--base", default="develop", help="revision to compare against (default develop)")
    p.add_argument("--head", default="WORKTREE",
                   help="revision to measure; WORKTREE is the checkout with its changes (default)")
    p.add_argument("--input", choices=INPUTS, default="subset",
                   help="small: 50,000 records; subset: 2 million (default); chr20: all 17.7 million")
    p.add_argument("--repeat", type=int, help="timed runs of every step that runs, instead of the input's defaults")
    p.add_argument("--repeat-convert", type=int, help="timed runs of each conversion, 0 to leave them out")
    p.add_argument("--no-samtools", action="store_true", help="leave out the samtools yardstick")
    p.add_argument("--results", help="where results go (default WORK/results)")
    p.add_argument("--keep", action="store_true", help="keep the converted files in WORK/runs")
    p.add_argument("--region", help=argparse.SUPPRESS)
    p = sub.add_parser("report", help="print a saved result again")
    p.add_argument("result", help="a .json file that run saved")
    args = ap.parse_args()

    if args.command == "check":
        sys.exit(0 if check(args) else 1)
    if args.command == "fetch":
        fetch(args)
    elif args.command == "run":
        run(args)
    else:
        report(args)


if __name__ == "__main__":
    main()
