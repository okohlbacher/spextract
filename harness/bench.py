#!/usr/bin/env python3
"""SpeXtractor benchmark harness -- the ONLY sanctioned way to benchmark this tool.

Ad-hoc shell scripts produced eleven distinct data-handling failures in a single session, all
bookkeeping rather than algorithmic. Each guard below exists because a specific failure
happened, and the failure is named in the comment. Guards ABORT rather than warn: a benchmark
that silently produces a plausible wrong number is worse than one that refuses to run.

  FAILURE 1  cross-sample comparison. Published "dataset A best 7,430 -> 66.2%"; 7,430 was dataset D and
             11,218 was dataset A's reference.
             GUARD: samples are IDs resolved from samples.yaml; a sample fingerprint must
             appear in the raw path AND the produced mzML; references are looked up per
             sample and cannot be crossed.

  FAILURE 2  mismatched baselines. "+36.7%/+37.4% confirmed" used three different baseline
             configs; the resulting dataset D "+4.1%" was pure artefact.
             GUARD: an arm declares `baseline`; the pair runs in ONE invocation, same binary,
             same node, same everything but the varied parameters, and the diff is printed.

  FAILURE 3  wrong FDR metric. spectrum_q inflated every ratio by 9-11 points vs the reference implementation.
             GUARD: peptide_q is THE metric. spectrum_q is reported only alongside fdr_loss.

  FAILURE 4  missing controls. The MS1 funnel's "97.8%" had a 91.9% chance floor and was
             retracted; co-isolation was ~half random matching.
             GUARD: every truth-set metric carries a decoy floor computed the same way, and
             uncorrected values are never emitted alone.

  FAILURE 5  no-op arms. Three arms passed flags that were already defaults (~2.5 h wasted).
             GUARD: resolved configs are diffed; an arm identical to its baseline ABORTS.

  FAILURE 6  scale mismatch. ambiguity_margin 0.1 against an INTEGER partner count could never
             fire; three gate thresholds were compared across incompatible statistics.
             GUARD: arms may declare `scale_note`; sweeps must give >=3 points and the harness
             refuses a single-point comparison of a rescaled parameter.

  FAILURE 7  stale binary. tar silently skipped overwriting a running executable; a 6 h run
             used a 2-day-old build.
             GUARD: binary sha256 + mtime recorded in every manifest and asserted identical
             across a comparison; every non-default parameter is probed before the run.

  FAILURE 8  unit mismatch. RT in seconds vs minutes matched 31 of ~300k spectra.
             GUARD: RT is read from the CV unit accessor, never the bare number.

  FAILURE 9  asserted absence. "No the reference implementation dataset D" twice; it was on ceph all along.
             GUARD: references are declared in samples.yaml and existence is CHECKED, so
             absence is a fact about the registry rather than about where someone looked.

  FAILURE 10 provenance loss. bestfree/best_25 had no record of its sample.
             GUARD: every run writes manifest.json (sample, fingerprint, binary hash, full
             resolved config, git ref, host, timestamps) beside its outputs.

  FAILURE 11 orphan process. A stale run perturbed a benchmark's memory-adaptive admission.
             GUARD: pre-flight refuses to start if another spextractor is running.

Usage:
    bench.py run   --plan plan.yaml            # execute arms + baselines
    bench.py truth --plan plan.yaml            # truth-set metrics on existing runs
    bench.py collate --plan plan.yaml          # write the vault report
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, re, shutil, subprocess, sys, time
from dataclasses import dataclass, field, asdict
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("need pyyaml: python3 -m pip install --user pyyaml")

HERE = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------- infrastructure
def die(msg: str) -> None:
    """Guards abort. A benchmark that produces a plausible wrong number is worse than one
    that refuses to run -- every retraction this project has issued came from the former."""
    print(f"\nABORT: {msg}\n", file=sys.stderr)
    sys.exit(2)


def sh(cmd: list[str] | str, capture=True, check=False, timeout=None) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=capture,
                          text=True, check=check, timeout=timeout)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()[:16]


# ----------------------------------------------------------------------------- config objects
@dataclass
class Arm:
    name: str
    params: dict                      # explicit, non-default parameters only
    baseline: str | None = None       # name of the arm this is compared against
    scale_note: str | None = None     # FAILURE 6: parameter whose scale differs from a sibling


@dataclass
class Plan:
    title: str
    samples: list[str]
    arms: list[Arm]
    binary: str
    git_ref: str | None = None
    threads: int = 120
    tcmalloc: bool = True


def load_plan(path: Path) -> tuple[Plan, dict]:
    reg = yaml.safe_load((HERE / "samples.yaml").read_text())
    raw = yaml.safe_load(Path(path).read_text())
    arms = [Arm(name=a["name"], params=a.get("params", {}) or {},
                baseline=a.get("baseline"), scale_note=a.get("scale_note"))
            for a in raw["arms"]]
    plan = Plan(title=raw["title"], samples=raw["samples"], arms=arms,
                binary=raw["binary"], git_ref=raw.get("git_ref"),
                threads=raw.get("threads", 120), tcmalloc=raw.get("tcmalloc", True))
    return plan, reg


# ----------------------------------------------------------------------------- pre-flight
def resolved_defaults(binary: Path) -> dict:
    """Every registered parameter and its default, from the tool's own -write_ini.

    FAILURE 5: three arms passed values that were already the defaults and were byte-identical
    to baseline. Knowing the defaults is what makes a no-op detectable before it costs an hour.

    Was parsed from --helphelp "(default: '...')" text; OpenMS 3.6 builds no longer print that,
    so the parser silently returned {} and the FAILURE-5 guard aborted. -write_ini is the
    authoritative source across versions: a full INI whose <ITEM> paths, minus the tool/instance
    wrapper NODEs, are exactly the CLI parameter names (assembly:rp_max, trace:ms2_split_valleys).
    """
    import tempfile, xml.etree.ElementTree as ET
    ini = Path(tempfile.mkdtemp()) / "defs.ini"
    sh([str(binary), "-write_ini", str(ini)])
    defs: dict[str, str] = {}
    try:
        root = ET.parse(ini).getroot()          # <PARAMETERS>
    except Exception:
        return defs
    inst = None                                  # NODE[tool] > NODE["1"] = the instance section
    for tool in root.findall("NODE"):
        cand = tool.find('NODE[@name="1"]')
        if cand is not None:
            inst = cand
            break
    if inst is None:
        return defs

    def walk(node, prefix):
        for ch in node:
            nm = ch.get("name")
            if ch.tag == "NODE":
                walk(ch, prefix + [nm])
            elif ch.tag in ("ITEM", "ITEMLIST"):
                defs[":".join(prefix + [nm])] = ch.get("value", "")

    walk(inst, [])
    return defs


def preflight(plan: Plan, reg: dict) -> dict:
    binary = Path(plan.binary)
    if not binary.is_file():
        die(f"binary not found: {binary}")

    # FAILURE 11: an orphaned run perturbed a benchmark's memory-adaptive admission and
    # competed for cores for 37 minutes without being noticed.
    ps = sh("ps -eo pid,etime,args").stdout
    running = [l for l in ps.splitlines()
               if "spextractor" in l and "grep" not in l and str(os.getpid()) not in l]
    if running:
        die("another spextractor is already running -- it will contend for cores and\n"
            "       perturb memory-adaptive admission. Kill it or wait:\n         "
            + "\n         ".join(x.strip()[:110] for x in running[:5]))

    # FAILURE 7: tar silently skips overwriting a running executable, so a 6 h benchmark ran
    # against a 2-day-old build. Identity is recorded, and asserted equal across the comparison.
    bmeta = {"path": str(binary), "sha256": sha256(binary),
             "mtime": time.strftime("%F %T", time.localtime(binary.stat().st_mtime))}
    print(f"[preflight] binary {bmeta['sha256']}  built {bmeta['mtime']}")

    defs = resolved_defaults(binary)
    if not defs:
        die("could not parse tool defaults -- no-op detection (FAILURE 5) would be disabled")

    # FAILURE 7 again: probe every non-default parameter BEFORE committing hours.
    for arm in plan.arms:
        for k, v in arm.params.items():
            probe = sh([str(binary), "-in", "/nonexistent.d", f"-{k}", str(v)]).stdout \
                  + sh([str(binary), "-in", "/nonexistent.d", f"-{k}", str(v)]).stderr
            if "Unknown option" in probe:
                die(f"arm '{arm.name}': parameter -{k} does not exist in this binary")

    # FAILURE 5: an arm whose resolved config equals its baseline's tests nothing.
    by_name = {a.name: a for a in plan.arms}
    for arm in plan.arms:
        if not arm.baseline:
            continue
        if arm.baseline not in by_name:
            die(f"arm '{arm.name}' references unknown baseline '{arm.baseline}'")
        base = by_name[arm.baseline]
        eff_a = {**{k: defs.get(k) for k in arm.params}, **{k: str(v) for k, v in arm.params.items()}}
        eff_b = {**{k: defs.get(k) for k in base.params}, **{k: str(v) for k, v in base.params.items()}}
        diff = {k: (eff_b.get(k, defs.get(k)), eff_a.get(k, defs.get(k)))
                for k in set(eff_a) | set(eff_b)
                if eff_a.get(k, defs.get(k)) != eff_b.get(k, defs.get(k))}
        if not diff:
            die(f"arm '{arm.name}' resolves identically to baseline '{arm.baseline}'.\n"
                f"       This is FAILURE 5: passing a value that is already the default.\n"
                f"       defaults seen: "
                + ", ".join(f"{k}={defs.get(k)}" for k in arm.params))
        print(f"[preflight] {arm.name} vs {arm.baseline}: " +
              ", ".join(f"{k} {a} -> {b}" for k, (a, b) in sorted(diff.items())))

    # FAILURE 9: absence must be a fact about the registry, not about where someone looked.
    for s in plan.samples:
        if s not in reg["samples"]:
            die(f"sample '{s}' not in samples.yaml")
        entry = reg["samples"][s]
        if not entry.get("raw"):
            die(f"sample '{s}' has no raw path in samples.yaml")
        if not Path(entry["raw"]).exists():
            die(f"sample '{s}' raw missing: {entry['raw']}")
        # FAILURE 1: the fingerprint must actually be in the path we are about to run.
        if entry["fingerprint"] not in entry["raw"]:
            die(f"sample '{s}' fingerprint '{entry['fingerprint']}' absent from its raw path")
        for rk, rv in (entry.get("references") or {}).items():
            print(f"[preflight] {s}.{rk}: {'OK' if rv and Path(rv).exists() else 'MISSING'}")
    return bmeta


# ----------------------------------------------------------------------------- execution
def verify_sample_identity(mzml: Path, fingerprint: str, sample: str) -> None:
    """FAILURE 1: bestfree/best_25 was an dataset D run reported as dataset A for a whole day. The sample
    is now read back out of the produced file rather than trusted from the invocation."""
    head = mzml.open("rb").read(400_000).decode("utf-8", "ignore")
    if fingerprint not in head:
        found = re.findall(r"S\d+-A-\d+", head)
        die(f"{mzml} does not contain fingerprint '{fingerprint}' for sample {sample}.\n"
            f"       found instead: {sorted(set(found))[:5] or 'nothing'}\n"
            f"       This is FAILURE 1 -- the run is not the sample it claims to be.")


def run_one(plan: Plan, reg: dict, arm: Arm, sample: str, outroot: Path, bmeta: dict) -> Path:
    entry = reg["samples"][sample]
    out = outroot / f"{sample}__{arm.name}"
    out.mkdir(parents=True, exist_ok=True)
    mzml = out / "pseudo.mzML"

    cmd = [plan.binary, "-in", entry["raw"], "-out", str(mzml), "-threads", str(plan.threads)]
    for k, v in arm.params.items():
        cmd += [f"-{k}", str(v)]

    env = dict(os.environ)
    if plan.tcmalloc:
        tc = "/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"
        if Path(tc).exists():
            env["LD_PRELOAD"] = tc

    t0 = time.time()
    with open(out / "run.log", "w") as log:
        proc = subprocess.run(["/usr/bin/time", "-v"] + cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    wall = time.time() - t0

    if not mzml.exists() or mzml.stat().st_size == 0:
        die(f"{sample}/{arm.name}: tool produced no output (exit {proc.returncode}); see {out/'run.log'}")
    verify_sample_identity(mzml, entry["fingerprint"], sample)

    rss = 0
    for line in (out / "run.log").read_text(errors="ignore").splitlines():
        if "Maximum resident set size" in line:
            rss = int(re.findall(r"(\d+)\s*$", line)[0]) // 1024

    # FAILURE 10: provenance travels WITH the artefact, so a number can never be orphaned
    # from the sample and config that produced it.
    manifest = {
        "sample": sample, "fingerprint": entry["fingerprint"], "raw": entry["raw"],
        "arm": arm.name, "baseline": arm.baseline, "params": arm.params,
        "binary": bmeta, "git_ref": plan.git_ref, "threads": plan.threads,
        "host": platform.node(), "wall_s": round(wall), "peak_rss_mb": rss,
        "started": time.strftime("%F %T", time.localtime(t0)),
        "harness_version": HARNESS_VERSION,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[run] {sample}/{arm.name}: {round(wall)}s {rss}MB -> {out}")
    return out


def run_search(reg: dict, out: Path) -> None:
    s = reg["search"]
    sh([s["sage_binary"], s["sage_config"], "-o", str(out),
        "--disable-telemetry-i-dont-want-to-improve-sage", str(out / "pseudo.mzML")],
       capture=True)


def score(out: Path) -> dict:
    """FAILURE 3: peptide_q is THE metric. spectrum_q is kept only to expose fdr_loss, which
    is itself diagnostic -- ours runs -10% to -23%, the reference implementation -5% to -7%, and that gap IS the
    redundancy signal."""
    import csv
    f = out / "results.sage.tsv"
    if not f.exists():
        return {"error": "no sage output"}
    sp, pp, npsm = set(), set(), 0
    for x in csv.DictReader(f.open(), delimiter="\t"):
        if float(x.get("spectrum_q", 1)) <= 0.01:
            sp.add(x["peptide"]); npsm += 1
        if float(x.get("peptide_q", 1)) <= 0.01:
            pp.add(x["peptide"])
    return {"psms": npsm, "spectrum_q_peptides": len(sp), "peptide_q_peptides": len(pp),
            "fdr_loss_pct": round(100.0 * (len(pp) - len(sp)) / max(len(sp), 1), 1)}


HARNESS_VERSION = "1.0"


def cmd_run(args) -> None:
    plan, reg = load_plan(Path(args.plan))
    outroot = Path(args.out or f"/path/to/scratch/harness/{Path(args.plan).stem}")
    outroot.mkdir(parents=True, exist_ok=True)
    bmeta = preflight(plan, reg)

    # FAILURE 2: baseline and treatment run in ONE invocation on ONE node with ONE binary,
    # so a delta can never be taken across mismatched configurations again.
    for sample in plan.samples:
        for arm in plan.arms:
            out = run_one(plan, reg, arm, sample, outroot, bmeta)
            run_search(reg, out)
            (out / "metrics.json").write_text(json.dumps(score(out), indent=2))
    print(f"\n[done] {outroot}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "truth", "collate"):
        p = sub.add_parser(name)
        p.add_argument("--plan", required=True)
        p.add_argument("--out")
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
    else:
        from collate import cmd_truth, cmd_collate
        (cmd_truth if args.cmd == "truth" else cmd_collate)(args)


if __name__ == "__main__":
    main()
