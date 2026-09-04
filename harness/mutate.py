#!/usr/bin/env python3
"""Mutation testing: the ONLY acceptance criterion for a guard.

Both harness generations shipped guards that were documented, "tested", and inert.
v1 had six such guards. v2 replaced them, added a 64/88-case selftest, and reproduced
the same failure -- including on G-B05, the guard for the cross-sample comparison that
was actually published as a wrong result. Disabling that guard's comparison by hand
left the suite reporting "64 passed, 0 FAILED".

A passing selftest is not evidence. The question a guard must answer is:

    IF I DELETE THIS GUARD, DOES THE SUITE NOTICE?

This tool answers it mechanically for every guard site. For each `die("G-XXX", ...)`
call it substitutes a no-op, runs the module's selftest, and records the outcome:

    KILLED    selftest failed  -> the guard is genuinely exercised
    SURVIVED  selftest passed  -> THEATRE: the guard can be deleted silently

A SURVIVED site is not necessarily wrong code -- it may be a real guard with no test.
It is a claim without evidence, which is exactly what produced both failures. Every
SURVIVED site must either gain a test that kills it or be struck from the guard table
and moved to the README's UNCOVERED section.

Usage:
    mutate.py bench2.py            # sweep every guard site
    mutate.py collate2.py --guard G-C14
    mutate.py bench2.py --json out.json
"""
import argparse, ast, json, shutil, subprocess, sys, tempfile
from pathlib import Path


def str_literal(node):
    """Value of a string-literal AST node, or None.

    Handles both ast.Constant and the legacy ast.Str: this interpreter still emits Str,
    and matching only Constant silently found ZERO guard sites -- a measuring instrument
    that reports perfect coverage because it measured nothing. Exactly the failure class
    this tool exists to catch, so it is guarded rather than assumed.
    """
    v = getattr(node, "value", None)
    if isinstance(v, str):
        return v
    s = getattr(node, "s", None)
    return s if isinstance(s, str) else None


def call_end_line(lines, start_lineno):
    """Last 1-based line of a call starting at start_lineno, by paren balance.

    ast.Call.end_lineno does not exist before Python 3.8 and this interpreter is 3.7.
    Balancing parentheses over the raw text is version-independent; string contents are
    skipped so a ')' inside a message cannot end the call early.
    """
    depth, i, started = 0, start_lineno - 1, False
    while i < len(lines):
        ln, j, instr, q = lines[i], 0, False, ""
        while j < len(ln):
            c = ln[j]
            if instr:
                if c == "\\":
                    j += 2; continue
                if c == q:
                    instr = False
            elif c in "\"'":
                instr, q = True, c
            elif c == "(":
                depth += 1; started = True
            elif c == ")":
                depth -= 1
                if started and depth == 0:
                    return i + 1
            j += 1
        i += 1
    return start_lineno


def die_sites(src, tree):
    """Every call to die(...) whose first argument is a "G-*" literal, with line span."""
    src_lines = src.splitlines(keepends=True)
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        if not (isinstance(fn, ast.Name) and fn.id == "die"):
            continue
        if not node.args:
            continue
        gid = str_literal(node.args[0])
        if not (gid and gid.startswith("G-")):
            continue

        end = call_end_line(src_lines, node.lineno)
        out.append({
            "guard": gid,
            "lineno": node.lineno,
            "end_lineno": end,
            "col": node.col_offset,
            # Raw source of the call. Criticals pin to a fragment of THIS, never to a line
            # number: adding a test above a guard shifts every line below it, and a
            # line-pinned critical then matches nothing and passes silently. That happened.
            "text": "".join(src_lines[node.lineno - 1:end]),
        })
    # deterministic order; dedupe identical spans (decorated/nested walks)
    seen, uniq = set(), []
    for s in sorted(out, key=lambda d: (d["lineno"], d["guard"])):
        k = (s["lineno"], s["end_lineno"])
        if k not in seen:
            seen.add(k)
            uniq.append(s)
    return uniq


def apply_mutant(lines, site):
    """Replace the die(...) call with a no-op of identical indentation.

    Replacing the CALL rather than the enclosing `if` is deliberate: it works uniformly
    for guards inside conditionals, loops, try blocks and comprehensions, and it isolates
    the abort itself -- which is the behaviour under test -- from the detection logic.
    """
    out = list(lines)
    i0, i1 = site["lineno"] - 1, site["end_lineno"] - 1
    indent = " " * site["col"]
    out[i0] = indent + "pass  # MUTANT: %s disabled\n" % site["guard"]
    for i in range(i0 + 1, i1 + 1):
        out[i] = ""
    return out


def run_selftest(path, timeout):
    try:
        p = subprocess.run([sys.executable, str(path), "selftest"],
                           capture_output=True, text=True, timeout=timeout)
        return p.returncode, (p.stdout + p.stderr)
    except subprocess.TimeoutExpired:
        # A hang is not a pass. Treat it as killed-by-timeout and flag it separately:
        # silently counting a timeout as SURVIVED would understate coverage, and counting
        # it as KILLED would overstate it.
        return None, "TIMEOUT"


def main():
    ap = argparse.ArgumentParser(description="Mutation-test harness guards.")
    ap.add_argument("module")
    ap.add_argument("--guard", help="restrict to one guard id")
    ap.add_argument("--json", help="write full results here")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--critical", metavar="FILE", help=(
        "newline-delimited 'GUARD:LINE' or 'GUARD' entries that MUST be killed. Exit 3 if any "
        "survives. This is the gate: guards covering failures that actually recurred are held "
        "to 100%% mutation coverage, while the rest are merely REPORTED. Claiming protection "
        "everywhere and verifying it nowhere is what produced two generations of theatre."))
    args = ap.parse_args()

    mod = Path(args.module).resolve()
    src = mod.read_text()
    tree = ast.parse(src)
    lines = src.splitlines(keepends=True)
    sites = die_sites(src, tree)
    if args.guard:
        sites = [s for s in sites if s["guard"] == args.guard]
    if not sites:
        raw = src.count('die("G-')
        print("no guard sites found by AST walk, but the source contains %d literal "
              "'die(\"G-' occurrences.\nThe extractor is broken -- refusing to report "
              "coverage from a measurement that measured nothing." % raw, file=sys.stderr)
        return 2

    # Baseline must be green, or every mutant "fails" for an unrelated reason and the
    # whole sweep reads as perfect coverage.
    rc, out = run_selftest(mod, args.timeout)
    if rc != 0:
        print("ABORT: baseline selftest is not green (rc=%s). Mutation results would be "
              "meaningless -- every mutant would appear KILLED.\n%s" % (rc, out[-2000:]),
              file=sys.stderr)
        return 2
    print("baseline selftest: PASS (%d guard sites to mutate)\n" % len(sites))

    backup = Path(tempfile.mkdtemp()) / mod.name
    shutil.copy2(mod, backup)
    results = []
    try:
        for n, site in enumerate(sites, 1):
            mod.write_text("".join(apply_mutant(lines, site)))
            rc, out = run_selftest(mod, args.timeout)
            if rc is None:
                verdict = "TIMEOUT"
            elif rc != 0:
                verdict = "KILLED"
            else:
                verdict = "SURVIVED"
            results.append({**site, "verdict": verdict})
            print("[%3d/%3d] %-8s line %-5d %s" % (n, len(sites), site["guard"],
                                                   site["lineno"], verdict), flush=True)
    finally:
        shutil.copy2(backup, mod)

    # --- summary, per guard id -------------------------------------------------------
    by = {}
    for r in results:
        b = by.setdefault(r["guard"], {"KILLED": 0, "SURVIVED": 0, "TIMEOUT": 0, "lines": []})
        b[r["verdict"]] += 1
        if r["verdict"] == "SURVIVED":
            b["lines"].append(r["lineno"])

    killed = sum(1 for r in results if r["verdict"] == "KILLED")
    survived = [r for r in results if r["verdict"] == "SURVIVED"]
    print("\n" + "=" * 88)
    print("%-10s %8s %9s %8s   %s" % ("GUARD", "KILLED", "SURVIVED", "COVER", "unverified lines"))
    print("-" * 88)
    for g in sorted(by):
        b = by[g]
        tot = b["KILLED"] + b["SURVIVED"] + b["TIMEOUT"]
        cov = 100.0 * b["KILLED"] / max(tot, 1)
        mark = "" if b["SURVIVED"] == 0 else "  <-- THEATRE" if b["KILLED"] == 0 else ""
        ln = ",".join(str(x) for x in b["lines"][:6]) + ("..." if len(b["lines"]) > 6 else "")
        print("%-10s %8d %9d %7.0f%%   %s%s" % (g, b["KILLED"], b["SURVIVED"], cov, ln, mark))
    print("-" * 88)
    print("%d/%d guard sites KILLED (%.0f%%). %d SURVIVED -- unverified claims."
          % (killed, len(results), 100.0 * killed / max(len(results), 1), len(survived)))
    whole = [g for g, b in by.items() if b["KILLED"] == 0]
    if whole:
        print("\nGUARDS THAT ARE PURE THEATRE (no site killed): %s" % ", ".join(sorted(whole)))

    if args.json:
        Path(args.json).write_text(json.dumps(
            {"module": str(mod), "results": results, "by_guard": by}, indent=2))
        print("\nwrote %s" % args.json)

    # --- critical gate ---------------------------------------------------------------
    # Without a declared critical set this is only a report, and a report is what let 121
    # unverified claims ship twice. Guards covering RECURRING failures must be killed.
    if args.critical:
        want = [l.split("#")[0].strip() for l in Path(args.critical).read_text().splitlines()]
        want = [w for w in want if w]
        breaches, unmatched = [], []
        for w in want:
            mod_name, _, rest = w.partition(" ")
            guard, _, frag = rest.partition(" ")
            if not (mod_name and guard and frag):
                unmatched.append("%s -- malformed (want: MODULE GUARD MESSAGE_FRAGMENT)" % w)
                continue
            if mod_name != mod.name:
                continue                       # belongs to the other module's sweep
            hits = [r for r in results
                    if r["guard"] == guard and frag.lower() in r["text"].lower()]
            if not hits:
                # An entry that matches NOTHING must be an ERROR. Silently passing an
                # unmatchable critical is how the first version of this gate reported
                # "PASSED (6 entries)" while checking two.
                unmatched.append("%s -- matched no guard site in %s" % (w, mod.name))
                continue
            for h in hits:
                if h["verdict"] != "KILLED":
                    breaches.append("%s %s (line %d) %s" % (mod_name, guard, h["lineno"], frag))
        print()
        if unmatched:
            print("CRITICAL LIST IS STALE -- these entries match nothing and were NOT checked:")
            for u in unmatched:
                print("    %s" % u)
        if breaches:
            print("CRITICAL GATE FAILED -- these guards claim to prevent a failure that has")
            print("ALREADY HAPPENED in this project, and deleting them changes nothing:")
            for b in breaches:
                print("    %s" % b)
        if breaches or unmatched:
            return 3
        checked = sum(1 for w in want if w.split(" ")[0] == mod.name)
        print("CRITICAL GATE PASSED for %s (%d entries verified by mutation)" % (mod.name, checked))

    # Non-zero while any site survives: this is a gate, not a report.
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())
