#!/usr/bin/env python3
"""Are within-cycle repeats the SAME SPECTRUM, or different spectra hitting the same peptide?

64.9% of our PSMs are within-cycle repeats (the reference implementation 20.3%), and 2.85 PSMs come from each
(peptide, cycle). Sharing a precursor is not the same as being redundant: two spectra can share
a precursor assignment and still carry different fragment content.

That distinction decides what may be done about it:

  HIGH cosine  -> near-duplicate spectra. Merging is safe and loses nothing.
                  `consolidate:delta_rt` is the right tool.
  LOW cosine   -> distinct fragment content that happens to match the same peptide. Merging
                  DESTROYS information, and both consolidation and any NMF/NNLS cleanup would
                  be actively harmful. The right response is to stop generating them.

This is the rank-1 test from the backlog, done directly: cosine is the rank-1 reconstruction
quality for a pair.

Controls, because a similarity number without a null means nothing here -- the MS1 funnel was
retracted for exactly that:
  * RANDOM pairs from DIFFERENT cycles, same peptide  -> how similar are two spectra of the
    same peptide that are NOT within-cycle repeats?
  * RANDOM pairs of unrelated spectra                 -> the floor for "these share nothing".
A within-cycle cosine only means something relative to those two.
"""
import argparse, collections, csv, math, random, sys
from pyteomics import mzml

CYCLE_S = 1.385
BIN = 0.02          # m/z bin for the cosine; 0.02 Da is ~20 ppm at m/z 1000


def binned(mzs, ints):
    v = collections.defaultdict(float)
    for m, i in zip(mzs, ints):
        v[int(m / BIN)] += float(i)
    n = math.sqrt(sum(x * x for x in v.values()))
    return {k: x / n for k, x in v.items()} if n > 0 else {}


def cosine(a, b):
    if not a or not b:
        return 0.0
    small, large = (a, b) if len(a) < len(b) else (b, a)
    return sum(x * large.get(k, 0.0) for k, x in small.items())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    ap.add_argument("mzml")
    ap.add_argument("--label", default="run")
    ap.add_argument("--max-groups", type=int, default=4000)
    a = ap.parse_args()
    random.seed(11)

    # ---- PSMs at 1% FDR, grouped by (peptide, cycle) -------------------------------------------
    rows = []
    with open(a.tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for x in rd:
            try:
                rows.append({"pep": x["peptide"], "sid": x["scannr"],
                             "cyc": int(float(x["rt"]) * 60 / CYCLE_S),
                             "score": float(x["sage_discriminant_score"]),
                             "decoy": all(p.startswith("rev_")
                                          for p in x["proteins"].split(";") if p)})
            except (TypeError, ValueError, KeyError):
                continue
    ranked = sorted(rows, key=lambda r: -r["score"])
    nt = nd = 0
    keep, acc = [], []
    for r in ranked:
        if r["decoy"]:
            nd += 1
        else:
            nt += 1
            acc.append(r)
        if nt and nd / nt <= 0.01:
            keep = list(acc)

    g = collections.defaultdict(list)
    for r in keep:
        g[(r["pep"], r["cyc"])].append(r["sid"])
    dup = {k: v for k, v in g.items() if len(v) > 1}
    bypep = collections.defaultdict(list)     # same peptide, DIFFERENT cycles -> control
    for (pep, cyc), sids in g.items():
        bypep[pep].append((cyc, sids[0]))

    sel = dict(random.sample(sorted(dup.items()), min(a.max_groups, len(dup))))
    print("[%s] %d duplicated (peptide,cycle) groups; sampling %d"
          % (a.label, len(dup), len(sel)), flush=True)

    # control pairs: same peptide, different cycle
    ctrl_pairs = []
    for pep, entries in bypep.items():
        if len(entries) > 1:
            (c1, s1), (c2, s2) = random.sample(entries, 2)
            if c1 != c2:
                ctrl_pairs.append((s1, s2))
    ctrl_pairs = ctrl_pairs[:len(sel)]

    need = set()
    for sids in sel.values():
        need.update(sids[:4])                  # cap: 4 members per group is plenty for a mean
    for s1, s2 in ctrl_pairs:
        need.add(s1); need.add(s2)

    # ---- single pass over the mzML -------------------------------------------------------------
    print("[%s] reading %d spectra from %s" % (a.label, len(need), a.mzml), flush=True)
    spec = {}
    for s in mzml.read(a.mzml):
        sid = s.get("id")
        if sid in need:
            spec[sid] = binned(s["m/z array"], s["intensity array"])
            if len(spec) == len(need):
                break

    def pairs_of(sids):
        out = []
        for i in range(len(sids)):
            for j in range(i + 1, len(sids)):
                if sids[i] in spec and sids[j] in spec:
                    out.append(cosine(spec[sids[i]], spec[sids[j]]))
        return out

    within = [c for sids in sel.values() for c in pairs_of(sids[:4])]
    across = [cosine(spec[x], spec[y]) for x, y in ctrl_pairs if x in spec and y in spec]
    allsp = [v for v in spec.values() if v]
    rnd = []
    for _ in range(min(3000, len(allsp) // 2)):
        x, y = random.sample(allsp, 2)
        rnd.append(cosine(x, y))

    def q(v, p):
        v = sorted(v)
        return v[int(p * (len(v) - 1))] if v else float("nan")

    print()
    print("%-34s %7s %8s %8s %8s %8s" % ("cosine similarity", "n", "median", "p25", "p75", ">0.9"))
    print("-" * 78)
    for name, v in (("WITHIN-cycle repeats (same pep)", within),
                    ("same peptide, DIFFERENT cycle", across),
                    ("random unrelated spectra", rnd)):
        if v:
            print("%-34s %7d %8.3f %8.3f %8.3f %7.1f%%"
                  % (name, len(v), q(v, .5), q(v, .25), q(v, .75),
                     100.0 * sum(1 for c in v if c > 0.9) / len(v)))
    print("-" * 78)
    print("high within-cycle cosine -> near-duplicates, merging is safe")
    print("low  within-cycle cosine -> distinct content, merging DESTROYS signal")


if __name__ == "__main__":
    sys.exit(main())
