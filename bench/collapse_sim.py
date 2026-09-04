#!/usr/bin/env python3
"""Does collapsing same-precursor-across-cycles spectra reduce the open-search testing load enough to
help? Simulated on EXISTING search results (no re-search): group rank-1 PSMs by precursor coordinate
ACROSS RT (charge + expmass ppm + ion_mobility, single-link chained over an RT gap), keep ONE per group,
then RECOMPUTE the entrapment FDR from scores (peptide_q is stale after subsetting) and count target
peptides at <=1% corrected FDR.

Selection rules:
  best  = highest sage_discriminant_score in group  -> ORACLE CEILING (uses the search result; not implementable)
  apex  = highest ms2_intensity in group            -> REALISTIC (ID-agnostic, an actual tool could do this)
  full  = no collapse (control)
"""
import sys, os, csv, argparse, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # entrapment.py is a sibling
from entrapment import peptide_ratio, _bootstrap_fdr

def load(tsv, tag="ENTRAP_"):
    rows = []
    with open(tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                if int(x["rank"]) != 1: continue
                prots = [p for p in x["proteins"].split(";") if p]
                if not prots or all(p.startswith("rev_") for p in prots): continue
                rows.append((
                    float(x["sage_discriminant_score"]), x["peptide"],
                    all(p.startswith(tag) or p.startswith("rev_" + tag) for p in prots),
                    float(x["expmass"]), int(x["charge"]), float(x["rt"]),
                    float(x.get("ion_mobility") or 0.0), float(x.get("ms2_intensity") or 0.0)))
            except (ValueError, KeyError):
                continue
    return rows

def _chain(rows, keyfn, tol_fn):
    """Single-link chain a sorted list: consecutive items within tolerance join the same group."""
    out, cur = [], [rows[0]]
    for r in rows[1:]:
        if keyfn(r) - keyfn(cur[-1]) <= tol_fn(r): cur.append(r)
        else: out.append(cur); cur = [r]
    out.append(cur)
    return out

def collapse(rows, rule, ppm=10.0, dim=0.01, rt_gap=10.0):
    """Same precursor ACROSS CYCLES: chain on charge+mass(ppm), then IM, then RT. Keep one per group."""
    if rule == "full" or not rows: return rows
    kept = []
    by_z = {}
    for r in rows: by_z.setdefault(r[4], []).append(r)          # r[4]=charge
    for grp in by_z.values():
        grp.sort(key=lambda r: r[3])                            # r[3]=expmass
        for mgrp in _chain(grp, lambda r: r[3], lambda r: r[3] * ppm * 1e-6):
            mgrp.sort(key=lambda r: r[6])                       # r[6]=ion_mobility
            for igrp in _chain(mgrp, lambda r: r[6], lambda r: dim):
                igrp.sort(key=lambda r: r[5])                   # r[5]=rt
                for rgrp in _chain(igrp, lambda r: r[5], lambda r: rt_gap):
                    kept.extend(pick(rgrp, rule))
    return kept

def pick(chain, rule):
    """Return the rows a collapsed group would still yield."""
    if rule == "apex": return [max(chain, key=lambda r: r[7])]      # ID-agnostic: highest-intensity spectrum
    if rule == "best": return [max(chain, key=lambda r: r[0])]      # oracle ceiling
    if rule.startswith("top"):                                       # merged spectrum + chimeric report_psms=N
        n = int(rule[3:]); out, seen = [], set()
        for r in sorted(chain, key=lambda r: -r[0]):
            if r[1] in seen: continue                                # distinct peptides only
            seen.add(r[1]); out.append(r)
            if len(out) >= n: break
        return out
    return [max(chain, key=lambda r: r[0])]

def targets_at_fdr(rows, ratio, fdr_cut=1.0):
    """Peptide-level rollup, rank by score, walk down, return the deepest set with corrected FDR<=cut."""
    best = {}
    for s, pep, ent, *_ in rows:
        if pep not in best or s > best[pep][0]: best[pep] = (s, ent)
    peps = sorted(best.values(), key=lambda v: -v[0])
    n_t = n_e = 0; best_t = 0; best_flags = []; flags = []
    for s, ent in peps:
        if ent: n_e += 1
        else: n_t += 1
        flags.append(1 if ent else 0)
        fdr = 100.0 * (n_e / ratio) / max(n_t, 1)
        if fdr <= fdr_cut and n_t > best_t:
            best_t = n_t; best_flags = list(flags)
    lo, hi = _bootstrap_fdr(best_flags, ratio) if best_flags else (float("nan"),) * 2
    return best_t, len(peps), lo, hi

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default="/path/to/scratch/bench/human_entrap.fasta")
    ap.add_argument("--rt-gap", type=float, default=10.0)
    ap.add_argument("arms", nargs="+", help="name=tsv")
    a = ap.parse_args()
    r = peptide_ratio(a.fasta)
    print(f"peptide ratio = {r:.4f}; RT chain gap = {a.rt_gap}s\n")
    print(f"{'arm':22s} {'PSMs':>9s} {'kept':>9s} {'reduce':>7s} {'targets@1%FDR':>14s}  {'95% CI':>14s}")
    for spec in a.arms:
        name, tsv = spec.split("=", 1)
        rows = load(tsv)
        for rule in ("full", "apex", "best", "top3", "top5", "top10"):
            k = collapse(rows, rule, rt_gap=a.rt_gap)
            t, npep, lo, hi = targets_at_fdr(k, r)
            red = len(rows) / max(len(k), 1)
            print(f"{name+'/'+rule:22s} {len(rows):9d} {len(k):9d} {red:6.2f}x {t:14d}  {lo:6.2f}-{hi:.2f}%")
