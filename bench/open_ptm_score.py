#!/usr/bin/env python3
"""Delta-mass-family scoring of open-search output, v3 — rebuilt to the 2026-09-01 step-01
adversarial review (kimi F1-F9, codex #1-#18).

What changed vs v2 and why:
- NEAREST-CANDIDATE assignment (candidates = 0, isotope lattice +-k*1.0033548 k=1..15, curated KNOWN
  masses) with a ppm-SCALED tolerance, instead of ordered fixed-window checks: kills the
  didehydro/iso-2 ordering bias [codex #8], the k<=9 cutoff arbitrariness [kimi F2, codex #10] and
  the Da-constant tolerance mass bias [kimi F4, codex #12].
- Pooled per-FAMILY walks only (unmod / nearzero / isotope / knownPTM / other), each reporting
  accepted targets AND entrapments, achieved FDR, score cutoff, and a 95% bootstrap CI, with a
  conservative (e+1)/ratio numerator [codex 'Direct decisions']; tie-GROUPED walk (a cutoff never
  splits a tied score) [codex #17].
- Per-bin numbers are DESCRIPTIVE ONLY and printed with a minimum-evidence rule: an @1% claim
  appears only if the accepted prefix contains >=10 entrapment observations; otherwise units + raw
  entrap fraction only [kimi F5/F6, codex #3].
- Peptide-level DEDUP union report across unmod+nearzero (a peptide may hold several precursor
  hypotheses) [kimi B4, codex #11].
- Provenance header: script+input sha256, invocation, estimator ratio [codex #16].
CAVEATS the output itself carries: deamidation (+0.9840, shared with citrullination) vs iso+1
(+1.0034) are 19 mDa apart and NOT robustly separable at this data's precursor accuracy — they are
reported as one ambiguous pair inside their families; entrapment FDR certifies the PEPTIDE, not the
DELTA [kimi F6] — no per-bin walk can legitimate a delta value; and no site localization is done, so
'knownPTM' counts mass-compatible peptide hypotheses, not localized PTM sites [codex #1].
"""
import sys, os, csv, argparse, hashlib, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from entrapment import peptide_ratio, _bootstrap_fdr

ISO = 1.0033548
KNOWN = [
    (15.9949, "Oxidation"), (31.9898, "Dioxidation"), (47.9847, "Trioxidation"),
    (42.0106, "Acetyl"), (79.9663, "Phospho"), (0.98402, "Deamid/Citrullination[iso+1-ambig]"),
    (-0.98402, "Amidation[iso-1-ambig]"), (14.0157, "Methyl"), (28.0313, "Dimethyl"),
    (42.0470, "Trimethyl"), (27.9949, "Formyl"), (43.0058, "Carbamyl"),
    (57.0215, "Carbamidomethyl(offtarget)"), (-57.0215, "UnalkylatedCys"),
    (-17.0265, "PyroGlu(Q)"), (-18.0106, "Dehydration/PyroGlu(E)"), (-2.01565, "Didehydro[iso-2-ambig]"),
    (21.9819, "Cation:Na"), (37.9559, "Cation:K"), (52.9115, "Cation:Fe"),
    (114.0429, "GG(ubiquitin)"), (119.0041, "Cysteinylation"), (100.0160, "Succinyl"),
    (26.0157, "Acetaldehyde"), (-48.0034, "Dethiomethyl"), (203.0794, "HexNAc"), (162.0528, "Hexose"),
    (-131.0405, "Met-loss"), (-89.0299, "Met-loss+Acetyl"),
]
CAND = [(0.0, "unmod", "unmod")] \
     + [(s * k * ISO, f"iso{'+' if s>0 else '-'}{k}", "isotope") for k in range(1, 16) for s in (1, -1)] \
     + [(m, n, "knownPTM") for m, n in KNOWN]

def tol_of(mass, ppm, floor):
    return max(floor, mass * ppm * 1e-6)

def classify(d, mass, ppm, floor):
    tol = tol_of(mass, ppm, floor)
    delta_best, name, fam = min(((abs(d - c), n, f) for c, n, f in CAND), key=lambda t: t[0])
    if delta_best <= tol: return (fam, name)
    if abs(d) < 0.75: return ("nearzero", "masserr")
    return ("other", f"{d:+.0f}Da")

def sha(path, cap=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()[:16]

def load(tsv, ppm, floor, tag="ENTRAP_"):
    best = {}   # (peptide, name) -> (score, ent, fam)
    pep_best = {}  # peptide -> (score, ent, fam, name)  for dedup reporting
    with open(tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                if int(x["rank"]) != 1: continue
                prots = [p for p in x["proteins"].split(";") if p]
                if not prots or all(p.startswith("rev_") for p in prots): continue
                s = float(x["sage_discriminant_score"])
                m = float(x["calcmass"]); d = float(x["expmass"]) - m
            except (ValueError, KeyError): continue
            fam, name = classify(d, m, ppm, floor)
            ent = all(p.startswith(tag) or p.startswith("rev_" + tag) for p in prots)
            k = (x["peptide"], name)
            if k not in best or s > best[k][0]: best[k] = (s, ent, fam)
            if x["peptide"] not in pep_best or s > pep_best[x["peptide"]][0]:
                pep_best[x["peptide"]] = (s, ent, fam, name)
    return best, pep_best

def walk(items, ratio, fdr_cut):
    """Tie-grouped conservative walk: items [(score, ent)]. FDR = ((e+1)/ratio)/t, evaluated only at
    tie-group boundaries. Returns (t, e, achieved, cutoff_score, flags_of_accepted)."""
    from itertools import groupby
    items = sorted(items, key=lambda v: -v[0])
    t = e = 0; acc = (0, 0, float("nan"), float("nan"), [])
    flags = []
    for score, grp in groupby(items, key=lambda v: v[0]):
        for _, ent in grp:
            if ent: e += 1
            else: t += 1
            flags.append(1 if ent else 0)
        fdr = 100.0 * ((e + 1) / ratio) / max(t, 1)
        if fdr <= fdr_cut and t > acc[0]: acc = (t, e, fdr, score, list(flags))
    return acc

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default="/path/to/scratch/bench/human_entrap.fasta")
    ap.add_argument("--ppm", type=float, default=12.0)
    ap.add_argument("--floor", type=float, default=0.008)
    ap.add_argument("--fdr", type=float, default=1.0)
    ap.add_argument("arms", nargs="+", help="name=tsv")
    a = ap.parse_args()
    r = peptide_ratio(a.fasta)
    print(f"# open_ptm_score v3  {datetime.datetime.now().isoformat(timespec='seconds')}")
    print(f"# invocation: {' '.join(sys.argv)}")
    print(f"# script sha256/16: {sha(os.path.abspath(__file__))}   estimator ratio: {r:.4f}")
    print(f"# tol = max({a.floor} Da, {a.ppm} ppm x calcmass); conservative (e+1)/r walk; tie-grouped")
    for spec in a.arms:
        name, tsv = spec.split("=", 1)
        print(f"\n=== {name} ===  input {tsv}  sha256/16 {sha(tsv)}")
        best, pep_best = load(tsv, a.ppm, a.floor)
        for fam in ("unmod", "nearzero", "isotope", "knownPTM", "other"):
            items = [(s, e) for (s, e, f) in best.values() if f == fam]
            t, e, ach, cut, flags = walk(items, r, a.fdr)
            lo, hi = _bootstrap_fdr(flags, r) if flags else (float("nan"),) * 2
            print(f"  {fam:9s} @{a.fdr}%: {t:>6} targets  (e={e}, achieved {ach:.2f}%, CI [{lo:.2f}-{hi:.2f}], cutoff {cut:.3f}, units {len(items)})")
        # peptide-level dedup unions (a peptide may carry several precursor hypotheses)
        uz = {p for p, (s, e, f, n) in pep_best.items() if f in ("unmod", "nearzero") and not e}
        print(f"  unique peptides with best hypothesis in unmod|nearzero (target strand): {len(uz)}")
        # descriptive per-name table with minimum-evidence rule
        from collections import defaultdict
        by_name = defaultdict(list)
        for (pep, n), (s, e, f) in best.items():
            if f in ("knownPTM", "isotope", "other"): by_name[f"{f}:{n}"].append((s, e))
        print("  per-name (DESCRIPTIVE; @1% shown only when accepted entrapment >= 10):")
        for n, items in sorted(by_name.items(), key=lambda kv: -len(kv[1]))[:20]:
            t, e, ach, cut, flags = walk(items, r, a.fdr)
            ef = 100.0 * sum(1 for _, x in items if x) / len(items)
            claim = f"@1%:{t:>5} (e={e})" if e >= 10 else f"@1%: n/a (insufficient entrapment evidence, e={e})"
            print(f"    {n:34s} {claim}   units {len(items):>6}, raw entrap {ef:.1f}%")
