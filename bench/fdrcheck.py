#!/usr/bin/env python3
"""Is the min_correlation gain REAL SIGNAL or an FDR RECALIBRATION artefact?

The adversarial review's top BLOCKER: peptides rise (+7.7%) while PSMs FALL (-6.6%) and
fdr_loss_pct shifts. Sage computes peptide_q from ITS OWN decoys per run, so q<=0.01 is a
DIFFERENT score cut in every arm -- the counts may not be comparable at all.

Decisive test, needing no new runs: apply ONE FIXED discriminant-score threshold to every arm.
  * gain persists at fixed score  -> REAL: the arm genuinely produces better-scoring spectra
  * gain vanishes at fixed score  -> ARTEFACT: same signal, looser effective cut

Also reports each arm's own 1%-FDR score threshold and its target/decoy counts, so a shift in the
decoy distribution is visible rather than inferred. A rising decoy count with a falling threshold
is the signature of recalibration.
"""
import csv, sys, os


def load(tsv):
    """best discriminant score per (peptide, is_decoy)."""
    best = {}
    with open(tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                pep = x["peptide"]
                s = float(x["sage_discriminant_score"])
                q = float(x["peptide_q"])
                prots = [p for p in x["proteins"].split(";") if p]
                dec = bool(prots) and all(p.startswith("rev_") for p in prots)
            except (ValueError, KeyError):
                continue
            k = (pep, dec)
            if k not in best or s > best[k][0]:
                best[k] = (s, q)
    return best


FIX = [1.0, 1.5, 2.0, 2.5]
print("%-9s %8s %8s %11s %10s   %s"
      % ("arm", "targets", "decoys", "q<=.01 thr", "own count", "peptides above FIXED score"))
rows = {}
for arm in sys.argv[1:]:
    tsv = os.path.join(arm, "results.sage.tsv")
    if not os.path.exists(tsv):
        continue
    b = load(tsv)
    name = os.path.basename(arm).replace("S30__", "")
    tgt = [(s, q) for (p, d), (s, q) in b.items() if not d]
    dec = [(s, q) for (p, d), (s, q) in b.items() if d]
    own = sum(1 for s, q in tgt if q <= 0.01)
    thr = min([s for s, q in tgt if q <= 0.01], default=float("nan"))
    fixed = [sum(1 for s, q in tgt if s >= f) for f in FIX]
    rows[name] = (len(tgt), len(dec), thr, own, fixed)
    print("%-9s %8d %8d %11.4f %10d   %s"
          % (name, len(tgt), len(dec), thr, own,
             "  ".join("s>=%.1f:%d" % (f, c) for f, c in zip(FIX, fixed))))

if "base" in rows:
    b = rows["base"]
    print("\nDELTA vs base --  q-based (each arm's own threshold)  vs  FIXED-score (common cut):")
    for n in rows:
        if n == "base":
            continue
        t, d, thr, own, fx = rows[n]
        dq = 100.0 * (own - b[3]) / b[3]
        df = "  ".join("s>=%.1f %+.1f%%" % (f, 100.0 * (c - bc) / bc)
                       for f, c, bc in zip(FIX, fx, b[4]))
        print("  %-9s q-based %+6.1f%%   |  %s" % (n, dq, df))
    print("\nREAD: fixed-score deltas tracking the q-based delta  -> the gain is REAL.")
    print("      fixed-score deltas ~0 while q-based is large    -> RECALIBRATION artefact.")
    print("      decoy count and threshold columns show whether the null distribution moved.")
