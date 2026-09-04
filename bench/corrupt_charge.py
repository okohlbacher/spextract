#!/usr/bin/env python3
"""Inject KNOWN charge errors into an mzML, so the open-search instrument has a real positive control.

Why
---
The entrapment open-search benchmark rests on a positive control that did not work. The
`envelope` arm has 10,199 charge confusions measured against a DIA-NN reference, and only **26**
surfaced as slope-2 ray mass in the delta-mass geometry. Two readings are possible and the data
cannot separate them:

  (a) the instrument cannot see charge errors -> every number from it is uninterpretable
  (b) a wrong precursor mass usually prevents identification entirely, so the errors are real
      but invisible by construction -- only 1,175 of 1.35M raw PSMs survived 1% entrapment FDR

Both reviewers independently said: stop arguing, corrupt a known fraction and measure recovery.

What this does
--------------
Takes a WORKING mzML and rewrites the precursor charge of a random `--fraction` of spectra from
z to z', leaving the m/z untouched. That is exactly the error the tool makes. The neutral mass
the search engine then computes is wrong by

    dM = (z' - z) * (m/z - PROTON)

so those spectra MUST appear on the slope-(z'-z) ray if the instrument works.

The recovery rate is the instrument's sensitivity, and it is measurable because the ground truth
is written by this script:
    recovered / injected  ->  what fraction of charge errors the assay can actually see

If recovery is high, reading (a) is dead and the assay is usable. If recovery is near zero, the
assay cannot detect charge errors at all and NO conclusion from it may be quoted -- including
the 1.05x modified-peptide ratio, which both reviewers flagged as resting on this control.

Deliberately NOT done: shifting m/z as well. That would change which peptide matches and confound
"can the assay see charge errors" with "does the spectrum still identify".
"""
import argparse, random, re, sys
from pathlib import Path

PROTON = 1.00727646
CHG = re.compile(r'(name="charge state" value=")(\d+)(")')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--fraction", type=float, default=0.10,
                    help="fraction of charge-bearing spectra to corrupt")
    ap.add_argument("--from-z", type=int, default=2)
    ap.add_argument("--to-z", type=int, default=4,
                    help="z'; the induced ray slope is (to_z - from_z)")
    ap.add_argument("--truth", help="TSV of injected corruptions (spectrum index, old z, new z)")
    ap.add_argument("--seed", type=int, default=17)
    a = ap.parse_args()
    random.seed(a.seed)

    n_spec = n_elig = n_hit = 0
    truth = open(a.truth, "w") if a.truth else None
    if truth:
        truth.write("spectrum_index\told_z\tnew_z\n")

    with open(a.src, "r", encoding="utf-8", errors="replace") as f, \
         open(a.dst, "w", encoding="utf-8") as g:
        for line in f:
            if "<spectrum " in line:
                n_spec += 1
            m = CHG.search(line)
            if m and int(m.group(2)) == a.from_z:
                n_elig += 1
                if random.random() < a.fraction:
                    line = CHG.sub(r"\g<1>%d\g<3>" % a.to_z, line, count=1)
                    n_hit += 1
                    if truth:
                        truth.write("%d\t%d\t%d\n" % (n_spec, a.from_z, a.to_z))
            g.write(line)
    if truth:
        truth.close()

    if n_elig == 0:
        sys.exit("ABORT: no spectra with charge %d found in %s -- wrong file or wrong --from-z"
                 % (a.from_z, a.src))
    print("spectra=%d  eligible(z=%d)=%d  corrupted=%d (%.1f%% of eligible)"
          % (n_spec, a.from_z, n_elig, n_hit, 100.0 * n_hit / n_elig))
    print("induced ray slope k = %+d, so dM = %+d * (m/z - %.5f)"
          % (a.to_z - a.from_z, a.to_z - a.from_z, PROTON))
    print("at m/z 600 that is dM = %+.1f Da -- OUTSIDE a +-500 window, so the WIDE config is required"
          % ((a.to_z - a.from_z) * (600.0 - PROTON)))
    if a.truth:
        print("ground truth -> %s" % a.truth)


if __name__ == "__main__":
    sys.exit(main())
