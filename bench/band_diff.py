#!/usr/bin/env python3
"""How much do two band counts actually disagree?

`perf:trace_bands` is documented as an exact partition ("traces cannot span more than
trace:mass_error_ppm in m/z, so banding with a halo is exact"), and the semantic digest says
otherwise. A digest mismatch is binary, so it cannot tell a one-peak rounding difference from a
systematic loss. This reports the magnitude at three levels:

  spectra    how many precursor coordinates exist in one run and not the other
  content    for the shared coordinates, how many have different peak lists
  peptides   the identified-peptide sets, if Sage results are given

Usage: band_diff.py A.tsv B.tsv [A.peps B.peps]
       (the .tsv files come from bench/iso_dup.py extract)
"""
import sys, collections


def load(p):
    """coordinate -> (npeaks, count). The coordinate is what a search engine keys a PSM on."""
    d = collections.Counter()
    npk = {}
    with open(p) as f:
        next(f)
        for line in f:
            q = line.rstrip("\n").split("\t")
            key = (round(float(q[2]), 5), round(float(q[1]), 3), int(q[3]), round(float(q[4]), 5))
            d[key] += 1
            npk[key] = int(q[5])
    return d, npk


def main(a, b, pa=None, pb=None):
    da, na = load(a)
    db, nb = load(b)
    only_a = sum(v for k, v in da.items() if k not in db)
    only_b = sum(v for k, v in db.items() if k not in da)
    shared = [k for k in da if k in db]
    diff_n = sum(1 for k in shared if na[k] != nb[k])
    ta, tb = sum(da.values()), sum(db.values())
    print(f"{'':22}{'A':>12}{'B':>12}")
    print(f"{'spectra':22}{ta:>12,}{tb:>12,}")
    print(f"{'distinct coordinates':22}{len(da):>12,}{len(db):>12,}")
    print(f"\nonly in A: {only_a:,} spectra ({100.0*only_a/max(ta,1):.4f}% of A)")
    print(f"only in B: {only_b:,} spectra ({100.0*only_b/max(tb,1):.4f}% of B)")
    print(f"shared coordinates: {len(shared):,}; of those {diff_n:,} "
          f"({100.0*diff_n/max(len(shared),1):.4f}%) differ in peak count")
    if pa and pb:
        sa = set(open(pa).read().split())
        sb = set(open(pb).read().split())
        print(f"\npeptides A {len(sa):,}  B {len(sb):,}  shared {len(sa & sb):,}  "
              f"only-A {len(sa - sb):,}  only-B {len(sb - sa):,}")
        u = len(sa | sb)
        print(f"Jaccard {100.0*len(sa & sb)/max(u,1):.3f}%  "
              f"net {len(sb) - len(sa):+,} ({100.0*(len(sb)-len(sa))/max(len(sa),1):+.3f}%)")


if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    main(*sys.argv[1:5])
