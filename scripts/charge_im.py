#!/usr/bin/env python3
"""Can (m/z, ion mobility) predict charge?

The funnel says 15.9% of precursors get the WRONG CHARGE and only 2.2% are undetected.
The dominant confusion is 4->2: at z=4 the isotope spacing is 0.25 Da vs 0.50 at z=2, so
a spurious partner halfway between real isotopes reads as a higher charge. Isotope
spacing alone cannot break that tie.

But diaPASEF measures mobility, and peptide charge is strongly determined by (m/z, 1/K0)
-- the acquisition scheme's own tiles track the charge-2 line. If charge is separable in
that plane, it is a physics constraint we are simply not using.

Ground truth: DIA-NN's 43,499 identified precursors, which carry m/z, IM and charge.
"""
import sys, numpy as np, pyarrow.parquet as pq

r = pq.read_table(sys.argv[1], columns=["Precursor.Mz","Precursor.Charge","IM","Q.Value"]).to_pandas()
r = r[r["Q.Value"] <= 0.01]
mz, z, im = r["Precursor.Mz"].values, r["Precursor.Charge"].values, r["IM"].values
print(f"{len(r)} precursors  charges {dict(zip(*np.unique(z, return_counts=True)))}\n")

# 1. Is there a clean charge-vs-mobility band at fixed m/z?
print("IM (1/K0) by charge, in m/z slices — overlap is what decides separability")
print(f"{'m/z slice':>12s} " + "".join(f"{'z='+str(c):>22s}" for c in (2,3,4)))
for lo, hi in ((400,500),(500,600),(600,700),(700,800),(800,1000)):
    m = (mz>=lo)&(mz<hi)
    row = f"{lo}-{hi:<7d}"
    for c in (2,3,4):
        s = im[m & (z==c)]
        row += f"{(f'{s.mean():.3f}+-{s.std():.3f} n={len(s)}' if len(s)>20 else '-'):>22s}"
    print(row)

# 2. How much does a charge error cost in NEUTRAL MASS?  dM = (z'-z)(mz - 1.00728)
print("\nneutral-mass error from a charge misassignment (open search sees a phantom mod):")
for a, b in ((4,2),(2,3),(3,2),(1,2)):
    d = (b-a)*(600.0-1.00728)
    print(f"  ours z={a} vs true z={b} at m/z 600: dM = {d:+.1f} Da")

# 3. Nearest-neighbour separability in the (m/z, IM) plane: for each precursor, does its
#    nearest neighbour in that plane share its charge?  A high rate means the plane alone
#    carries charge information -- no model needed to establish that.
from collections import Counter
idx = np.argsort(mz)
mzs, zs, ims = mz[idx], z[idx], im[idx]
agree = Counter()
step = max(1, len(mzs)//20000)          # subsample; O(n^2) otherwise
for i in range(0, len(mzs), step):
    j0, j1 = np.searchsorted(mzs, mzs[i]-2.0), np.searchsorted(mzs, mzs[i]+2.0)
    if j1-j0 < 2: continue
    d = np.abs(ims[j0:j1] - ims[i]); d[i-j0] = 1e9      # exclude self
    nn = j0 + int(np.argmin(d))
    agree[zs[nn] == zs[i]] += 1
tot = sum(agree.values())
print(f"\nnearest neighbour in (m/z +-2 Da, closest IM) shares charge: "
      f"{100*agree[True]/max(tot,1):.1f}%  (n={tot})")
# NULL CORRECTION. This test asks whether a precursor's nearest neighbour SHARES its charge --
# a coincidence rate between two draws, whose null is sum(p_z^2), NOT max(p_z). max(p_z) is the
# null for PREDICTING one label. Using it inflated the baseline from 55.8% to 69.6% and turned a
# +19.2-point signal into +5.4, which is the sentence that falsified the whole ion-mobility
# lever in docs/charge-inference.md. Both nulls are printed so the error cannot recur silently.
_c = np.bincount(z)[np.bincount(z) > 0]
_p = _c / _c.sum()
print(f"null, PREDICT one label   (max p_z) : {100*_p.max():.1f}%   <- WRONG for this test")
print(f"null, SHARE a label     (sum p_z^2) : {100*(_p**2).sum():.1f}%   <- correct")
