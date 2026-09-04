#!/usr/bin/env python3
"""Attribute precursor loss to a stage.

Measured recall is 74.3% vs the reference implementation's 84.4%, but that single number cannot say
WHERE a precursor dies. Three stages, each independently fixable by different code:

  A  MS1 TRACE      is there any MS1 trace at the precursor's (m/z, RT, IM)?
                    failure => detectTraces_ sensitivity (noise/snr/min_length)
  B  HYPOTHESIS     did inferPrecursors_ turn that trace into a precursor with the
                    right mono m/z?  and with the right CHARGE?
                    failure => isotope/charge inference
  C  SPECTRUM       did the window loop emit a spectrum for it?  (measured separately)
                    failure => assembly/scoring

A precursor is "seen" at stage A if a trace sits at the mono m/z OR at the +1/+2
isotope m/z — a trace on the wrong isotope still proves the ion was detected, and
distinguishes "invisible" from "misassigned".
"""
import sys, collections, numpy as np, pyarrow.parquet as pq

REPORT, PREFIX = sys.argv[1], sys.argv[2]
PPM, D_RT, D_IM = 20.0, 0.30, 0.05
PROTON, NEUTRON = 1.007276, 1.003355

r = pq.read_table(REPORT, columns=["Precursor.Id","Precursor.Mz","Precursor.Charge",
      "RT","RT.Start","RT.Stop","IM","Q.Value"]).to_pandas()
r = r[r["Q.Value"] <= 0.01].reset_index(drop=True)
print(f"[truth] {len(r)} precursors q<=1%", flush=True)

def load(path, cols):
    import csv
    out = []
    with open(path) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd: out.append(tuple(float(row[c]) for c in cols))
    return np.array(out, dtype=float)

tr = load(PREFIX + ".traces.tsv",     ["mz","rt","im"])
pc = load(PREFIX + ".precursors.tsv", ["mono_mz","charge","rt","im"])
print(f"[ours]  {len(tr)} MS1 traces, {len(pc)} precursor hypotheses", flush=True)

# RT-bucketed m/z-sorted index; RT in the dumps is SECONDS, DIA-NN reports MINUTES
bkeys = {}
def build(arr, rt_col):
    idx = collections.defaultdict(list)
    for i, row in enumerate(arr):
        idx[int(row[rt_col] / 60.0 / D_RT)].append(i)
    kk = {}
    for k in idx:
        idx[k] = sorted(idx[k], key=lambda i: arr[i][0])
        kk[k] = [arr[i][0] for i in idx[k]]     # m/z keys for bisect
    bkeys[id(arr)] = kk
    return idx
ti, pi = build(tr, 1), build(pc, 2)

import bisect
def near(arr, idx, mz, rt_min, im, ppm=PPM, want_charge=None):
    """any row within ppm of mz, in an RT bucket covering rt_min, and within D_IM."""
    tol = mz * ppm * 1e-6
    for b in (int(rt_min/D_RT)-1, int(rt_min/D_RT), int(rt_min/D_RT)+1):
        bucket = idx.get(b)
        if not bucket: continue
        keys = bkeys[id(arr)].get(b)
        lo = bisect.bisect_left(keys, mz - tol)
        hi = bisect.bisect_right(keys, mz + tol)
        for i in bucket[lo:hi]:
            row = arr[i]
            rt_i = row[1 if arr is tr else 2] / 60.0
            im_i = row[2 if arr is tr else 3]
            if abs(rt_i - rt_min) > D_RT or abs(im_i - im) > D_IM: continue
            if want_charge is not None and int(row[1]) != want_charge: continue
            return True
    return False

st = collections.Counter()
for _, p in r.iterrows():
    mz, z, rt, im = p["Precursor.Mz"], int(p["Precursor.Charge"]), p["RT"], p["IM"]
    iso = [mz + k*NEUTRON/z for k in (0,1,2)]
    has_trace = any(near(tr, ti, m, rt, im) for m in iso)
    has_mono  = near(pc, pi, mz, rt, im)
    has_zmono = near(pc, pi, mz, rt, im, want_charge=z)
    if   not has_trace: st["A_no_ms1_trace"] += 1
    elif not has_mono:  st["B_trace_but_no_hypothesis"] += 1
    elif not has_zmono: st["B_hypothesis_wrong_charge"] += 1
    else:               st["C_correct_hypothesis"] += 1

n = sum(st.values())
print(f"\n===== MS1 FUNNEL, {n} DIA-NN precursors =====")
for k in ("A_no_ms1_trace","B_trace_but_no_hypothesis","B_hypothesis_wrong_charge","C_correct_hypothesis"):
    print(f"  {k:28s} {st[k]:7d}  {100*st[k]/n:5.1f}%")
seen = n - st["A_no_ms1_trace"]
print(f"\nMS1 detection ceiling   : {100*seen/n:.1f}%  (trace exists at mono or +1/+2)")
print(f"survives charge inference: {100*st['C_correct_hypothesis']/n:.1f}%")
print(f"lost in inference, not detection: {100*(st['B_trace_but_no_hypothesis']+st['B_hypothesis_wrong_charge'])/n:.1f}%")
