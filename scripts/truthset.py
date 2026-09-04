#!/usr/bin/env python3
"""DIA-NN-anchored truth set for pseudo-spectra (v2).

Fixes vs v1, each of which would have biased the comparison:
  * the reference implementation writes the PRECURSOR as its isolation window (+-0.01 Da), we write
    the real acquisition window.  Trusting either file's window makes the reference implementation
    look artificially pure.  Both tools are now scored against ONE common window
    scheme recovered from the acquisition method.
  * the reference implementation emits 500 peaks/spectrum, we emit ~23.  Intensity fractions over
    unequal depth are not comparable, so every number is also reported at equal
    depth (top-N by intensity).
  * precursor m/z and ion mobility live under different keys in the two writers.

Unmatched spectra are reported as UNEXPLAINED, never as noise: DIA-NN's library
is tryptic/2-mod/z2-4 and is structurally blind to the non-canonical species
this tool exists to find.
"""
import sys, collections, numpy as np, pyarrow.parquet as pq
from pyteomics import mzml

REPORT, LIB, SPECTRA, TAG = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
NBINS   = int(sys.argv[5]) if len(sys.argv) > 5 else 200
TOPN    = 20                       # equal-depth comparison
PPM_PREC, PPM_FRAG = 20.0, 20.0
D_RT, D_IM = 0.30, 0.05
DECOY_SHIFT = 11.003

# 26 variable-width diaPASEF windows recovered from the acquisition method.
WINDOWS = [(327.5,468.8),(467.8,500.6),(499.6,524.3),(523.3,545.8),(544.8,564.6),
           (563.6,582.9),(581.8,602.8),(601.8,622.4),(621.4,642.6),(641.6,663.3),
           (662.3,685.3),(684.3,708.3),(707.3,733.0),(732.0,759.9),(758.9,789.5),
           (788.5,822.6),(821.6,860.3),(859.3,903.9),(902.9,955.6),(954.6,1018.6),
           (1017.6,1098.9),(1097.9,1206.0),(1205.0,1364.6),(1363.6,1650.0)]
def window_of(mz):
    for lo, hi in WINDOWS:
        if lo <= mz <= hi: return (lo, hi)
    return None

r = pq.read_table(REPORT, columns=["Precursor.Id","Precursor.Mz","Precursor.Charge",
      "RT","RT.Start","RT.Stop","IM","Q.Value"]).to_pandas()
r = r[r["Q.Value"] <= 0.01].reset_index(drop=True)
lib = pq.read_table(LIB, columns=["Precursor.Id","Product.Mz","Decoy"]).to_pandas()
lib = lib[lib["Decoy"] == 0]
frags = collections.defaultdict(list)
for pid, mz in zip(lib["Precursor.Id"], lib["Product.Mz"]): frags[pid].append(mz)
frags = {k: np.sort(np.array(v)) for k, v in frags.items()}

rt_lo, rt_hi = r["RT.Start"].min(), r["RT.Stop"].max()
NB = 400; edges = np.linspace(rt_lo, rt_hi, NB+1)
dens = np.zeros(NB, int)
for a,b in zip(r["RT.Start"], r["RT.Stop"]):
    dens[max(0,np.searchsorted(edges,a)-1):np.searchsorted(edges,b)] += 1
order = np.argsort(dens); per = max(1, NBINS//10)
sel = np.sort(np.concatenate([order[i*NB//10:(i+1)*NB//10][
        np.linspace(0, max(0,(NB//10)-1), per).astype(int)] for i in range(10)]))
keep = [(edges[i], edges[i+1]) for i in sel]
ka = np.array([k[0] for k in keep]); kb = np.array([k[1] for k in keep])
def in_sel(rt):
    i = np.searchsorted(ka, rt) - 1
    return 0 <= i < len(ka) and rt <= kb[i]
print(f"[truth] {len(r)} precursors q<=1%  [strat] {len(sel)} bins, "
      f"density {dens[sel].min()}-{dens[sel].max()} (global {dens.min()}-{dens.max()})", flush=True)

P = r.to_dict("records")
RTBIN = 0.1
pidx = collections.defaultdict(list)
for i,p in enumerate(P):
    for b in range(int((p["RT.Start"]-D_RT)/RTBIN), int((p["RT.Stop"]+D_RT)/RTBIN)+1):
        pidx[b].append(i)

def hits(arr, mz, ppm):
    if len(arr)==0: return False
    i = np.searchsorted(arr, mz)
    for j in (i-1, i):
        if 0 <= j < len(arr) and abs(arr[j]-mz)/mz*1e6 <= ppm: return True
    return False

def attribute(mzs, ints, fo, others, fo_dec, oth_dec):
    tot = ints.sum() or 1.0; io=ib=iu=idc=iod=0.0
    for mz, it in zip(mzs, ints):
        if   hits(fo, mz, PPM_FRAG):     io += it
        elif hits(others, mz, PPM_FRAG): ib += it
        else:                            iu += it
        if hits(fo_dec, mz, PPM_FRAG):   idc += it
        if hits(oth_dec, mz, PPM_FRAG):  iod += it
    return io/tot, ib/tot, iu/tot, idc/tot, iod/tot

stat=collections.Counter(); nco=[]; A=[[],[],[],[],[]]; T=[[],[],[],[],[]]
zc=collections.Counter(); mppm=[]; matched=set(); n_spec=0; npk=[]
matched_dec=set()
for s in mzml.read(SPECTRA):
    if s.get("ms level") != 2: continue
    pre = s["precursorList"]["precursor"][0]
    ion = pre["selectedIonList"]["selectedIon"][0]
    iso = pre.get("isolationWindow", {})
    pmz = float(ion.get("selected ion m/z",
                iso.get("isolation window target m/z", 0)) or 0)
    if pmz <= 0: continue
    z   = int(ion.get("charge state", 0) or 0)
    sc  = s["scanList"]["scan"][0]
    # DIA-NN reports minutes; OUR mzML is in seconds and the reference implementation's in minutes.
    # Trust the CV unit, never the bare number.
    _rt = sc["scan start time"]
    rt  = float(_rt) / (60.0 if getattr(_rt, "unit_info", "minute") == "second" else 1.0)
    im  = float(ion.get("inverse reduced ion mobility",
                sc.get("inverse reduced ion mobility", -1)) or -1)
    if not in_sel(rt): continue
    n_spec += 1
    mzs, ints = s["m/z array"], s["intensity array"]
    if len(mzs) == 0: continue
    npk.append(len(mzs))
    w = window_of(pmz)
    if w is None: stat["outside_window_scheme"] += 1; continue
    co = [P[i] for i in pidx.get(int(rt/RTBIN), ())
          if w[0] <= P[i]["Precursor.Mz"] <= w[1]
          and (im < 0 or abs(P[i]["IM"]-im) <= D_IM)]
    nco.append(len(co))
    # RECALL DECOY: ~50 true precursors are co-isolated per spectrum, so a 20 ppm
    # precursor match can land by chance. Shifting our m/z off the real grid gives
    # the recall floor; without it 74% vs 84% is uninterpretable.
    dmz = pmz + DECOY_SHIFT
    for p in co:
        if abs(p["Precursor.Mz"]-dmz)/dmz*1e6 <= PPM_PREC: matched_dec.add(p["Precursor.Id"])
    cand = [p for p in co if abs(p["Precursor.Mz"]-pmz)/pmz*1e6 <= PPM_PREC]
    if not cand: stat["unexplained"] += 1; continue
    best = min(cand, key=lambda p: abs(p["Precursor.Mz"]-pmz))
    matched.add(best["Precursor.Id"])
    if z: zc[(z, best["Precursor.Charge"])] += 1
    mppm.append((pmz-best["Precursor.Mz"])/best["Precursor.Mz"]*1e6)
    fo = frags.get(best["Precursor.Id"], np.array([]))
    ol = [frags[p["Precursor.Id"]] for p in co
          if p["Precursor.Id"] != best["Precursor.Id"] and p["Precursor.Id"] in frags]
    others = np.sort(np.concatenate(ol)) if ol else np.array([])
    # SIZE-MATCHED DECOY for CO-ISOLATED: `fo` holds ~9 m/z but `others` holds ~450,
    # so its random-match rate is ~50x higher. The OWN decoy does not control it.
    others_dec = np.sort(others + DECOY_SHIFT) if len(others) else others
    for dst, (m2, i2) in ((A, (mzs, ints)),
                          (T, (lambda k: (mzs[k], ints[k]))(np.argsort(ints)[-TOPN:]))):
        for lst, v in zip(dst, attribute(m2, i2, fo, others, fo+DECOY_SHIFT, others_dec)):
            lst.append(v)
    o, b = A[0][-1], A[1][-1]
    stat["pure" if o>=0.5 and b<0.10 else "chimeric" if b>=0.10 else "low_signal"] += 1

present = {p["Precursor.Id"] for p in P if in_sel(p["RT"])}
print(f"\n===== {TAG} =====")
print(f"spectra in selected cycles : {n_spec}   median peaks/spectrum {int(np.median(npk)) if npk else 0}")
print(f"RECALL  {len(matched & present)}/{len(present)} = "
      f"{100*len(matched&present)/max(1,len(present)):.1f}% of DIA-NN precursors")
print(f"RECALL decoy floor (our m/z +{DECOY_SHIFT} Da): "
      f"{len(matched_dec & present)}/{len(present)} = "
      f"{100*len(matched_dec&present)/max(1,len(present)):.1f}%")
print(f"class   {dict(stat)}")
if nco:
    n = np.array(nco)
    print(f"CO-ISOLATED truth precursors/spectrum: mean {n.mean():.1f} median {np.median(n):.0f} "
          f"p90 {np.percentile(n,90):.0f} max {n.max()} | {100*(n<=1).mean():.1f}% singleton")
for name, D in (("ALL PEAKS", A), (f"TOP-{TOPN} (equal depth)", T)):
    if D[0]:
        o,b,u,d,od = map(np.array, D)
        print(f"{name:22s} OWN {o.mean()*100:5.2f}% (decoy {d.mean()*100:.2f}%, "
              f"net {(o.mean()-d.mean())*100:5.2f}%, {o.mean()/max(d.mean(),1e-9):4.1f}x)  |  "
              f"CO-ISO {b.mean()*100:5.2f}% (size-matched decoy {od.mean()*100:.2f}%, "
              f"net {(b.mean()-od.mean())*100:+5.2f}%)  |  UNEXPL {u.mean()*100:5.1f}%")
if zc:
    ag = sum(v for (a,bq),v in zc.items() if a==bq); tt = sum(zc.values())
    print(f"CHARGE agreement {ag}/{tt} = {100*ag/tt:.1f}%   "
          f"confusions {sorted(((v,f'{a}->{bq}') for (a,bq),v in zc.items() if a!=bq), reverse=True)[:4]}")
if mppm:
    m = np.array(mppm)
    print(f"PRECURSOR MASS error vs truth: median {np.median(m):+.2f} ppm  "
          f"IQR [{np.percentile(m,25):+.2f},{np.percentile(m,75):+.2f}]  "
          f"|err|>10ppm {100*(abs(m)>10).mean():.1f}%")
