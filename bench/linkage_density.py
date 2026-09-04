#!/usr/bin/env python3
"""Per-spectrum ID linkage (#2) -> cycle-density cross (#4) + long/high-mass diagnosis (#1).

v2: dumps the parsed population to a TSV (reused, so the 14.7 GB mzML is streamed once), and
computes co-isolation density in O(N) by binning into (RT-cycle, IM-band, m/z-tile) cells instead
of an O(N*window) sweep. The cell IS the acquisition unit: one diaPASEF cycle x one IM ramp step x
one ~25 Th isolation tile. Spectra in the same cell were packed together by the method -- that is
the co-isolation the "high-density cycles" question is about.

Same honest limits as v1: "identified" = Sage peptide_q<=0.01 on that spectrum id; density is a
proxy for MS1 co-isolation, not a measurement of it; ID-rate-vs-density is confounded by elution
zone (busy cells are mid-gradient where the real peptides are), so the effect is read WITHIN a
fixed m/z band to strip most of that confound.
"""
import csv, sys, os, collections

MZML = sys.argv[1]
SAGE = sys.argv[2]
POP = "/path/to/scratch/pop_s30.tsv"

# ---- 1. population: parse once, cache to TSV ------------------------------------------------
pop = {}
if os.path.exists(POP):
    with open(POP) as f:
        for ln in f:
            idx, rt, mz, im, npk = ln.split("\t")
            pop[int(idx)] = [None if rt == "" else float(rt), None if mz == "" else float(mz),
                             None if im == "" else float(im), int(npk)]
    print("[pop] %d emitted spectra (from cache)" % len(pop), flush=True)
else:
    cur = rt = mz = im = npk = None
    n = 0
    out = open(POP, "w")
    def flush(idx):
        out.write("%d\t%s\t%s\t%s\t%d\n" % (idx, rt if rt is not None else "",
                  mz if mz is not None else "", im if im is not None else "", npk or 0))
    with open(MZML) as f:
        for line in f:
            if "<spectrum " in line and 'id="spectrum=' in line:
                if cur is not None:
                    pop[cur] = [rt, mz, im, npk]; flush(cur)
                i = line.index('id="spectrum=') + 13
                cur = int(line[i:line.index('"', i)])
                j = line.find('defaultArrayLength="')
                npk = int(line[j+20:line.index('"', j+20)]) if j >= 0 else 0
                rt = mz = im = None
                n += 1
                if n % 500000 == 0:
                    print("  ...%d parsed" % n, flush=True)
            elif cur is not None:
                if mz is None and 'name="selected ion m/z"' in line:
                    k = line.index('value="') + 7; mz = float(line[k:line.index('"', k)])
                elif rt is None and 'name="scan start time"' in line:
                    k = line.index('value="') + 7; rt = float(line[k:line.index('"', k)])
                elif im is None and 'name="inverse reduced ion mobility"' in line:
                    k = line.index('value="') + 7; im = float(line[k:line.index('"', k)])
        if cur is not None:
            pop[cur] = [rt, mz, im, npk]; flush(cur)
    out.close()
    print("[pop] %d emitted spectra (parsed, cached to %s)" % (len(pop), POP), flush=True)

# ---- 2. identified subset -------------------------------------------------------------------
ident = {}
with open(SAGE) as f:
    for x in csv.DictReader(f, delimiter="\t"):
        try:
            if float(x["peptide_q"]) > 0.01:
                continue
            ident[int(x["scannr"].split("=")[-1])] = int(x["peptide_len"])
        except (ValueError, KeyError):
            continue
print("[ident] %d identified (%.1f%% of emitted)" % (len(ident), 100.0*len(ident)/max(len(pop),1)),
      flush=True)

# ---- 3. co-isolation density in O(N): count spectra per acquisition cell ---------------------
RTW, IMW, MZW = 0.03, 0.03, 12.5
cell = collections.Counter()
key = {}
for idx, v in pop.items():
    if v[0] is None or v[1] is None:
        continue
    c = (round(v[0] / RTW), round((v[2] or 0) / IMW), int(v[1] // MZW))
    key[idx] = c
    cell[c] += 1
dens = {idx: cell[c] - 1 for idx, c in key.items()}   # siblings in the same cell
print("[density] %d cells, max occupancy %d" % (len(cell), max(cell.values())), flush=True)


def rate_by(keyfn, order, label, sub=None):
    tot = collections.Counter(); hit = collections.Counter()
    for idx, v in pop.items():
        if v[0] is None or v[1] is None or (sub and not sub(idx, v)):
            continue
        b = keyfn(idx, v)
        if b is None:
            continue
        tot[b] += 1
        if idx in ident:
            hit[b] += 1
    print("\n=== ID rate by %s ===" % label)
    print("%-18s %10s %10s %8s" % ("bin", "emitted", "id'd", "ID rate"))
    for nm in order:
        t, h = tot[nm], hit[nm]
        print("%-18s %10d %10d %7.1f%%" % (nm, t, h, 100.0*h/t if t else 0.0))


def dens_bin(idx, v):
    c = dens.get(idx, 0)
    for lo, hi, nm in [(0,1,"0(alone)"),(1,3,"1-2"),(3,6,"3-5"),(6,11,"6-10"),
                       (11,21,"11-20"),(21,10**9,"21+")]:
        if lo <= c < hi:
            return nm
    return None


ORD = ["0(alone)","1-2","3-5","6-10","11-20","21+"]
rate_by(dens_bin, ORD, "co-isolation density (siblings in acquisition cell)")
rate_by(dens_bin, ORD, "co-isolation density | m/z 700-1200 (fixed band, strips elution confound)",
        sub=lambda idx, v: v[1] is not None and 700 <= v[1] < 1200)
rate_by(dens_bin, ORD, "co-isolation density | m/z 500-700 (fixed band)",
        sub=lambda idx, v: v[1] is not None and 500 <= v[1] < 700)
print("\n(ID rate falling with density AT FIXED m/z = co-isolation hurts beyond elution zone.)")
