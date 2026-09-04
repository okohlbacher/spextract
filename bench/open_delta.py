#!/usr/bin/env python3
# Open-search delta-mass diagnostic. delta = expmass - calcmass (the mass shift the open search found).
# For a TRUSTWORTHY extractor the histogram is dominated by 0 (unmodified) + sharp peaks at real PTM
# masses; PHANTOM structure (isotope miscalls at +-1.00335, charge-error rays) = extraction artifacts
# that fabricate modifications in blind search. Compares one arm; run per file and eyeball side by side.
import sys, csv, collections
KNOWN = {0.0:"unmod", 15.9949:"oxid(+16)", 42.0106:"acetyl(+42)", 79.9663:"phospho(+80)",
         0.9840:"deamid(+1)", -17.0265:"pyroGlu(-17)", 57.0215:"carbamidomethyl", 1.00335:"iso+1",
         -1.00335:"iso-1", 2.00671:"iso+2", 21.9819:"Na-H", 43.0058:"carbamyl(+43)", 100.0:"other"}
def load(path, qcut=0.01):
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t"); ci={c:i for i,c in enumerate(h)}
    exp,calc = ci["expmass"], ci["calcmass"]
    qf = ci.get("spectrum_q", ci.get("peptide_q"))
    lbl = ci.get("label"); iso = ci.get("isotope_error")
    deltas=[]; n=0
    with open(path) as f:
        r=csv.reader(f,delimiter="\t"); next(r)
        for row in r:
            if lbl is not None and row[lbl]!="1": continue
            if float(row[qf])>qcut: continue
            d=float(row[exp])-float(row[calc]); deltas.append(d); n+=1
    return deltas
def summarize(name, deltas):
    n=len(deltas)
    if not n: print(f"{name}: 0 PSMs"); return
    near=lambda d,t,w=0.05: abs(d-t)<w
    unmod=sum(1 for d in deltas if near(d,0.0))
    iso1 =sum(1 for d in deltas if near(d,1.00335) or near(d,-1.00335) or near(d,2.00671) or near(d,-2.00671))
    modded=n-unmod
    print(f"\n{name}: {n} PSMs @1%FDR | unmod {100*unmod/n:.1f}% | isotope-shift {100*iso1/n:.2f}% | modified {100*modded/n:.1f}%")
    # coarse 1-Da histogram over -50..+200 for the modified fraction
    hist=collections.Counter(round(d) for d in deltas if not near(d,0.0))
    top=sorted(hist.items(), key=lambda kv:-kv[1])[:12]
    print("  top delta-mass bins (Da: count):", ", ".join(f"{k:+d}:{v}" for k,v in top))
if __name__=="__main__":
    for spec in sys.argv[1:]:
        name,tsv=spec.split("=",1)
        try: summarize(name, load(tsv))
        except Exception as e: print(f"{name}: ERR {e}")
