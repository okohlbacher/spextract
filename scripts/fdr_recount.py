import csv, glob, os
# The benchmark counts peptides at spectrum_q<=0.01 (per-SPECTRUM FDR). That does NOT
# control peptide-level error: searching 2.4x more spectra gives each peptide more
# chances to pass. Recount at peptide_q to see which gains survive.
hdr = ("run", "PSMs", "spec_q pep", "pep_q pep", "delta")
print("%-22s %11s %11s %10s %8s" % hdr)
for d in sorted(glob.glob("/path/to/scratch/bench/out/*/results.sage.tsv")):
    n = os.path.basename(os.path.dirname(d))
    if not any(k in n for k in ("way2", "s23", "basepin", "mcw", "wide")): continue
    sp, pp, npsm = set(), set(), 0
    try:
        for x in csv.DictReader(open(d), delimiter="\t"):
            pep = x.get("peptide", "")
            if float(x.get("spectrum_q", 1)) <= 0.01: sp.add(pep); npsm += 1
            if float(x.get("peptide_q", 1)) <= 0.01: pp.add(pep)
    except Exception as e:
        print("%-22s ERROR %s" % (n, e)); continue
    dl = 100.0 * (len(pp) - len(sp)) / max(len(sp), 1)
    print("%-22s %11d %11d %10d %+7.1f%%" % (n, npsm, len(sp), len(pp), dl))
