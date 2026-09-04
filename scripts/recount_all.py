import csv, glob, os
# Every lever was judged at spectrum_q. Levers that REDUCE emitted spectra (competitive,
# rank-pruning, min_corr, dedup) look like pure losses there but should gain at peptide_q,
# where redundancy is penalised rather than rewarded. Re-judge them all.
rows = []
for d in glob.glob("/path/to/scratch/bench/**/results.sage.tsv", recursive=True):
    n = "/".join(d.split("/")[-3:-1])
    sp, pp = set(), set()
    try:
        for x in csv.DictReader(open(d), delimiter="\t"):
            pep = x.get("peptide", "")
            if float(x.get("spectrum_q", 1)) <= 0.01: sp.add(pep)
            if float(x.get("peptide_q", 1)) <= 0.01: pp.add(pep)
    except Exception:
        continue
    if not sp: continue
    rows.append((len(pp), len(sp), 100.0*(len(pp)-len(sp))/len(sp), n))
rows.sort(reverse=True)
print("%-38s %10s %10s %8s" % ("run", "spec_q", "PEP_Q", "FDR loss"))
for pp, sp, dl, n in rows:
    print("%-38s %10d %10d %+7.1f%%" % (n, sp, pp, dl))
