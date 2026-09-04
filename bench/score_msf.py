#!/usr/bin/env python3
"""Score MSFragger .tsv arms under a common target-decoy FDR + entrapment, for the cross-engine
confirmation of the min_correlation tightening. MSFragger hyperscore is RAW (no q-value), so FDR
is applied here identically across arms -- the same procedure joint_bench.py uses, so Sage and
MSFragger are compared on one FDR machine, not each engine's own.

Peptide-level FDR: best hyperscore per distinct peptide, then target-decoy walk. Entrapment: among
accepted TARGET peptides, the fraction whose proteins are all ENTRAP_ (a known false positive).
"""
import csv, sys, glob, math

TAG = "ENTRAP_"


def load(tsv):
    rows = []
    with open(tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames or []
        def pick(*c):
            for x in c:
                if x in cols: return x
            return None
        cp = pick("peptide", "Peptide", "Modified Peptide", "modified_peptide")
        cs = pick("hyperscore", "Hyperscore")
        cpr = pick("proteins", "Proteins", "protein", "Protein")
        cprot_alt = pick("alternative_proteins", "Alternative Proteins")
        for x in rd:
            try:
                s = float(x[cs]); pep = x[cp]
            except (TypeError, ValueError, KeyError):
                continue
            prot = (x.get(cpr) or "") + ";" + (x.get(cprot_alt) or "")
            plist = [p for p in prot.replace(",", ";").split(";") if p.strip()]
            dec = bool(plist) and all(p.startswith("rev_") for p in plist)
            ent = bool(plist) and all(p.startswith(TAG) or p.startswith("rev_" + TAG) for p in plist)
            rows.append((pep, s, dec, ent))
    return rows


def score(rows, ratio, fdr=0.01):
    best = {}
    for pep, s, dec, ent in rows:
        k = (pep, dec)
        if k not in best or s > best[k][0]:
            best[k] = (s, dec, ent)
    ranked = sorted(best.values(), key=lambda r: -r[0])
    nt = nd = 0; acc = []
    thr = None
    for s, dec, ent in ranked:
        if dec: nd += 1
        else: nt += 1; acc.append(ent)
        if nt and nd / nt <= fdr:
            thr = len(acc)
    acc = acc[:thr] if thr else []
    ne = sum(1 for e in acc if e)
    ntg = len(acc) - ne
    raw = 100.0 * ne / max(len(acc), 1)
    ent_fdr = 100.0 * (ne / ratio) * 2.0 / max(len(acc), 1) if ratio > 0 else float("nan")
    return len(acc), ntg, ne, raw, ent_fdr


ratio = float(sys.argv[1])          # entrap/target protein ratio (0.800)
print("%-10s %9s %9s %8s %11s %11s" % ("arm", "peptides", "target", "entrap", "raw entrap%", "entrap FDR%"))
for spec in sys.argv[2:]:
    name, d = spec.split("=", 1)
    t = glob.glob(d + "/*.tsv")
    if not t:
        print("%-10s (no tsv)" % name); continue
    n, ntg, ne, raw, ef = score(load(t[0]), ratio)
    print("%-10s %9d %9d %8d %10.2f%% %10.2f%%" % (name, n, ntg, ne, raw, ef))
