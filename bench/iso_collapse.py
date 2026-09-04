#!/usr/bin/env python3
"""Given the isotope-linked groups found by iso_dup.py, can ANY id-agnostic rule collapse them
without losing peptides?

Each group is a connected component of the "one isotope step apart, same charge, co-eluting,
co-mobile" graph. Collapsing keeps ONE member per group. The rules differ only in which member:

  lightest   lowest precursor m/z              -- the naive "the mono is the lightest" rule
  niso       most isotope partners (spx_n_isotopes), ties -> lightest
  richest    most fragment peaks               -- content, not mass, decides
  oracle     the member Sage identified        -- NOT implementable, the ceiling

Reported against the do-nothing control: spectra kept, peptides still covered.
Usage: iso_collapse.py <iso.tsv> <results.sage.tsv> [--rt 3.0] [--im 0.02] [--ppm 20]
"""
import sys, csv, collections

ISO = 1.00335


def main(tsv, sage, rt_tol=3.0, im_tol=0.02, ppm=20.0):
    rows = []
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append((float(p[1]), float(p[2]), int(p[3]), float(p[4]), int(p[5]),
                         int(p[7] or 0)))
    n = len(rows)

    buckets = collections.defaultdict(list)
    for i, (rt, mz, z, im, npk, ni) in enumerate(rows):
        buckets[(z, int(im / im_tol), int(rt / rt_tol))].append((mz, i))

    # union-find over isotope links
    par = list(range(n))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb

    links = 0
    for i, (rt, mz, z, im, npk, ni) in enumerate(rows):
        if z <= 0: continue
        tol = mz * ppm * 1e-6
        ib, rb = int(im / im_tol), int(rt / rt_tol)
        for k in (1, 2, 3):
            want = mz - k * ISO / z
            for di in (-1, 0, 1):
                for dr in (-1, 0, 1):
                    for (pmz, j) in buckets.get((z, ib + di, rb + dr), ()):
                        if j == i or abs(pmz - want) > tol: continue
                        if abs(rows[j][0] - rt) <= rt_tol and abs(rows[j][3] - im) <= im_tol:
                            union(i, j); links += 1

    groups = collections.defaultdict(list)
    for i in range(n): groups[find(i)].append(i)
    multi = {g: m for g, m in groups.items() if len(m) > 1}
    print(f"{n:,} spectra, {links:,} isotope links, {len(groups):,} groups "
          f"({len(multi):,} with >1 member, covering {sum(len(m) for m in multi.values()):,} spectra)")
    sz = collections.Counter(len(m) for m in multi.values())
    print("group sizes:", dict(sorted(sz.items())[:8]))

    idp = {}
    with open(sage) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            if x["label"] != "1" or float(x["spectrum_q"]) > 0.01: continue
            idp[int(x["scannr"].split("=")[-1])] = x["peptide"]
    allpep = set(idp.values())
    print(f"\ncontrol: {n:,} spectra, {len(allpep):,} peptides at 1% spectrum-FDR\n")

    rules = {
        "lightest": lambda m: min(m, key=lambda i: rows[i][1]),
        "niso":     lambda m: max(m, key=lambda i: (rows[i][5], -rows[i][1])),
        "richest":  lambda m: max(m, key=lambda i: (rows[i][4], -rows[i][1])),
        "oracle":   lambda m: next((i for i in m if i in idp), m[0]),
    }
    print(f"{'rule':<10}{'spectra kept':>14}{'vs control':>12}{'peptides kept':>15}{'lost':>8}")
    for name, pick in rules.items():
        keep = set()
        for g, m in groups.items():
            keep.add(pick(m) if len(m) > 1 else m[0])
        pep = {idp[i] for i in keep if i in idp}
        print(f"{name:<10}{len(keep):>14,}{100.0*len(keep)/n:>11.1f}%{len(pep):>15,}"
              f"{len(allpep)-len(pep):>8,}")


if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    kw = {}
    for a in sys.argv[3:]:
        k, v = a.lstrip("-").split("=")
        kw[{"rt": "rt_tol", "im": "im_tol", "ppm": "ppm"}[k]] = float(v)
    main(sys.argv[1], sys.argv[2], **kw)
