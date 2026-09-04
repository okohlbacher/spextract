#!/usr/bin/env python3
"""Are isotope-linked pseudo-spectra "basically identical" in CONTENT, or only in precursor lattice?

The emission-level test (iso_dup.py) shows ~28% of our spectra sit one to three isotope steps
above another co-eluting spectrum of the same charge. That alone does not say they are the same
spectrum: sharing a precursor lattice position is not sharing fragments.

  HIGH cosine -> same content under different mass labels. The duplication is a MASS-HYPOTHESIS
                 hedge; one spectrum plus a wider isotope_errors would do, and merging is safe.
  LOW cosine  -> different fragment sets that merely sit on the isotope lattice. Merging destroys
                 content, which is exactly why every previous merge variant lost peptides.

Controls (a cosine without a null is meaningless):
  lattice pair : the isotope-linked partner                     (the thing under test)
  neighbour    : nearest co-eluting spectrum, same z/IM, NOT on the isotope lattice
  random       : two unrelated spectra                          (the floor)

Usage: iso_sim.py <iso.tsv> <pseudo.mzML> [--n 2000] [--rt 3.0] [--im 0.02] [--ppm 20]
"""
import sys, random, collections, math, base64, struct, re

ISO, BIN = 1.00335, 0.02


def cosine(a, b):
    """Binned normalised dot product, sqrt-intensity (standard for spectral similarity)."""
    va, vb = collections.defaultdict(float), collections.defaultdict(float)
    for mz, it in a: va[round(mz / BIN)] += math.sqrt(it)
    for mz, it in b: vb[round(mz / BIN)] += math.sqrt(it)
    na = math.sqrt(sum(v * v for v in va.values())); nb = math.sqrt(sum(v * v for v in vb.values()))
    if na == 0 or nb == 0: return 0.0
    return sum(v * vb[k] for k, v in va.items() if k in vb) / (na * nb)


def fetch(mzmlp, want):
    """Pull the peak arrays for the wanted spectrum indices. The mzML is uncompressed base64,
    64-bit m/z and 32-bit intensity, one tag per line -- a plain stream is simpler than a library."""
    out, idx, keep, arr, prec, cur = {}, -1, False, None, 8, {}
    with open(mzmlp, "r", errors="replace") as f:
        for line in f:
            if "<spectrum " in line:
                m = re.search(r'index="(\d+)"', line)
                idx = int(m.group(1)) if m else idx + 1
                keep, cur = idx in want, {}
                continue
            if not keep: continue
            if 'name="m/z array"' in line: arr = "mz"
            elif 'name="intensity array"' in line: arr = "it"
            elif 'name="64-bit float"' in line: prec = 8
            elif 'name="32-bit float"' in line: prec = 4
            elif "<binary>" in line:
                b = line[line.find("<binary>") + 8: line.find("</binary>")]
                if b and arr:
                    raw = base64.b64decode(b)
                    cur[arr] = struct.unpack(f"<{len(raw)//prec}{'d' if prec == 8 else 'f'}", raw)
                arr = None
            elif "</spectrum>" in line:
                if "mz" in cur and "it" in cur: out[idx] = list(zip(cur["mz"], cur["it"]))
                want.discard(idx); keep = False
                if not want: break
    return out


def main(tsv, mzmlp, N=2000, rt_tol=3.0, im_tol=0.02, ppm=20.0):
    rows = []
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append((float(p[1]), float(p[2]), int(p[3]), float(p[4])))
    n = len(rows)
    buckets = collections.defaultdict(list)
    for i, (rt, mz, z, im) in enumerate(rows):
        buckets[(z, int(im / im_tol), int(rt / rt_tol))].append((mz, i))

    def probe(i, offset):
        rt, mz, z, im = rows[i]
        if z <= 0: return -1
        tol, want = mz * ppm * 1e-6, mz - offset
        ib, rb = int(im / im_tol), int(rt / rt_tol)
        for di in (-1, 0, 1):
            for dr in (-1, 0, 1):
                for (pmz, j) in buckets.get((z, ib + di, rb + dr), ()):
                    if j != i and abs(pmz - want) <= tol and abs(rows[j][0] - rt) <= rt_tol \
                       and abs(rows[j][3] - im) <= im_tol:
                        return j
        return -1

    random.seed(7)
    order = list(range(n)); random.shuffle(order)
    lattice, neigh = [], []
    for i in order:
        if len(lattice) < N:
            j = probe(i, ISO / max(rows[i][2], 1))
            if j >= 0: lattice.append((i, j))
        if len(neigh) < N:
            j = probe(i, 0.5 / max(rows[i][2], 1))          # same size step, off the lattice
            if j >= 0: neigh.append((i, j))
        if len(lattice) >= N and len(neigh) >= N: break
    rnd = [(random.randrange(n), random.randrange(n)) for _ in range(N)]

    need = sorted({x for p in (lattice, neigh, rnd) for pr in p for x in pr})
    print(f"pairs: lattice={len(lattice)} neighbour={len(neigh)} random={len(rnd)}; "
          f"fetching {len(need):,} spectra")
    peaks = fetch(mzmlp, set(need))
    print(f"fetched {len(peaks):,}")

    for name, prs in (("lattice (isotope)", lattice), ("neighbour (off-lattice)", neigh),
                      ("random", rnd)):
        cs = sorted(cosine(peaks[i], peaks[j]) for i, j in prs if i in peaks and j in peaks)
        if not cs: continue
        q = lambda f: cs[min(len(cs) - 1, int(f * len(cs)))]
        print(f"{name:<26} n={len(cs):>5}  median={q(.5):.3f}  p25={q(.25):.3f}  p75={q(.75):.3f}"
              f"  >0.9: {100.0*sum(1 for c in cs if c > 0.9)/len(cs):.1f}%")


if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    kw = {}
    for a in sys.argv[3:]:
        k, v = a.lstrip("-").split("=")
        kw[{"n": "N", "rt": "rt_tol", "im": "im_tol", "ppm": "ppm"}[k]] = int(v) if k == "n" else float(v)
    main(sys.argv[1], sys.argv[2], **kw)
