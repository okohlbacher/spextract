#!/usr/bin/env python3
"""Is our emission inflated by the SAME peptide being emitted once per isotope peak?

Hypothesis (user, 2026-09-03): we select M+1/M+2 as a precursor in its own right, so one
peptide yields several near-identical pseudo-spectra that differ by ~1.00335/z in reported
precursor m/z. If true, collapsing isotope-linked precursors would cut emission with no loss.

This measures the CEILING directly on an emitted mzML -- every spectrum, not just the
identified ones -- and against a DECOY OFFSET, because at ~10^6 spectra in a narrow m/z range
isotope-spacing coincidences happen by chance and an uncontrolled count would be meaningless.

  isotope offset   dm = k * 1.00335 / z     k = 1..3
  decoy offset     dm = k * 1.00335 / z + 0.5 / z    (same magnitude, not a real spacing)

A spectrum is a REMOVABLE duplicate if another emitted spectrum sits one (or two, or three)
isotope steps BELOW it at the same charge, co-eluting in RT and co-located in IM.

Usage:  iso_dup.py extract <pseudo.mzML> <out.tsv>
        iso_dup.py analyse <out.tsv> [--rt 3.0] [--im 0.02] [--ppm 20]
"""
import sys, re, collections

ISO = 1.00335


def extract(mzml, out):
    """Stream the mzML and write one row per spectrum. Header fields are plain XML."""
    pat = {
        "rt":  re.compile(r'name="scan start time" value="([^"]+)"'),
        "mz":  re.compile(r'name="selected ion m/z" value="([^"]+)"'),
        "z":   re.compile(r'name="charge state" value="([^"]+)"'),
        "im":  re.compile(r'name="inverse reduced ion mobility" value="([^"]+)"'),
    }
    cur, n = {}, 0
    with open(mzml, "r", errors="replace") as f, open(out, "w") as w:
        w.write("index\trt\tmz\tz\tim\tnpeaks\tguessed\tniso\n")
        for line in f:
            if "<spectrum " in line:
                cur = {}
                m = re.search(r'index="(\d+)"', line)
                if m: cur["index"] = m.group(1)
                m = re.search(r'defaultArrayLength="(\d+)"', line)
                cur["npeaks"] = m.group(1) if m else "0"
                continue
            if not cur:
                continue
            if "spx_guessed" in line:
                cur["guessed"] = re.search(r'value="(\d+)"', line).group(1); continue
            if "spx_n_isotopes" in line:
                cur["niso"] = re.search(r'value="(\d+)"', line).group(1); continue
            for k, p in pat.items():
                if k not in cur:
                    m = p.search(line)
                    if m: cur[k] = m.group(1)
            if "</precursorList>" in line:
                if {"index", "rt", "mz", "z", "im"} <= cur.keys():
                    w.write("\t".join((cur["index"], cur["rt"], cur["mz"], cur["z"], cur["im"],
                                       cur.get("npeaks", "0"), cur.get("guessed", ""), cur.get("niso", ""))) + "\n")
                    n += 1
                cur = {}
    print(f"wrote {n:,} spectra to {out}")


def analyse(tsv, rt_tol=3.0, im_tol=0.02, ppm=20.0, sage=None):
    rows = []
    with open(tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append((float(p[1]), float(p[2]), int(p[3]), float(p[4]), int(p[5]),
                         p[6], p[7]))
    n = len(rows)
    print(f"{n:,} emitted spectra\n")

    # bucket on (charge, IM cell, RT cell); a partner may sit in a neighbouring cell, so probe +-1
    buckets = collections.defaultdict(list)
    for i, (rt, mz, z, im, npk, g, ni) in enumerate(rows):
        buckets[(z, int(im / im_tol), int(rt / rt_tol))].append((mz, i))

    def scan(decoy):
        """Mark every spectrum that has a partner k isotope steps BELOW it. Returns marks, khist, partner."""
        mark = bytearray(n)
        partner = [-1] * n
        khist = collections.Counter()
        for i, (rt, mz, z, im, npk, g, ni) in enumerate(rows):
            if z <= 0: continue
            tol = mz * ppm * 1e-6
            ib, rb = int(im / im_tol), int(rt / rt_tol)
            hit = 0
            for k in (1, 2, 3):
                want = mz - k * ISO / z - (0.5 / z if decoy else 0.0)
                for di in (-1, 0, 1):
                    for dr in (-1, 0, 1):
                        for (pmz, j) in buckets.get((z, ib + di, rb + dr), ()):
                            if j == i or abs(pmz - want) > tol: continue
                            prt, pim = rows[j][0], rows[j][3]
                            if abs(prt - rt) <= rt_tol and abs(pim - im) <= im_tol:
                                hit = k; partner[i] = j; break
                        if hit: break
                    if hit: break
                if hit: break
            if hit:
                mark[i] = 1; khist[hit] += 1
        return mark, khist, partner

    iso_mark, iso_k, iso_partner = scan(False)
    dec_mark, dec_k, _ = scan(True)
    ni, nd = sum(iso_mark), sum(dec_mark)
    excess = ni - nd
    print(f"{'offset':<18}{'spectra with a partner below':>30}{'% of emission':>16}")
    print(f"{'isotope k*1.00335/z':<18}{ni:>30,}{100.0*ni/n:>15.1f}%")
    print(f"{'decoy  +0.5/z':<18}{nd:>30,}{100.0*nd/n:>15.1f}%")
    print(f"{'EXCESS (real)':<18}{excess:>30,}{100.0*excess/n:>15.1f}%\n")
    print("isotope hits by step k:", dict(iso_k), " decoy:", dict(dec_k))
    print(f"\nceiling: collapsing every isotope-linked duplicate leaves {n-excess:,} spectra "
          f"({100.0*(n-excess)/n:.1f}% of current emission)")

    # who are they? a duplicate that is itself a rich spectrum is a different problem from a thin one
    for name, mk in (("isotope-dup", iso_mark), ("all", bytearray(b"\x01" * n))):
        sel = [rows[i] for i in range(n) if mk[i]]
        if not sel: continue
        npk = sorted(r[4] for r in sel)
        gz = sum(1 for r in sel if r[5] == "1")
        print(f"{name:<12} n={len(sel):>9,}  median peaks={npk[len(npk)//2]:>4}  "
              f"guessed={100.0*gz/len(sel):>5.1f}%")

    if sage:
        # Does the HEAVY member carry an identification the LIGHT one does not? That is the only
        # way collapsing can lose a peptide, and it is what killed every previous merge variant.
        import csv as _csv
        idp, off = {}, {}
        with open(sage) as f:
            for x in _csv.DictReader(f, delimiter="\t"):
                if x["label"] != "1" or float(x["spectrum_q"]) > 0.01: continue
                i = int(x["scannr"].split("=")[-1])
                idp[i] = x["peptide"]
                # Sage's own isotope correction: how far off was the mass WE reported?
                off[i] = round((float(x["expmass"]) - float(x["calcmass"])) / ISO)
        print(f"\nSage 1% spectrum-FDR identifies {len(idp):,} of {n:,} spectra "
              f"({100.0*len(idp)/n:.1f}%)")
        c = collections.Counter()
        lost = set()
        for i in range(n):
            j = iso_partner[i]
            if j < 0: continue
            a, b = idp.get(i), idp.get(j)
            if a is None and b is None: c["neither identified"] += 1
            elif a is None: c["only the LIGHT one identified"] += 1
            elif b is None:
                c["only the HEAVY one identified (collapse risk)"] += 1; lost.add(a)
            elif a == b: c["both, SAME peptide (redundant)"] += 1
            else: c["both, DIFFERENT peptides"] += 1
        tot = sum(c.values())
        print(f"\nisotope-linked pairs, identification of heavy vs light ({tot:,} pairs):")
        for k2, v in c.most_common():
            print(f"  {k2:<46}{v:>9,}{100.0*v/tot:>8.1f}%")
        # Which member did Sage think was the real monoisotope? offset 0 = the mass we reported
        # was right; +1/+2 = we reported an M+k peak; -1 = we reported a peak BELOW the true mono.
        oh, ol = collections.Counter(), collections.Counter()
        for i in range(n):
            j = iso_partner[i]
            if j < 0: continue
            if i in off and j not in off: oh[off[i]] += 1
            if j in off and i not in off: ol[off[j]] += 1
        print("\nSage isotope offset of the identified member (0 = the mass we reported was right):")
        print("  heavy-only identified:", dict(sorted(oh.items())))
        print("  light-only identified:", dict(sorted(ol.items())))

        allp = set(idp.values())
        keep = {idp[i] for i in range(n) if idp.get(i) and not iso_mark[i]}
        print(f"\npeptides seen anywhere: {len(allp):,};  still seen after dropping every heavy "
              f"member: {len(keep):,}  (would lose {len(allp-keep):,})")


if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    if sys.argv[1] == "extract": extract(sys.argv[2], sys.argv[3])
    else:
        kw = {}
        for a in sys.argv[3:]:
            k, v = a.lstrip("-").split("="); kw["sage"] = v if k == "sage" else kw.get("sage")
            if k != "sage": kw[{"rt": "rt_tol", "im": "im_tol", "ppm": "ppm"}[k]] = float(v)
        analyse(sys.argv[2], **kw)
