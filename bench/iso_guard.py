#!/usr/bin/env python3
"""Can a guard available AT COLLAPSE TIME separate the free removals from the harmful ones?

Every removal is classed by what it cost (dark / covered / LOST, see iso_loss.py). The decision log
carries the features the extractor actually had when it decided: isotope step, the MS1 elution
correlation between removed and survivor, their intensities and isotope-peak counts. If the harmful
removals are lattice coincidences rather than duplicates, the elution correlation should separate
them -- two hypotheses of one peptide share a profile, an independent faint neighbour does not.

Usage: iso_guard.py <collapse.tsv> <base.tsv> <base.sage.tsv> <collapsed.sage.tsv>
"""
import sys, csv, collections


def ids(p):
    o = {}
    for x in csv.DictReader(open(p), delimiter="\t"):
        if x["label"] != "1" or float(x["spectrum_q"]) > 0.01: continue
        o[int(x["scannr"].split("=")[-1])] = x["peptide"]
    return o


def main(log, bt, bs, cs):
    bid, cid = ids(bs), ids(cs)
    lost = set(bid.values()) - set(cid.values())
    coord = {}
    with open(bt) as f:
        next(f)
        for line in f:
            q = line.split("\t")
            coord[(round(float(q[2]), 6), round(float(q[1]), 4), int(q[3]))] = int(q[0])

    rows = []
    with open(log) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            i = coord.get((round(float(x["rm_mz"]), 6), round(float(x["rm_rt"]), 4), int(x["rm_z"])))
            pep = bid.get(i) if i is not None else None
            cls = "unemitted" if i is None else ("dark" if pep is None else
                                                 ("LOST" if pep in lost else "covered"))
            si, ri = float(x["sv_int"]), float(x["rm_int"])
            rows.append({"cls": cls, "pep": pep, "k": abs(int(x["k"])), "corr": float(x["corr"]),
                         "ratio": ri / si if si > 0 else 0.0, "rm_niso": int(x["rm_niso"]),
                         "sv_niso": int(x["sv_niso"])})
    n = len(rows)
    c = collections.Counter(r["cls"] for r in rows)
    print(f"{n:,} logged removals: " + ", ".join(f"{k} {v:,}" for k, v in c.most_common()) + "\n")

    def table(key, bins, label):
        print(f"{label:<26}{'removals':>10}{'dark':>10}{'covered':>10}{'LOST':>8}{'harm%':>8}")
        for lo, hi in bins:
            sel = [r for r in rows if lo <= r[key] < hi]
            if not sel: continue
            cc = collections.Counter(r["cls"] for r in sel)
            harm = 100.0 * cc["LOST"] / len(sel)
            print(f"  [{lo:>6.2f}, {hi:>6.2f})       {len(sel):>10,}{cc['dark']:>10,}"
                  f"{cc['covered']:>10,}{cc['LOST']:>8,}{harm:>7.2f}%")
        print()

    table("corr", [(-2.01, -1.99), (-1.0, 0.0), (0.0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.85),
                   (0.85, 0.95), (0.95, 1.01)], "MS1 elution correlation")
    table("ratio", [(0.0, 0.02), (0.02, 0.1), (0.1, 0.3), (0.3, 0.7), (0.7, 1.5), (1.5, 1e9)],
          "intensity removed/survivor")

    print(f"{'guard (keep only these removals)':<44}{'removals':>10}{'dark kept':>11}"
          f"{'peptides at risk':>18}")
    base_pep = {r["pep"] for r in rows if r["cls"] == "LOST"}
    print(f"{'no guard (K<=3)':<44}{n:>10,}{c['dark']:>11,}{len(base_pep):>18,}")
    for name, fn in (
        ("k == 1", lambda r: r["k"] == 1),
        ("corr >= 0.7", lambda r: r["corr"] >= 0.7),
        ("corr >= 0.85", lambda r: r["corr"] >= 0.85),
        ("corr >= 0.95", lambda r: r["corr"] >= 0.95),
        ("k == 1 and corr >= 0.7", lambda r: r["k"] == 1 and r["corr"] >= 0.7),
        ("k == 1 and corr >= 0.85", lambda r: r["k"] == 1 and r["corr"] >= 0.85),
        ("k == 1 and corr >= 0.95", lambda r: r["k"] == 1 and r["corr"] >= 0.95),
        ("k == 1, corr >= 0.85, ratio >= 0.1",
         lambda r: r["k"] == 1 and r["corr"] >= 0.85 and r["ratio"] >= 0.1),
    ):
        sel = [r for r in rows if fn(r)]
        cc = collections.Counter(r["cls"] for r in sel)
        # a peptide is only truly at risk if EVERY removal that carried it survives the guard
        kept_lost = collections.Counter(r["pep"] for r in sel if r["cls"] == "LOST")
        all_lost = collections.Counter(r["pep"] for r in rows if r["cls"] == "LOST")
        risk = sum(1 for p in kept_lost if kept_lost[p] == all_lost[p])
        print(f"{name:<44}{len(sel):>10,}{cc['dark']:>11,}{risk:>18,}")


if __name__ == "__main__":
    if len(sys.argv) < 5: sys.exit(__doc__)
    main(*sys.argv[1:5])
