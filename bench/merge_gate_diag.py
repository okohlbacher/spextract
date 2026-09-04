#!/usr/bin/env python3
"""Should two adjacent spectra be merged? Compare precursor deltas for SAME vs DIFFERENT peptide.

The merge gate is currently (RT window, 20 ppm m/z, 0.02 1/K0, same charge) -- purely precursor
coordinates. Merging loses peptides monotonically, and PSMs per SPECTRUM falls 8% (0.432 ->
0.398), so merged spectra identify LESS often than unmerged ones. That is the signature of
fusing unrelated precursors into chimeras, not of recovering split halves.

This measures the gate's discriminating power directly. For every pair of spectra adjacent in
RT that BOTH identified, split by whether they got the SAME peptide:

  SAME peptide      -> these SHOULD be merged (same feature, split across emissions)
  DIFFERENT peptide -> these MUST NOT be merged (merging fuses two peptides into a chimera)

and plot the distribution of |d m/z| (ppm), |d IM|, |d RT| for each class. Where the two
distributions overlap, no coordinate threshold can separate them -- that overlap IS the
chimera rate the current gate cannot avoid.

Deliberately restricted to pairs where BOTH spectra identified: those are the only pairs where
"same or different peptide" is known. This biases toward identifiable spectra, which is stated
rather than hidden -- but the bias runs AGAINST the gate looking bad, since identifiable spectra
are the cleaner population.
"""
import argparse, collections, csv, json, math, sys

PROTON = 1.00727646


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", help="Sage results.sage.tsv")
    ap.add_argument("--rt-window", type=float, default=2.6,
                    help="seconds; pairs closer than this in RT are merge candidates")
    ap.add_argument("--json", help="write the distributions here for plotting")
    a = ap.parse_args()

    rows = []
    with open(a.tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            if any(p.startswith("rev_") for p in x["proteins"].split(";") if p):
                continue
            try:
                if float(x.get("peptide_q", 1)) > 0.01:
                    continue
                z = int(x["charge"])
                rows.append({"pep": x["peptide"], "rt": float(x["rt"]) * 60.0,
                             "mz": (float(x["expmass"]) + z * PROTON) / z,
                             "im": float(x.get("ion_mobility") or 0), "z": z})
            except (TypeError, ValueError, KeyError):
                continue
    if not rows:
        sys.exit("ABORT: no identified PSMs parsed from %s" % a.tsv)
    rows.sort(key=lambda r: r["rt"])
    print("[data] %d identified PSMs at peptide_q<=0.01" % len(rows), flush=True)

    same, diff = collections.defaultdict(list), collections.defaultdict(list)
    n_pairs = 0
    for i, r in enumerate(rows):
        for j in range(i + 1, len(rows)):
            s = rows[j]
            if s["rt"] - r["rt"] > a.rt_window:
                break
            # the current gate also requires equal charge; keep that so this measures the gate
            # as it actually is, not an idealised version
            if s["z"] != r["z"]:
                continue
            d = same if s["pep"] == r["pep"] else diff
            d["ppm"].append(abs(s["mz"] - r["mz"]) / min(s["mz"], r["mz"]) * 1e6)
            if r["im"] > 0 and s["im"] > 0:
                d["im"].append(abs(s["im"] - r["im"]))
            d["rt"].append(abs(s["rt"] - r["rt"]))
            n_pairs += 1

    def q(v, p):
        v = sorted(v)
        return v[int(p * (len(v) - 1))] if v else float("nan")

    print("[pairs] %d adjacent same-charge pairs within %.1f s: %d SAME peptide, %d DIFFERENT"
          % (n_pairs, a.rt_window, len(same["ppm"]), len(diff["ppm"])))
    print()
    print("%-22s %10s %10s %10s %10s %10s" % ("", "n", "median", "p75", "p90", "p99"))
    print("-" * 76)
    for key, unit in (("ppm", "ppm"), ("im", "1/K0"), ("rt", "s")):
        for name, d in (("SAME peptide", same), ("DIFF peptide", diff)):
            v = d[key]
            if not v:
                continue
            print("%-14s %-7s %10d %10.4f %10.4f %10.4f %10.4f"
                  % (name, unit, len(v), q(v, .5), q(v, .75), q(v, .90), q(v, .99)))
        print("-" * 76)

    # The number that matters: at the CURRENT gate, what fraction of admitted pairs are
    # DIFFERENT peptides? Those are the chimeras the gate cannot refuse.
    print()
    print("chimera rate admitted by a coordinate gate (pairs passing BOTH m/z and IM cuts):")
    print("%-12s %-12s %10s %10s %12s" % ("m/z ppm", "d IM", "SAME", "DIFF", "% chimera"))
    out = []
    for ppm_cut in (5.0, 10.0, 20.0):
        for im_cut in (0.005, 0.01, 0.02):
            ns = sum(1 for k in range(len(same["ppm"]))
                     if same["ppm"][k] <= ppm_cut and k < len(same["im"]) and same["im"][k] <= im_cut)
            nd = sum(1 for k in range(len(diff["ppm"]))
                     if diff["ppm"][k] <= ppm_cut and k < len(diff["im"]) and diff["im"][k] <= im_cut)
            frac = 100.0 * nd / max(ns + nd, 1)
            out.append({"ppm": ppm_cut, "im": im_cut, "same": ns, "diff": nd, "chimera_pct": frac})
            print("%-12.1f %-12.3f %10d %10d %11.1f%%" % (ppm_cut, im_cut, ns, nd, frac))
    print()
    print("A gate can only work if tightening it drops DIFF much faster than SAME.")
    print("If %chimera stays flat as the cuts tighten, precursor coordinates cannot separate")
    print("the two populations and merging on them will always fuse peptides.")

    if a.json:
        json.dump({"same": {k: v for k, v in same.items()},
                   "diff": {k: v for k, v in diff.items()},
                   "gate_scan": out}, open(a.json, "w"))
        print("\nwrote %s" % a.json)


if __name__ == "__main__":
    sys.exit(main())
