#!/usr/bin/env python3
"""Are within-cycle repeats SPLIT SPECTRA (one feature cut in two) or DISTINCT features?

Two competing explanations for 2.85 PSMs per (peptide, cycle):

  (a) MISANNOTATION -- distinct features wrongly given the same precursor/charge, so they land in
      the same group without being the same thing.
  (b) SPLITNESS     -- ONE true spectrum whose peaks were divided across two or more emitted
      spectra. The opposite of a chimera.

They make opposite predictions, and the data can separate them:

                        (a) distinct features        (b) split spectrum
  d(IM), d(RT)          SEPARATED                    CLUSTERED at ~0
  peak overlap          moderate, incidental         LOW (peaks partitioned)
  union vs parts        union = mixture, no better   union RICHER than either part
  cosine                intermediate                 intermediate  <- cannot distinguish alone

The measured cosine of 0.845 was read earlier as "partially distinct spectra". Splitness predicts
the SAME intermediate cosine -- two halves of one spectrum share some peaks and each carry unique
ones. Cosine alone cannot tell them apart; (IM, RT) geometry and peak complementarity can.

Controls, since a spread without a null is meaningless -- the MS1 funnel was retracted for that:
  * same peptide, DIFFERENT cycle -> genuine elution, must be SEPARATED in RT by construction
  * the tolerance used to build a group bounds how close members can be, so the null is
    reported alongside every spread.
"""
import argparse, collections, csv, math, random, sys
from pyteomics import mzml

CYCLE_S = 1.385
BIN = 0.02


def peaks(s):
    return {int(m / BIN) for m in s["m/z array"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv")
    ap.add_argument("mzml")
    ap.add_argument("--label", default="run")
    ap.add_argument("--max-groups", type=int, default=3000)
    a = ap.parse_args()
    random.seed(5)

    rows = []
    with open(a.tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                rows.append({"pep": x["peptide"], "sid": x["scannr"],
                             "rt": float(x["rt"]) * 60.0,
                             "im": float(x.get("ion_mobility") or 0),
                             "mz": float(x["expmass"]), "z": int(x["charge"]),
                             "score": float(x["sage_discriminant_score"]),
                             "decoy": all(p.startswith("rev_")
                                          for p in x["proteins"].split(";") if p)})
            except (TypeError, ValueError, KeyError):
                continue
    ranked = sorted(rows, key=lambda r: -r["score"])
    nt = nd = 0; keep, acc = [], []
    for r in ranked:
        if r["decoy"]: nd += 1
        else: nt += 1; acc.append(r)
        if nt and nd / nt <= 0.01: keep = list(acc)

    g = collections.defaultdict(list)
    for r in keep:
        g[(r["pep"], int(r["rt"] / CYCLE_S))].append(r)
    dup = {k: v for k, v in g.items() if len(v) > 1}
    sel = dict(random.sample(sorted(dup.items()), min(a.max_groups, len(dup))))

    # control: same peptide, different cycle
    bypep = collections.defaultdict(list)
    for (pep, cyc), v in g.items():
        bypep[pep].append((cyc, v[0]))
    ctrl = []
    for pep, e in bypep.items():
        if len(e) > 1:
            (c1, r1), (c2, r2) = random.sample(e, 2)
            if c1 != c2:
                ctrl.append((r1, r2))
    ctrl = ctrl[:len(sel)]

    print("[%s] %d duplicated groups, sampling %d" % (a.label, len(dup), len(sel)), flush=True)

    need = set()
    for v in sel.values():
        need.update(r["sid"] for r in v[:3])
    for r1, r2 in ctrl:
        need.add(r1["sid"]); need.add(r2["sid"])
    P = {}
    for s in mzml.read(a.mzml):
        if s.get("id") in need:
            P[s["id"]] = peaks(s)
            if len(P) == len(need): break

    d_im, d_rt, d_mz, same_z = [], [], [], 0
    jac, uniq_frac, union_gain = [], [], []
    for v in sel.values():
        for i in range(min(len(v), 3)):
            for j in range(i + 1, min(len(v), 3)):
                A, B = v[i], v[j]
                d_rt.append(abs(A["rt"] - B["rt"]))
                if A["im"] > 0 and B["im"] > 0:
                    d_im.append(abs(A["im"] - B["im"]))
                d_mz.append(abs(A["mz"] - B["mz"]) / min(A["mz"], B["mz"]) * 1e6)
                same_z += (A["z"] == B["z"])
                pa, pb = P.get(A["sid"]), P.get(B["sid"])
                if pa and pb:
                    inter, union = len(pa & pb), len(pa | pb)
                    jac.append(inter / union if union else 0)
                    uniq_frac.append(1 - inter / min(len(pa), len(pb)))
                    union_gain.append(union / max(len(pa), len(pb)))

    c_im = [abs(x["im"] - y["im"]) for x, y in ctrl if x["im"] > 0 and y["im"] > 0]
    c_rt = [abs(x["rt"] - y["rt"]) for x, y in ctrl]

    def q(v, p):
        v = sorted(v); return v[int(p * (len(v) - 1))] if v else float("nan")

    print()
    print("%-38s %8s %8s %8s" % ("", "median", "p25", "p75"))
    print("-" * 66)
    print("%-38s %8.4f %8.4f %8.4f" % ("d(IM) within-cycle pairs, 1/K0", q(d_im,.5), q(d_im,.25), q(d_im,.75)))
    print("%-38s %8.4f %8.4f %8.4f" % ("  CONTROL: same pep, diff cycle", q(c_im,.5), q(c_im,.25), q(c_im,.75)))
    print("%-38s %8.3f %8.3f %8.3f" % ("d(RT) within-cycle pairs, s", q(d_rt,.5), q(d_rt,.25), q(d_rt,.75)))
    print("%-38s %8.1f %8.1f %8.1f" % ("  CONTROL: same pep, diff cycle", q(c_rt,.5), q(c_rt,.25), q(c_rt,.75)))
    print("%-38s %8.1f %8.1f %8.1f" % ("d(precursor m/z), ppm", q(d_mz,.5), q(d_mz,.25), q(d_mz,.75)))
    print("-" * 66)
    print("%-38s %8.3f %8.3f %8.3f" % ("Jaccard peak overlap", q(jac,.5), q(jac,.25), q(jac,.75)))
    print("%-38s %8.3f %8.3f %8.3f" % ("fraction of peaks UNIQUE to one", q(uniq_frac,.5), q(uniq_frac,.25), q(uniq_frac,.75)))
    print("%-38s %8.3f %8.3f %8.3f" % ("union size / larger member", q(union_gain,.5), q(union_gain,.25), q(union_gain,.75)))
    print("-" * 66)
    print("identical charge in pair: %.1f%%" % (100.0 * same_z / max(len(d_rt), 1)))
    print()
    print("SPLIT  predicts: d(IM),d(RT) ~ 0 (clustered); LOW Jaccard; HIGH unique fraction;")
    print("                 union/larger >> 1 (the parts are complementary)")
    print("DISTINCT predicts: d(IM),d(RT) separated; overlap incidental; union/larger ~ 1")


if __name__ == "__main__":
    sys.exit(main())
