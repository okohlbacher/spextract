#!/usr/bin/env python3
"""FWHM of MS1 features as a function of retention time.

Why this matters for the across-cycle merge
-------------------------------------------
The across-cycle merge was killed on the grounds that "every co-elution statistic decays
smoothly to its null with no breakpoint, so no rt_window can be read off the data". That
conclusion assumes ONE global window. If chromatographic peak width varies with RT, a single
window is wrong by construction and a pooled analysis smears exactly the breakpoint it looks for.

A first pass using DIA-NN's RT.Stop - RT.Start said width is flat (slope -0.004 s/min), but that
proxy is quantised to 1.385 s -- one cycle -- so it cannot resolve a change below ~14%. This uses
FWHM measured per MS1 trace from its own XIC.

Two things come out of it:
  1. Does width vary with RT?  -> is a single rt_window defensible at all?
  2. What IS the width, in cycles?  -> the principled bound the review said the data could not give.

Context: we currently spread a peptide over 2.66 cycles. The coarse proxy put the elution window
at ~7.5 cycles, i.e. we may be UNDER-covering the peak rather than over-spreading it.

Robustness: medians and quantiles throughout, never means -- FWHM distributions are heavy-tailed
and a mean is dominated by a handful of merged or truncated traces. Traces whose XIC never falls
to half maximum on both sides are reported as -1 by the extractor and EXCLUDED here: a truncated
peak is not a narrow one, and imputing it would bias the curve downward exactly where peaks are
broadest.
"""
import argparse, csv, json, math, sys

CYCLE_S = 1.385


def theil_sen(xs, ys, cap=400):
    """Median of pairwise slopes -- robust to the heavy tail, no dependency on numpy/scipy."""
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan")
    step = max(1, n // cap)
    idx = list(range(0, n, step))
    sl = []
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            a, b = idx[i], idx[j]
            if xs[b] != xs[a]:
                sl.append((ys[b] - ys[a]) / (xs[b] - xs[a]))
    if not sl:
        return float("nan"), float("nan")
    sl.sort()
    m = sl[len(sl) // 2]
    ic = sorted(y - m * x for x, y in zip(xs, ys))
    return m, ic[len(ic) // 2]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traces_tsv")
    ap.add_argument("--min-xic", type=int, default=5,
                    help="minimum XIC points; below this an FWHM is grid noise")
    ap.add_argument("--bins", type=int, default=12)
    ap.add_argument("--json", help="write bin medians here (for plotting)")
    a = ap.parse_args()

    rt, fw = [], []
    n_all = n_trunc = n_short = 0
    with open(a.traces_tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        if "fwhm_s" not in (rd.fieldnames or []):
            sys.exit("ABORT: %s has no fwhm_s column -- this is the OLD dump format, rerun the "
                     "tool with the FWHM patch" % a.traces_tsv)
        for x in rd:
            n_all += 1
            try:
                w, r, npts = float(x["fwhm_s"]), float(x["rt"]), int(x["n_xic"])
            except (TypeError, ValueError):
                continue
            if w < 0:
                n_trunc += 1; continue          # never reached half max both sides
            if npts < a.min_xic:
                n_short += 1; continue
            rt.append(r / 60.0); fw.append(w)

    if not fw:
        sys.exit("ABORT: no usable FWHM values")
    order = sorted(range(len(rt)), key=lambda i: rt[i])
    rt = [rt[i] for i in order]; fw = [fw[i] for i in order]

    def q(v, p):
        v = sorted(v); return v[int(p * (len(v) - 1))]

    print("traces: %d total, %d truncated (excluded), %d too few XIC points (excluded), %d used"
          % (n_all, n_trunc, n_short, len(fw)))
    print("FWHM overall: median %.2f s (%.2f cycles)  IQR [%.2f, %.2f]"
          % (q(fw, .5), q(fw, .5) / CYCLE_S, q(fw, .25), q(fw, .75)))
    print()
    print("%-16s %8s %9s %9s %9s %9s" % ("RT bin (min)", "n", "med FWHM", "cycles", "p25", "p75"))
    print("-" * 68)
    out = []
    per = max(1, len(fw) // a.bins)
    for b in range(a.bins):
        lo, hi = b * per, (b + 1) * per if b < a.bins - 1 else len(fw)
        if hi - lo < 30: continue
        seg = fw[lo:hi]; rseg = rt[lo:hi]
        m = q(seg, .5)
        out.append({"rt": sum(rseg) / len(rseg), "fwhm": m,
                    "p25": q(seg, .25), "p75": q(seg, .75), "n": len(seg)})
        print("%6.1f-%-9.1f %8d %9.2f %9.2f %9.2f %9.2f"
              % (rseg[0], rseg[-1], len(seg), m, m / CYCLE_S, q(seg, .25), q(seg, .75)))
    print("-" * 68)

    xs = [o["rt"] for o in out]; ys = [o["fwhm"] for o in out]
    slope, icept = theil_sen(xs, ys)
    lo_rt, hi_rt = xs[0], xs[-1]
    w_lo, w_hi = slope * lo_rt + icept, slope * hi_rt + icept
    print("Theil-Sen: FWHM(s) = %.4f * RT(min) + %.2f" % (slope, icept))
    print("  at RT=%.1f min: %.2f s (%.2f cycles)" % (lo_rt, w_lo, w_lo / CYCLE_S))
    print("  at RT=%.1f min: %.2f s (%.2f cycles)" % (hi_rt, w_hi, w_hi / CYCLE_S))
    print("  ratio across the gradient: %.2fx" % (w_hi / w_lo if w_lo else float("nan")))
    print()
    spread = (max(ys) - min(ys)) / q(fw, .5)
    if abs(w_hi / w_lo - 1.0) < 0.15 and spread < 0.25:
        print("VERDICT: width is effectively FLAT with RT. A single rt_window is defensible,")
        print("         and the review's 'no breakpoint' finding is NOT explained by RT-dependence.")
    else:
        print("VERDICT: width VARIES with RT. A single global rt_window is wrong by construction,")
        print("         and the pooled across-cycle analysis would smear the breakpoint. The merge")
        print("         window must be a function of RT, and that analysis must be redone.")
    print()
    print("We currently spread a peptide over 2.66 cycles; median FWHM is %.2f cycles."
          % (q(fw, .5) / CYCLE_S))

    if a.json:
        open(a.json, "w").write(json.dumps(
            {"bins": out, "slope_s_per_min": slope, "intercept_s": icept,
             "median_fwhm_s": q(fw, .5), "cycle_s": CYCLE_S}, indent=2))
        print("\nwrote %s" % a.json)


if __name__ == "__main__":
    sys.exit(main())
