#!/usr/bin/env python3
"""Golden test for the TDF MzCalibration (ModelType 1) TOF->m/z model.

Guards the calibration that produced the 2026-09-01 +6-11% closed-search gain. The 60 golden cases
were produced by Bruker's own timsdata library on three diaPASEF files (dataset D/dataset A/dataset B, 4 frames x 5
TOF positions each) and are stored in tests/calibration_golden.json, so this test needs NO vendor
library, NO cluster and NO raw data.

Model (derived numerically against the vendor oracle). NOTE this is an INDEPENDENT re-implementation
of the model, not a copy of the shipped C++: it uses the textbook quadratic root where
src/TdfMzCalibration.h uses the cancellation-free form. That is deliberate -- two independent
implementations agreeing on the vendor numbers is a stronger check than one implementation tested
twice. tests/test_calibration_cpp.cpp is the one that exercises the code that actually ships.

    t_ns   = tof * DigitizerTimebase + DigitizerDelay
    C1_eff = C1 * (1 + dC1 * (T1_ref - T1_frame) / 1e6)          # temperature compensation
    t_ns   = C0 + (1e6 / sqrt(C1_eff)) * sqrt(m) + C2 * m        # solve the quadratic for sqrt(m)

Two regressions this catches:
  * dropping the C2*m term -- the bug in the widely-copied open implementation (timsrust-calibration,
    adapted then disabled by mzdata): -11..-40 ppm on this data;
  * dropping the dC1 temperature term: ~0.67 ppm. (m is proportional to C1_eff, NOT its square
    root, so dC1*dT/1e6 lands on m/z undiminished -- an earlier note here claimed ~0.33 ppm by
    applying a spurious halving; the measured guard value 0.669 confirms 0.67. The golden files
    span only ~0.03 K, so the guard threshold stays loose at 0.1 ppm.)
Run: python3 tests/test_calibration.py
"""
import json, math, os, sys

TOL_PPM = 0.0001         # measured max 2.6e-5 ppm; matches the C++ test's tolerance exactly
NO_C2_MIN_PPM = 5.0      # dropping C2 must be caught: it is worth >>5 ppm
NO_TEMP_MIN_PPM = 0.1    # dropping the temperature term must be visible

def tof_to_mz(tof, p, t1_frame, use_c2=True, use_temp=True):
    t_ns = tof * p["timebase"] + p["delay"]
    c1 = p["C1"] * (1.0 + p["dC1"] * (p["T1_ref"] - t1_frame) / 1e6) if use_temp else p["C1"]
    b = 1e6 / math.sqrt(c1)
    c2 = p["C2"] if use_c2 else 0.0
    if abs(c2) < 1e-12:                      # degenerate: pure sqrt law
        return ((t_ns - p["C0"]) / b) ** 2
    disc = b * b - 4.0 * c2 * (p["C0"] - t_ns)
    u = (-b + math.sqrt(max(disc, 0.0))) / (2.0 * c2)
    return u * u

def ppm(a, b):
    return abs(a - b) / b * 1e6

def check_linear_c2_zero():
    """C2 == 0 stored in the file (PXD017703, 2020 timsTOF Pro) is a real, purely linear-in-sqrt
    calibration and must be accepted; the quadratic collapses to sqrt(m) = (t - C0) / b."""
    p = {"timebase": 0.2, "delay": 25131.0, "C0": 315.70325869866065,
         "C1": 154272.1271422364, "C2": 0.0, "T1": 25.63315397685876, "dC1": -0.2}
    import math
    for tof in (60000.0, 120000.0, 240000.0, 400000.0):
        got = tof_to_mz(tof, p, p["T1"])
        t = tof * p["timebase"] + p["delay"]
        b = 1e6 / math.sqrt(p["C1"])
        expect = ((t - p["C0"]) / b) ** 2
        assert abs(got - expect) / expect < 1e-12, f"C2==0 closed form: tof {tof} got {got} want {expect}"
    print("OK  C2 == 0 (linear law, PXD017703 row) accepted and matches the closed form")

def main():
    check_linear_c2_zero()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration_golden.json")
    files = json.load(open(path))
    n = 0
    worst = 0.0
    worst_noc2 = 0.0
    worst_notemp = 0.0
    for f in files:
        assert f["model_type"] == 1, f"{f['file']}: golden set must be ModelType 1"
        for case in f["cases"]:
            got = tof_to_mz(case["tof"], f, case["t1"])
            e = ppm(got, case["mz"])
            assert e < TOL_PPM, (
                f"{f['file']} frame {case['frame']} tof {case['tof']:.0f}: "
                f"model {got:.6f} vs vendor {case['mz']:.6f} = {e:.5f} ppm > {TOL_PPM}")
            worst = max(worst, e)
            worst_noc2 = max(worst_noc2, ppm(tof_to_mz(case["tof"], f, case["t1"], use_c2=False), case["mz"]))
            worst_notemp = max(worst_notemp, ppm(tof_to_mz(case["tof"], f, case["t1"], use_temp=False), case["mz"]))
            n += 1
    # the ablations must actually be detectable, or this test proves nothing
    assert worst_noc2 > NO_C2_MIN_PPM, (
        f"dropping C2 changed the result by only {worst_noc2:.3f} ppm -- the golden set cannot "
        f"distinguish the known-bad open implementation; regenerate it over a wider m/z range")
    assert worst_notemp > NO_TEMP_MIN_PPM, (
        f"dropping the temperature term changed the result by only {worst_notemp:.4f} ppm -- "
        f"the golden set has no temperature spread; regenerate across more frames")
    print(f"OK  {n} golden cases, {len(files)} files: max error {worst:.6f} ppm (tol {TOL_PPM})")
    print(f"    ablation guards: without C2 {worst_noc2:.2f} ppm, without temperature {worst_notemp:.4f} ppm")
    return 0

if __name__ == "__main__":
    sys.exit(main())
