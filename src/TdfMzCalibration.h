// Copyright (c) 2026, SpeXtract authors. BSD-3-Clause.
//
// Exact TOF -> m/z conversion for Bruker TDF (timsTOF) data, ModelType 1.
//
// Header-only and dependency-free ON PURPOSE: this is the single source of truth for the
// calibration, shared by the OpenMS loader patch (patches/openms-brukertims-mz-calibration.patch,
// which supplies the constants from the file's MzCalibration table) and by the C++ unit test
// (tests/test_calibration_cpp.cpp), so the code that ships is the code that is tested.
//
// Model, derived numerically against Bruker's timsdata library used as a local oracle (no vendor
// code was read, no vendor binary is redistributed; every constant comes from the user's own file):
//
//     t_ns   = tof * DigitizerTimebase + DigitizerDelay
//     C1_eff = C1 * (1 + dC1 * (T1_ref - T1_frame) / 1e6)        // digitizer temperature drift
//     t_ns   = C0 + (1e6 / sqrt(C1_eff)) * sqrt(m) + C2 * m      // solved for sqrt(m)
//
// Verified to 2.5e-5 ppm max against the vendor library (tests/calibration_golden.json).
// SCOPE OF THAT CLAIM, stated precisely because it is narrower than the file count suggests: the
// three files share ONE identical MzCalibration vector (same C0/C1/C2/timebase/delay/T1_ref/dC1),
// so independent coverage of the parameter space is ONE vector -- 12 frames spanning only 0.034 K,
// probed at five TOF positions covering m/z 133.7-1573.7. ModelType 1, single-row table, dC2 == 0,
// one instrument cohort. Everything outside that is REJECTED rather than
// approximated -- see isSupported(): silently emitting wrong masses is the failure mode this whole
// file exists to prevent (the previous two-point chord was -5..-11 ppm and nobody noticed for months).
//
// Note for implementers: the widely-copied open implementation (timsrust-calibration, adapted then
// disabled by mzdata as "not consistently better") DROPS the C2*m term, which is worth -11..-40 ppm
// on this data. That omission -- not translation subtlety -- is why the port underperformed.

#pragma once
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>

namespace spextract
{

/// Constants as stored in the TDF `MzCalibration` row plus the per-frame digitizer temperature.
struct TdfMzCalibration
{
  int model_type = 0;
  double digitizer_timebase = 0.0;
  double digitizer_delay = 0.0;
  double C0 = 0.0;
  double C1 = 0.0;      ///< must be > 0
  double C2 = 0.0;      ///< the quadratic term the open implementations drop
  double T1_ref = 0.0;  ///< reference digitizer temperature for the calibration
  double dC1 = 0.0;     ///< ppm/K drift of C1
  double dC2 = 0.0;     ///< drift of C2 -- NOT modelled; must be 0 (see isSupported)
  double C3 = 0.0, C4 = 0.0;   ///< higher-order terms -- NOT modelled; must be 0

  /// Why this calibration cannot be used, or empty if it can. Fail closed, never approximate.
  std::string unsupportedReason() const
  {
    // NaN defeats ordered comparisons (every < and > below is false for NaN), so screen first.
    for (double v : {digitizer_timebase, digitizer_delay, C0, C1, C2, T1_ref, dC1, dC2, C3, C4})
      if (std::isnan(v)) return "MzCalibration contains NaN";
    if (model_type != 1) return "MzCalibration ModelType " + std::to_string(model_type) + " != 1";
    // A merely-positive C1 is not enough: C1 < ~5.6e-297 overflows b*b to +inf and yields m/z 0.0
    // for every peak, silently, under an "exact" log line. Real values are ~1e5. [claude review]
    if (!(C1 > 1.0) || !(C1 < 1e12)) return "MzCalibration C1 outside the plausible range (1, 1e12)";
    if (dC2 != 0.0) return "MzCalibration dC2 != 0 (temperature drift of C2 is not modelled)";
    if (C3 != 0.0 || C4 != 0.0) return "MzCalibration C3/C4 != 0 (higher-order terms are not modelled)";
    // C2 == 0 is NOT a benign degenerate case: the pure-sqrt law it selects is exactly the
    // known-bad open implementation (2-39 ppm on real vectors). The vendor schema declares
    // C0..C4 with no type and no NOT NULL, so a NULL or text C2 arrives here as 0.0 through
    // sqlite3_column_double -- indistinguishable from a real zero. Reject both. [claude review]
    if (!(C2 > 0.0)) return "MzCalibration C2 <= 0 or missing (the quadratic term is required)";
    if (!(digitizer_timebase > 0.0)) return "DigitizerTimebase <= 0";
    return std::string();
  }
  bool isSupported() const { return unsupportedReason().empty(); }

  /// Temperature-corrected 1e6/sqrt(C1_eff) for one frame; hoist out of per-peak loops.
  double frameFactor(double t1_frame) const
  {
    const double cf = 1.0 + (dC1 * (T1_ref - t1_frame)) / 1e6;
    return 1e6 / std::sqrt(C1 * cf);
  }

  /// TOF index -> m/z. @p b is frameFactor() for this frame.
  double tofToMz(double tof, double b) const
  {
    const double t = tof * digitizer_timebase + digitizer_delay;
    // Stable root of C2*u^2 + b*u + (C0 - t) = 0: the textbook form (-b + sqrt(disc))/(2*C2)
    // subtracts two nearly equal numbers (b ~ 2.5e3, C2 ~ 1e-3), throwing away ~4 digits.
    // u = 2*(t - C0) / (b + sqrt(disc)) is algebraically identical and cancellation-free.
    // t < C0 is unphysical (arrival before the calibration zero) and would yield a NEGATIVE u whose
    // square is a plausible-looking but WRONG mass -- the exact silent-wrongness this file exists to
    // prevent. Both it and a negative discriminant return 0.0, which callers see as an out-of-range
    // peak rather than a credible m/z.
    const double disc = b * b - 4.0 * C2 * (C0 - t);
    // Out-of-model inputs return NaN, never 0.0: a zero m/z is a plausible-looking number that
    // sorts first and propagates silently, which is the failure shape this file exists to remove.
    if (!(disc >= 0.0) || !(t > C0)) return std::numeric_limits<double>::quiet_NaN();
    const double denom = b + std::sqrt(disc);
    if (!(denom > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double u = 2.0 * (t - C0) / denom;
    return u * u;
  }

  /// m/z -> TOF index (exact inverse of tofToMz; the model is closed-form in this direction).
  double mzToTof(double mz, double b) const
  {
    const double t = C0 + b * std::sqrt(mz > 0.0 ? mz : 0.0) + C2 * mz;
    return (t - digitizer_delay) / digitizer_timebase;
  }
};

} // namespace spextract
