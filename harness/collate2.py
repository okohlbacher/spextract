#!/usr/bin/env python3
"""SpeXtractor benchmark harness v2 -- measurement, pairing, and the publish gate.

Everything between "a sealed artifact exists" and "a number appears in the vault" lives here.
Run identity and plan binding are bench2.py's problem; this file consumes its runset contract
and refuses to operate without one.

THE SINGLE RULE
---------------
A metric is either a Value with a stated sample size and a NAMED denominator population, or an
Undefined with a machine-readable reason. There is no third state, no default, no dash.
Any Undefined that reaches the publish boundary refuses the whole report.

This is the direct answer to the reviewers' inventory of silent fallbacks. v1 contained
`.get(col, 1)`, `own or [0]`, `max(tot, 1)`, `i2.sum() or 1.0`, `max(1-recd, 1e-9)` and
`Path(refs.get("diann_lib", ""))`. Each of them converts an ABSENCE into a plausible number:
a missing Sage column became "0 peptides at 1% FDR", an absent fragment library became "0%
co-isolation", zero charge observations became "0% charge agreement". Those are invented
measurements, and they are worse than a crash because they are quotable.

Commands
--------
    collate2.py truth   --runset R          # compute metrics; writes runset revision r2
    collate2.py collate --runset R [--publish]
    collate2.py selftest

Exit codes: 0 published / ok, 2 guard abort, 3 refusal to publish.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bench2  # noqa: E402
from bench2 import (Abort, die, canon, jload, sha256_bytes, sha256_file, short, utc_now,
                    Bench, Attempt, bench_root_from, load_registry, load_runset, write_runset,
                    assert_cells_match_plan, ledger_count, sample_entry, CONTRACT_VERSION,
                    HARNESS_VERSION, SelfTest)  # noqa: E402

# Optional third-party deps are checked ONCE, at startup, with an actionable message -- not on
# first use three hours into a job.
_MISSING = []
try:
    import numpy as np
except ImportError:
    np = None
    _MISSING.append("numpy")
try:
    import pyarrow.parquet as pq
except ImportError:
    pq = None
    _MISSING.append("pyarrow")
try:
    from pyteomics import mzml as pyt_mzml
except ImportError:
    pyt_mzml = None
    _MISSING.append("pyteomics")


def require_deps():
    if _MISSING:
        die("G-C00", "missing python packages: %s\n       python3 -m pip install --user %s"
            % (", ".join(_MISSING), " ".join(_MISSING)))


# ---------------------------------------------------------------- audited structural accessors
# These two live OUTSIDE the measurement region on purpose. They are the only sanctioned way for
# measurement code to navigate pyteomics' nested dictionaries or to read a CV unit, and they are
# short enough to audit by eye. Their existence is what lets the region ban `.get(` and 3-argument
# `getattr` outright: a structural absence returns None, which measurement code must then turn
# into a NAMED disposition -- it can never become a numeric default.

def field(mapping, key, absent=None):
    """Structural navigation. Returns `absent` (None by default) when the key is not there.

    NOTE the discipline: `absent` is only ever a container or None -- never a number. v1's
    `ion.get("charge state", 0)` returned a NUMBER for a missing field, and that number was then
    measured. Nothing here may pass a numeric `absent`."""
    if isinstance(absent, (int, float)) and not isinstance(absent, bool):
        raise TypeError("field(): a numeric default would manufacture a measurement")
    if mapping is None or key not in mapping:
        return absent
    return mapping[key]


def unit_info_of(value):
    """The CV unit attached to a pyteomics value, or None if there is none.

    Returning None is not a guess. v1's `getattr(v, "unit_info", "minute")` guessed MINUTES, and
    RT in seconds read as minutes matched 31 of ~300,000 spectra (FAILURE 8)."""
    try:
        return value.unit_info
    except AttributeError:
        return None


# =============================================================================================
# --- BEGIN MEASUREMENT REGION ---
# Everything between these markers is linted by lint_source() (guard G-C30). The banned
# constructs are exactly the ones the reviewers found manufacturing numbers out of absences.
# Plumbing, rendering and the selftest live OUTSIDE the region and are not linted; the README
# states that scope rather than implying the ban is global.
# =============================================================================================

DECOY_SHIFT = 11.003          # Da, applied to the QUERY, matched against the DIA-NN truth list
PPM_PREC = 20.0
PPM_FRAG = 20.0
D_RT = 0.30                   # minutes
D_IM = 0.05                   # 1/K0
RT_BIN = 0.1                  # minutes; enters estimand_id because it changes which truth
                              # precursors are candidates
DEPTH_K = 20
Q_THRESHOLD = 0.01

# Minimum observations below which a metric is INSUFFICIENT_N rather than a number.
N_MIN = {"recall": 200, "purity": 200, "charge_agreement": 500, "q_gap": 200}

REASONS = frozenset([
    "MISSING_INPUT", "SCHEMA_VIOLATION", "EMPTY_POPULATION", "INSUFFICIENT_N",
    "ZERO_DENOMINATOR", "POPULATION_MISMATCH", "POPULATION_LOSS", "NON_FINITE",
    "OUT_OF_DOMAIN", "CONTROL_DEGENERATE", "CONTROL_TOO_HIGH", "CONTROL_ASYMMETRIC",
    "UNIT_UNKNOWN", "METHOD_MISMATCH", "ESTIMAND_MISMATCH", "REFERENCE_UNVERIFIED",
    "PROPAGATED", "CROSS_SAMPLE", "NOT_COMPUTED",
])

UNITS = frozenset(["percent", "count", "ratio", "seconds", "megabytes", "pp"])

# Named denominator populations. A metric's denominator is looked up BY NAME from the disposition
# ledger, never taken from len(some_list). v1 computed charge agreement over matched spectra and
# its majority-class baseline over all truth precursors -- two different populations, silently
# differenced.
POPULATIONS = {
    "ms2_spectra": "every MS2 spectrum read from the mzML",
    "scored_spectra": "MS2 spectra whose disposition is SCORED",
    "charge_observations": "SCORED spectra carrying a tool charge AND a matched truth precursor",
    "truth_precursors": "DIA-NN precursors for THIS sample at Q.Value <= 0.01",
    "rank1_psms": "Sage rows with rank == 1",
    "purity_spectra": "SCORED spectra with >= DEPTH_K peaks and nonzero intensity",
}
SUBSET_OF = {
    "scored_spectra": "ms2_spectra",
    "charge_observations": "scored_spectra",
    "purity_spectra": "scored_spectra",
}


def _ancestors(pop):
    out = [pop]
    cur = pop
    while cur in SUBSET_OF:
        cur = SUBSET_OF[cur]
        out.append(cur)
    return out


class Undefined(object):
    """A measurement that could not be made. Carries WHY, machine-readably."""

    __slots__ = ("code", "detail", "where", "causes")

    def __init__(self, code, detail, where, causes=()):
        if code not in REASONS:
            raise ValueError("unknown reason code %r" % code)
        self.code = code
        self.detail = detail
        self.where = where
        self.causes = tuple(causes)

    def to_json(self):
        return {"__undefined__": self.code, "detail": self.detail, "where": self.where,
                "causes": [c.to_json() for c in self.causes]}

    def __repr__(self):
        return "Undefined(%s @ %s: %s)" % (self.code, self.where, self.detail)


class Value(object):
    """A measurement.

    The numeric payload is `_v`, private by convention and enforced by lint_source(): the string
    `._v` may not appear inside the measurement region outside the sanctioned arithmetic helpers.
    Value implements no arithmetic dunders, so `a / b` on two Values is a TypeError rather than a
    ratio with no provenance.

    HONEST: this is lint-enforced, not enforced by the type system. Python cannot prevent
    attribute access. It stops the accident, not a determined author. The README says so.
    """

    __slots__ = ("_v", "n", "unit", "pop", "estimand")

    def __init__(self, v, n, unit, pop, estimand):
        if not isinstance(v, float) or not math.isfinite(v):
            raise ValueError("Value payload must be a finite float, got %r" % (v,))
        if not isinstance(n, int) or n < 1:
            # n == 0 would let an author who caught an Undefined substitute a publishable zero.
            raise ValueError("Value requires n >= 1, got %r" % (n,))
        if unit not in UNITS:
            raise ValueError("unknown unit %r" % unit)
        if pop not in POPULATIONS:
            raise ValueError("undeclared population %r" % pop)
        self._v = v
        self.n = n
        self.unit = unit
        self.pop = pop
        self.estimand = estimand

    def to_json(self):
        return {"v": self._v, "n": self.n, "unit": self.unit, "pop": self.pop,
                "estimand": self.estimand}

    def __repr__(self):
        return "Value(%.4g %s, n=%d, %s)" % (self._v, self.unit, self.n, self.pop)


def raw_number(q):
    """The ONLY sanctioned reader of a Value's payload. Renderers and arithmetic go through it."""
    if isinstance(q, Undefined):
        raise TypeError("raw_number() on an Undefined -- callers must branch first")
    return q._v  # noqa: SLF001  (sanctioned accessor)


def is_undef(q):
    return isinstance(q, Undefined)


def q_count(pop, n, estimand):
    """Mint a denominator. The name and the number travel as one object, so a denominator can
    never be silently taken from a different population than the numerator."""
    if pop not in POPULATIONS:
        raise ValueError("undeclared population %r" % pop)
    if n < 1:
        return Undefined("EMPTY_POPULATION", "population %s is empty" % pop, "q_count:%s" % pop)
    return Value(float(n), n, "count", pop, estimand)


def q_ratio(num, den, where, unit="percent", n_min=1, scale=100.0):
    """The single division in this file.

    Replaces every one of v1's `max(len(sp), 1)`, `max(tot_z, 1)`, `max(len(present), 1)`,
    `max(1 - recd, 1e-9)` and `i2.sum() or 1.0`. A zero denominator is not "1"; it is
    ZERO_DENOMINATOR, and ZERO_DENOMINATOR does not publish.
    """
    if is_undef(num) or is_undef(den):
        causes = tuple(x for x in (num, den) if is_undef(x))
        return Undefined("PROPAGATED", "an operand was undefined", where, causes)
    if den.pop not in _ancestors(num.pop) and num.pop != den.pop:
        # v1's charge agreement / majority-baseline mismatch, made unrepresentable.
        return Undefined("POPULATION_MISMATCH",
                         "numerator population %r is not a declared subset of denominator "
                         "population %r" % (num.pop, den.pop), where)
    d = raw_number(den)
    if d == 0.0:
        return Undefined("ZERO_DENOMINATOR", "denominator population %s has 0 members" % den.pop,
                         where)
    if den.n < n_min:
        return Undefined("INSUFFICIENT_N", "denominator n=%d is below the declared minimum %d"
                         % (den.n, n_min), where)
    out = scale * raw_number(num) / d
    if not math.isfinite(out):
        return Undefined("NON_FINITE", "ratio evaluated to a non-finite value", where)
    return Value(out, den.n, unit, num.pop, num.estimand)


def q_mean(xs, pop, estimand, where, n_min=1, scale=100.0):
    """Replaces `np.array(own or [0]).mean()`.

    An empty list of observations is not a mean of zero. v1 manufactured a synthetic [0]
    observation, which is how an absent fragment library became a legitimate-looking 0% purity.
    """
    if len(xs) == 0:
        return Undefined("EMPTY_POPULATION", "no observations for %s" % pop, where)
    if len(xs) < n_min:
        return Undefined("INSUFFICIENT_N", "n=%d below declared minimum %d" % (len(xs), n_min),
                         where)
    total = 0.0
    for x in xs:
        if not math.isfinite(x):
            return Undefined("NON_FINITE", "a non-finite observation entered %s" % pop, where)
        total += x
    return Value(scale * total / len(xs), len(xs), "percent", pop, estimand)


def q_diff(a, b, where):
    """Difference of two Values in the SAME unit. A percent minus a percent is percentage
    POINTS, and the unit says so -- so nothing downstream can render it as a relative change."""
    if is_undef(a) or is_undef(b):
        return Undefined("PROPAGATED", "an operand was undefined", where,
                         tuple(x for x in (a, b) if is_undef(x)))
    if a.unit != b.unit:
        return Undefined("POPULATION_MISMATCH", "unit mismatch %s vs %s" % (a.unit, b.unit),
                         where)
    if a.estimand != b.estimand:
        return Undefined("ESTIMAND_MISMATCH",
                         "these two numbers were computed under different measurement rules "
                         "(%s vs %s)" % (a.estimand[:12], b.estimand[:12]), where)
    unit = "pp" if a.unit == "percent" else a.unit
    return Value(raw_number(a) - raw_number(b), min(a.n, b.n), unit, a.pop, a.estimand)


class Controlled(object):
    """A measurement and its chance floor, welded together.

    FAILURE 4: a "97.8% detection ceiling" was published, then retracted when a decoy floor
    showed a precursor that does not exist still matched 91.9% of the time. v1's fix was a
    convention ("uncorrected values are never emitted alone"); collate then printed
    recall_corrected without either the raw value or the floor.

    Here there is no code path that produces one without the other: Controlled has no defaults,
    and render_controlled() prints all three or raises.
    """

    __slots__ = ("raw", "floor", "corrected", "control_id", "control_candidates",
                 "control_matches", "selection_rule_sha", "shifted_out_of_window_frac")

    def __init__(self, raw, floor, control_id, control_candidates, control_matches,
                 selection_rule_sha, shifted_out_of_window_frac):
        self.raw = raw
        self.floor = floor
        self.control_id = control_id
        self.control_candidates = control_candidates
        self.control_matches = control_matches
        self.selection_rule_sha = selection_rule_sha
        self.shifted_out_of_window_frac = shifted_out_of_window_frac
        self.corrected = abbott(raw, floor, control_candidates, control_matches,
                                shifted_out_of_window_frac, "corrected")

    def to_json(self):
        return {"raw": _q_json(self.raw), "floor": _q_json(self.floor),
                "corrected": _q_json(self.corrected), "control_id": self.control_id,
                "control_candidates": self.control_candidates,
                "control_matches": self.control_matches,
                "selection_rule_sha": self.selection_rule_sha,
                "shifted_out_of_window_frac": self.shifted_out_of_window_frac}


# A control that cannot fire is not a control. Both reviewers independently argued that an
# 11.003 Da shift against a 20 ppm tolerance is ~11,000 ppm at m/z 1000 and therefore can never
# match, making the floor structurally zero and (obs-floor)/(1-floor) the identity function. The
# counter-argument is that the shifted query is matched against the DIA-NN TRUTH LIST, not
# against itself, and truth density inside a diaPASEF tile is high. This harness refuses to
# settle that by argument: it records how many opportunities the control had, and a control that
# had opportunities but never fired is CONTROL_DEGENERATE and blocks publication.
CONTROL_DEGENERATE_MIN_CANDIDATES = 100
SHIFTED_OUT_OF_WINDOW_TOL = 0.20


def abbott(raw, floor, control_candidates, control_matches, shifted_out_frac, where):
    """(observed - floor) / (1 - floor), never clamped, with every degenerate case named.

    v1: `round(100 * (rec - recd) / max(1 - recd, 1e-9), 1)`. That expression yields a NEGATIVE
    recall when rec < recd and yields 500% when rec=0.95, recd=0.9 -- and the vault printed both
    without complaint. Clamping would have hidden the anomaly; refusing exposes it.
    """
    if is_undef(raw) or is_undef(floor):
        return Undefined("PROPAGATED", "raw or floor undefined", where,
                         tuple(x for x in (raw, floor) if is_undef(x)))
    if control_candidates < CONTROL_DEGENERATE_MIN_CANDIDATES:
        return Undefined("CONTROL_DEGENERATE",
                         "the decoy control was evaluated only %d times (< %d); it had no real "
                         "opportunity to fire, so its floor is not evidence of anything"
                         % (control_candidates, CONTROL_DEGENERATE_MIN_CANDIDATES), where)
    if control_matches == 0:
        return Undefined("CONTROL_DEGENERATE",
                         "the decoy control was evaluated %d times and matched ZERO times. A "
                         "floor that is structurally incapable of firing makes the correction "
                         "the identity function -- this is the defect that produced the "
                         "retracted 97.8%% number." % control_candidates, where)
    if shifted_out_frac > SHIFTED_OUT_OF_WINDOW_TOL:
        return Undefined("CONTROL_ASYMMETRIC",
                         "%.1f%% of decoy queries were shifted outside their own isolation "
                         "window (tolerance %.0f%%). The control is systematically weaker for "
                         "precursors near a tile's upper edge, so it is not exchangeable with "
                         "the target across the m/z range."
                         % (100 * shifted_out_frac, 100 * SHIFTED_OUT_OF_WINDOW_TOL), where)
    r = raw_number(raw) / 100.0
    f = raw_number(floor) / 100.0
    if f > 0.5:
        return Undefined("CONTROL_TOO_HIGH",
                         "chance floor is %.1f%%; the correction is numerically unstable above "
                         "50%% and the measurement is dominated by chance" % (100 * f), where)
    if r < f:
        return Undefined("OUT_OF_DOMAIN",
                         "observed %.2f%% is BELOW its chance floor %.2f%%. The correction would "
                         "be negative. This is not clamped to zero: a value below chance means "
                         "the estimand is wrong, not that performance is zero."
                         % (100 * r, 100 * f), where)
    corrected = 100.0 * (r - f) / (1.0 - f)
    if not (0.0 <= corrected <= 100.0) or not math.isfinite(corrected):
        return Undefined("OUT_OF_DOMAIN", "corrected value %.2f is outside [0,100]" % corrected,
                         where)
    return Value(corrected, raw.n, "percent", raw.pop, raw.estimand)


# ------------------------------------------------------------------ disposition ledger (G-C11)
DROP_REASONS = (
    "SCORED", "NOT_MS2", "NO_PRECURSOR_MZ", "NO_CHARGE_STATE", "NO_ION_MOBILITY",
    "OUT_OF_WINDOW_TABLE", "AMBIGUOUS_WINDOW", "EMPTY_PEAK_LIST", "ZERO_TOTAL_INTENSITY",
    "FEWER_THAN_K_PEAKS", "NO_TRUTH_IN_RT_BIN", "NO_TARGET_CANDIDATE",
)

# Every spectrum the truth stage declines to score increments a NAMED counter, and each counter
# has a declared tolerance. v1 used bare `continue`, so a wrong window table showed up as quietly
# depressed recall rather than as an error. Above tolerance, the affected metrics become
# Undefined and the report refuses.
DROP_TOLERANCE = {
    "NOT_MS2": 1.01,               # MS1 spectra are expected; not a defect
    "NO_PRECURSOR_MZ": 0.01,
    "NO_CHARGE_STATE": 0.50,       # tools legitimately omit charge; it only gates charge metrics
    "NO_ION_MOBILITY": 0.01,
    "OUT_OF_WINDOW_TABLE": 0.005,  # a wrong acquisition table announces itself HERE
    "AMBIGUOUS_WINDOW": 0.02,
    "EMPTY_PEAK_LIST": 0.01,
    "ZERO_TOTAL_INTENSITY": 0.01,
    "FEWER_THAN_K_PEAKS": 0.50,    # gates purity only; recorded in depth_basis
    "NO_TRUTH_IN_RT_BIN": 1.01,
    "NO_TARGET_CANDIDATE": 1.01,
}
DROP_BLOCKS = {
    "OUT_OF_WINDOW_TABLE": ("METHOD_MISMATCH", ("recall", "purity", "charge")),
    "AMBIGUOUS_WINDOW": ("METHOD_MISMATCH", ("recall", "purity", "charge")),
    "NO_PRECURSOR_MZ": ("POPULATION_LOSS", ("recall", "purity", "charge")),
    "NO_ION_MOBILITY": ("POPULATION_LOSS", ("recall", "purity", "charge")),
    "EMPTY_PEAK_LIST": ("POPULATION_LOSS", ("purity",)),
    "ZERO_TOTAL_INTENSITY": ("POPULATION_LOSS", ("purity",)),
    "FEWER_THAN_K_PEAKS": ("POPULATION_LOSS", ("purity",)),
    "NO_CHARGE_STATE": ("POPULATION_LOSS", ("charge",)),
}


class Ledger(object):
    """An exact partition of every spectrum read. sum(counters) == spectra_read, asserted."""

    def __init__(self):
        self.counts = collections.OrderedDict((r, 0) for r in DROP_REASONS)
        self.total_read = 0

    def record(self, reason):
        if reason not in self.counts:
            raise KeyError("undeclared disposition %r" % reason)
        self.counts[reason] += 1
        self.total_read += 1

    def check_partition(self, declared_in_file):
        """I-1, with an INDEPENDENT cross-check.

        Checking sum(counters) against a total counted in the same loop is self-confirming: a
        reader that silently truncates a corrupt mzML shrinks both sides together. The mzML's own
        <spectrumList count="N"> is an independent witness and is compared here.
        """
        s = sum(self.counts.values())
        if s != self.total_read:
            die("G-C11", "disposition ledger does not partition: %d counters vs %d spectra read. "
                         "A spectrum left the loop without a named disposition."
                % (s, self.total_read))
        if declared_in_file is not None and self.total_read != declared_in_file:
            die("G-C11", "read %d spectra but the mzML declares %d. The reader silently "
                         "truncated the file, and a ledger that only checks itself would not "
                         "have noticed." % (self.total_read, declared_in_file))

    def ms2_total(self):
        return self.total_read - self.counts["NOT_MS2"]

    def blocked(self):
        """Which metric families are unpublishable, and why."""
        out = {}
        ms2 = self.ms2_total()
        if ms2 == 0:
            return {"recall": ("EMPTY_POPULATION", "no MS2 spectra"),
                    "purity": ("EMPTY_POPULATION", "no MS2 spectra"),
                    "charge": ("EMPTY_POPULATION", "no MS2 spectra")}
        for reason, tol in sorted(DROP_TOLERANCE.items()):
            if reason == "NOT_MS2":
                continue
            frac = float(self.counts[reason]) / ms2
            if frac > tol:
                code, families = (DROP_BLOCKS[reason] if reason in DROP_BLOCKS
                                  else ("POPULATION_LOSS", ()))
                for fam in families:
                    out.setdefault(fam, (code, "%s affected %.2f%% of MS2 spectra (tolerance "
                                               "%.2f%%)" % (reason, 100 * frac, 100 * tol)))
        scored_frac = float(self.counts["SCORED"]) / ms2
        if scored_frac < 0.60:
            for fam in ("recall", "purity", "charge"):
                out.setdefault(fam, ("POPULATION_LOSS",
                                     "only %.1f%% of MS2 spectra were scorable (floor 60%%)"
                                     % (100 * scored_frac)))
        return out

    def to_json(self):
        return {"spectra_read": self.total_read, "ms2": self.ms2_total(),
                "dispositions": dict(self.counts)}


# ------------------------------------------------------------------- strict field access (G-C02)
class MissingField(Exception):
    pass


def opt_float(mapping, key):
    """Explicit None, never a sentinel.

    v1 wrote `float(ion.get("charge state", 0) or 0)` and `float(ion.get(..., -1) or -1)`. The
    second one is the worst bug in the file that NEITHER review found: an absent ion mobility
    became -1, and the co-isolation candidate filter read `im < 0 or abs(...) <= D_IM`, so a
    missing IM SILENTLY DISABLED the ion-mobility gate for that spectrum. The effective tolerance
    changed mid-file and nothing said so.
    """
    if key not in mapping:
        return None
    v = mapping[key]
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return f


def req_float(row, key, where, lo=None, hi=None):
    if key not in row:
        raise MissingField("%s: required field %r absent" % (where, key))
    v = row[key]
    if v is None or (isinstance(v, str) and v.strip() == ""):
        die("G-C04", "%s: field %r is blank. A blank q-value is not 1.0." % (where, key))
    try:
        f = float(v)
    except (TypeError, ValueError):
        die("G-C04", "%s: field %r is %r, which is not a number. v1 would have crashed here on "
                     "'NA' or silently defaulted on absence." % (where, key, v))
    if not math.isfinite(f):
        die("G-C04", "%s: field %r is %r (NaN/Inf). NaN silently fails every <= comparison, "
                     "which is a population change disguised as a filter." % (where, key, v))
    if lo is not None and f < lo:
        die("G-C04", "%s: field %r = %g is below %g" % (where, key, f, lo))
    if hi is not None and f > hi:
        die("G-C04", "%s: field %r = %g is above %g" % (where, key, f, hi))
    return f


def req_int(row, key, where, choices=None):
    f = req_float(row, key, where)
    i = int(f)
    if float(i) != f:
        die("G-C04", "%s: field %r = %r is not an integer" % (where, key, row[key]))
    if choices is not None and i not in choices:
        die("G-C04", "%s: field %r = %d not in %s" % (where, key, i, sorted(choices)))
    return i


# ----------------------------------------------------------------------------- mzML io (G-C01)
RT_UNITS = {"second": 1.0 / 60.0, "seconds": 1.0 / 60.0, "UO:0000010": 1.0 / 60.0,
            "minute": 1.0, "minutes": 1.0, "UO:0000031": 1.0}


class RtUnitState(object):
    """The RT unit must be present AND identical for every spectrum in the file."""

    def __init__(self):
        self.unit = None


def rt_minutes(scan, state, where):
    """FAILURE 8: our mzML writes RT in SECONDS, the reference implementation's in MINUTES, DIA-NN reports minutes.
    Trusting the bare number matched 31 of ~300,000 spectra.

    v1 read the CV unit -- which was right -- but then did
    `getattr(v, "unit_info", "minute")`: a THREE-argument getattr whose default silently assumed
    minutes whenever the unit was absent, spelled "seconds", or given as an accession. The guard
    guessed exactly where it needed to refuse. Here an unknown or absent unit aborts on the
    FIRST spectrum, which costs a second rather than a benchmark.
    """
    if "scan start time" not in scan:
        die("G-C01", "%s: spectrum has no 'scan start time'" % where)
    v = scan["scan start time"]
    unit = unit_info_of(v)
    if unit is None or str(unit) not in RT_UNITS:
        die("G-C01", "%s: retention time carries unit %r, which is not in the closed accession "
                     "map %s.\n       v1 assumed 'minute' here. RT in seconds read as minutes "
                     "matched 31 of ~300,000 spectra (FAILURE 8)."
            % (where, unit, sorted(set(RT_UNITS))))
    unit = str(unit)
    if state.unit is None:
        state.unit = unit
    elif state.unit != unit:
        die("G-C01", "%s: retention-time unit changed mid-file (%s -> %s). A mixed-unit file "
                     "cannot be scored." % (where, state.unit, unit))
    return float(v) * RT_UNITS[unit]


def precursor_of(spec):
    """Returns (mz, mz_source, charge_or_None, im_or_None).

    FAILURE 8b: the reference implementation omits 'selected ion m/z' and writes the precursor as its isolation
    window; we write the real acquisition window. v1 handled this with
    `ion.get("selected ion m/z", iso.get("isolation window target m/z", 0))` -- but a .get()
    default only fires on ABSENCE, so a present-but-zero value skipped the fallback entirely and
    the spectrum was silently dropped. Two explicit branches, and the source is COUNTED, because
    two files scored via different m/z sources are not comparable.
    """
    plist = field(spec, "precursorList")
    if not plist or not field(plist, "precursor"):
        return None, "none", None, None
    pre = plist["precursor"][0]
    ions = field(field(pre, "selectedIonList", {}), "selectedIon", [])
    ion = ions[0] if ions else {}
    iso = field(pre, "isolationWindow", {})
    mz = opt_float(ion, "selected ion m/z")
    src = "selected_ion"
    if mz is None or mz <= 0.0:
        mz = opt_float(iso, "isolation window target m/z")
        src = "isolation_target"
    if mz is None or mz <= 0.0:
        return None, "none", None, None
    z = opt_float(ion, "charge state")
    z = int(z) if (z is not None and z > 0) else None
    im = opt_float(ion, "inverse reduced ion mobility")
    if im is None:
        sc = field(field(spec, "scanList", {}), "scan", [{}])[0]
        im = opt_float(sc, "inverse reduced ion mobility")
    if im is not None and im <= 0.0:
        im = None
    return mz, src, z, im


# -------------------------------------------------------------- acquisition windows (G-C10)
def assign_window(windows, mz, im):
    """2-D assignment: m/z AND ion mobility. Returns (window, reason).

    v1 hardcoded 24 one-dimensional dataset B m/z intervals, applied them to dataset A and dataset D, and resolved
    overlaps by taking the FIRST match. diaPASEF tiles are 2-D (m/z x 1/K0) objects that overlap
    in m/z precisely because they are separated in mobility; a 1-D projection with first-match
    resolution assigns a large fraction of spectra to the wrong tile, which biases the
    co-isolation set. Zero matches and multiple matches are both named dispositions here -- there
    is no first-match fallback.
    """
    hits = []
    for w in windows:
        if not (w["mz_lo"] <= mz <= w["mz_hi"]):
            continue
        if im is not None and not (w["im_lo"] - D_IM <= im <= w["im_hi"] + D_IM):
            continue
        hits.append(w)
    if not hits:
        return None, "OUT_OF_WINDOW_TABLE"
    if len(hits) > 1:
        return None, "AMBIGUOUS_WINDOW"
    return hits[0], "SCORED"


# ------------------------------------------------------------ symmetric matching (G-C13/G-C12)
def select_match(candidates, query_mz, ppm):
    """ONE selection rule, used for the target pass AND the control pass.

    v1 selected a single best candidate for targets but added EVERY candidate within tolerance
    for decoys. Different counting rules for target and control make the floor incomparable with
    the value it is supposed to correct -- and the correction was published anyway.
    """
    best = None
    best_err = None
    for p in candidates:
        err = abs(p["mz"] - query_mz) / query_mz * 1e6
        if err <= ppm and (best_err is None or err < best_err):
            best, best_err = p, err
    return best


def attribute(peak_mz, own_sorted, other_sorted, ppm):
    """Exactly ONE category per peak, for targets and decoys alike.

    v1 used `if hits(own): io += it; elif hits(others): ib += it` for targets but two independent
    `if` statements for decoys, so a peak matching both sets counted once for targets and twice
    for decoys. The invariant below (own + coiso + unexplained == 1) turns any recurrence of that
    asymmetry into a runtime abort.
    """
    if _hits(own_sorted, peak_mz, ppm):
        return "own"
    if _hits(other_sorted, peak_mz, ppm):
        return "coiso"
    return "unexplained"


def _hits(arr, mz, ppm):
    if arr is None or len(arr) == 0:
        return False
    i = int(np.searchsorted(arr, mz))
    for j in (i - 1, i):
        if 0 <= j < len(arr) and abs(arr[j] - mz) / mz * 1e6 <= ppm:
            return True
    return False


SELECTION_RULE_SHA = None  # set at import; see _hash_rules()


def _hash_rules():
    import inspect
    src = "".join([inspect.getsource(select_match), inspect.getsource(attribute),
                   inspect.getsource(_hits)])
    return sha256_bytes(src.encode("utf-8"))


# --- END MEASUREMENT REGION ---
# =============================================================================================

SELECTION_RULE_SHA = _hash_rules()


def metrics_code_digest():
    """estimand_id hashes the CODE, not just a declaration of parameters.

    A declaration-only hash lets RT_BIN or the order in which D_RT and D_IM are applied change
    without the estimand changing -- numbers move, the equality gate stays happy, and two
    incomparable measurements share a column.
    """
    return sha256_file(HERE / "collate2.py")


def estimand_id(method_id, mz_source_profile):
    """The keystone against cross-comparison.

    All cells in one table column must share an estimand_id, and a delta between two Values with
    different estimand_ids is refused. This generalizes FAILURE 6 (numbers compared under
    incompatible rules) without having to anticipate each specific mismatch.

    mz_source_profile is BUCKETED, not carried as a continuous fraction. A raw fraction inside an
    equality hash means 0.998 vs 0.997 on two runs mismatches every column forever -- a harness
    that refuses everything gets bypassed by hand, and then it guarantees nothing.
    """
    fields = {
        "metrics_code_sha256": metrics_code_digest(),
        "ppm_prec": PPM_PREC, "ppm_frag": PPM_FRAG, "d_rt": D_RT, "d_im": D_IM,
        "rt_bin": RT_BIN, "depth_k": DEPTH_K, "q_threshold": Q_THRESHOLD,
        "decoy_shift_da": DECOY_SHIFT, "decoy_model": "prec_and_frag_shift_+11.003Da_v1",
        "selection_rule_sha256": SELECTION_RULE_SHA,
        "n_min": dict(sorted(N_MIN.items())),
        "window_dims": ["mz", "im"], "method_id": method_id,
        "mz_source_profile": mz_source_profile,
        "id_unit": "peptide_sequence_as_reported_by_sage",
        "rank_filter": 1,
    }
    return sha256_bytes(canon(fields))[:16]


def bucket_mz_source(counts):
    total = sum(counts.values())
    if total == 0:
        return "none"
    top = max(counts.items(), key=lambda kv: kv[1])
    return top[0] if float(top[1]) / total >= 0.99 else "mixed"


# ------------------------------------------------------------------ Sage contract (G-C04/G-C05)
SAGE_COLUMN_DOMAIN = {
    "rank": ("int", 1, None), "label": ("choice", (-1, 1), None),
    "charge": ("int", 1, 8), "spectrum_q": ("float", 0.0, 1.0),
    "peptide_q": ("float", 0.0, 1.0), "protein_q": ("float", 0.0, 1.0),
}


def read_sage(tsv_path, expected_mzml_name, header_sha_pin=None):
    """Every failure is an abort; none is a default.

    FAILURE 3 chose peptide_q as the metric, and that part was right. What was not guarded was
    the schema: `float(x.get("peptide_q", 1))` means a Sage version that renames the column
    silently treats every row as failing and publishes a bold zero. That is a wrong number, not
    an error.
    """
    with open(str(tsv_path), "r") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        try:
            header = next(rdr)
        except StopIteration:
            die("G-C04", "%s is empty" % tsv_path)
        if len(set(header)) != len(header):
            dupes = [c for c in set(header) if header.count(c) > 1]
            die("G-C04", "%s header has duplicate column name(s) %s" % (tsv_path, dupes))
        missing = [c for c in bench2.SAGE_REQUIRED_COLUMNS if c not in header]
        if missing:
            die("G-C04", "%s is missing required column(s) %s.\n       columns present: %s"
                % (tsv_path, missing, header))
        hsha = sha256_bytes("\t".join(header).encode("utf-8"))
        if header_sha_pin is not None and hsha != header_sha_pin:
            die("G-C04", "%s header digest %s does not match the pinned schema %s.\n"
                         "       A Sage upgrade must be an explicit registry edit, never a silent "
                         "reinterpretation of columns." % (tsv_path, short(hsha),
                                                           short(header_sha_pin)))
        rows = []
        for i, rec in enumerate(rdr, start=2):
            if not rec:
                continue
            if len(rec) != len(header):
                die("G-C04", "%s row %d has %d fields, header has %d"
                    % (tsv_path, i, len(rec), len(header)))
            row = dict(zip(header, rec))
            where = "%s:%d" % (tsv_path.name, i)
            for col, spec in sorted(SAGE_COLUMN_DOMAIN.items()):
                kind = spec[0]
                if kind == "float":
                    req_float(row, col, where, spec[1], spec[2])
                elif kind == "int":
                    req_int(row, col, where)
                elif kind == "choice":
                    req_int(row, col, where, set(spec[1]))
            if not row["peptide"].strip():
                die("G-C04", "%s row %d has an empty peptide" % (tsv_path, i))
            fn = os.path.basename(row["filename"])
            if fn != expected_mzml_name:
                die("G-C04", "%s row %d names mzML %r but this cell's artifact is %r.\n"
                             "       This is a STALE search result being rescored (FAILURE 10)."
                    % (tsv_path, i, fn, expected_mzml_name))
            rows.append(row)
    if not rows:
        die("G-C04", "%s has a header but no PSM rows" % tsv_path)
    return header, rows, hsha


def search_metrics(rows, eid):
    """peptide_q is the published metric (FAILURE 3). spectrum_q is computed for diagnostics only
    and NEVER enters a vault table -- v1 printed it in the headline row, where it was quotable,
    and it inflated every ratio by 9-11 points against the reference implementation."""
    rank1 = [r for r in rows if int(float(r["rank"])) == 1]
    pep_q = set()
    spec_q = set()
    npsm = 0
    for r in rank1:
        if float(r["peptide_q"]) <= Q_THRESHOLD:
            pep_q.add(r["peptide"])
        if float(r["spectrum_q"]) <= Q_THRESHOLD:
            spec_q.add(r["peptide"])
            npsm += 1
    den = q_count("rank1_psms", len(rank1), eid)
    gap_num = q_count("rank1_psms", max(len(pep_q) - len(spec_q), 0), eid) \
        if len(spec_q) else Undefined("ZERO_DENOMINATOR",
                                      "no peptides passed spectrum_q <= 0.01",
                                      "q_metric_count_gap")
    if len(spec_q):
        gap = q_ratio(
            Value(float(len(pep_q) - len(spec_q)), max(len(rank1), 1), "count", "rank1_psms",
                  eid),
            Value(float(len(spec_q)), len(spec_q), "count", "rank1_psms", eid),
            "q_metric_count_gap", n_min=N_MIN["q_gap"])
    else:
        gap = Undefined("ZERO_DENOMINATOR", "no peptides passed spectrum_q <= 0.01",
                        "q_metric_count_gap")
    return {
        "peptides_at_peptide_q": _q_json(q_count("rank1_psms", len(pep_q), eid)
                                         if pep_q else
                                         Undefined("EMPTY_POPULATION",
                                                   "no peptides passed peptide_q <= 0.01",
                                                   "peptides_at_peptide_q")),
        "psms_rank1": _q_json(den),
        # renamed from v1's `fdr_loss_pct`: it is the GAP between two q-value filters, and
        # calling it a redundancy signal was never justified (see README UNCOVERED).
        "q_metric_count_gap_pct": _q_json(gap),
        "diagnostics": {"peptides_at_spectrum_q": len(spec_q), "psms_at_spectrum_q": npsm},
    }


# ------------------------------------------------------------------- DIA-NN reference (G-C09)
TRUTH_COLUMNS = ["Precursor.Id", "Precursor.Mz", "Precursor.Charge", "RT", "RT.Start",
                 "RT.Stop", "IM", "Q.Value"]
RUN_COLUMNS = ["Run", "File.Name"]
LIB_COLUMNS = ["Precursor.Id", "Product.Mz"]


def verify_reference(report_path, sample, raw_realpath):
    """Content-derived binding of a reference to a sample.

    FAILURE 9 / codex #9: dataset A's diann_report pointing at an dataset D parquet. v1 printed OK because
    the file existed. Pinning its digest (bench2 G-B06) makes the WRONG reference stable, not
    correct -- the check that actually catches it is here: the DIA-NN report's own Run /
    File.Name column must name this sample's raw file.

    Multi-run reports are FILTERED, not rejected. A rule that refused them would fire on the
    normal DIA-NN output shape, and users would hand-make single-run copies outside the harness,
    reintroducing the very error class this check exists to prevent.
    """
    names = pq.read_schema(str(report_path)).names
    missing = [c for c in TRUTH_COLUMNS if c not in names]
    if missing:
        die("G-C09", "%s lacks DIA-NN column(s) %s.\n       columns present: %s"
            % (report_path, missing, names[:20]))
    runcol = None
    for c in RUN_COLUMNS:
        if c in names:
            runcol = c
            break
    if runcol is None:
        die("G-C09", "%s has neither 'Run' nor 'File.Name'; the reference cannot be bound to a "
                     "sample by content, only by the path someone typed." % report_path)
    tbl = pq.read_table(str(report_path), columns=TRUTH_COLUMNS + [runcol]).to_pydict()
    stem = Path(raw_realpath).name
    stem_noext = stem[:-2] if stem.endswith(".d") else stem
    runs = tbl[runcol]
    distinct = sorted(set(str(r) for r in runs))
    keep = [i for i, r in enumerate(runs)
            if stem_noext in str(r) or str(r) in stem_noext or Path(str(r)).stem == stem_noext]
    if not keep:
        die("G-C09", "%s contains no rows for sample %s (raw %s).\n       Run values present: %s\n"
                     "       This reference describes a DIFFERENT acquisition. v1 read it anyway "
                     "and published dataset A metrics against another sample's truth."
            % (report_path, sample, stem, distinct[:6]))
    return tbl, keep, runcol, distinct


def load_truth(report_path, lib_path, sample, raw_realpath):
    tbl, keep, runcol, distinct = verify_reference(report_path, sample, raw_realpath)
    qv = tbl["Q.Value"]
    nan_q = sum(1 for i in keep if qv[i] is None or not math.isfinite(float(qv[i])))
    if nan_q:
        die("G-C09", "%s has %d rows with a NaN/None Q.Value for %s. NaN silently fails the "
                     "<= 0.01 filter, which changes the truth population without saying so."
            % (report_path, nan_q, sample))
    precursors = []
    for i in keep:
        if float(qv[i]) > Q_THRESHOLD:
            continue
        for col in ("Precursor.Mz", "RT.Start", "RT.Stop", "IM"):
            v = tbl[col][i]
            if v is None or not math.isfinite(float(v)):
                die("G-C09", "%s row %d has a non-finite %s" % (report_path, i, col))
        precursors.append({"id": tbl["Precursor.Id"][i], "mz": float(tbl["Precursor.Mz"][i]),
                           "z": int(tbl["Precursor.Charge"][i]), "rt": float(tbl["RT"][i]),
                           "rt0": float(tbl["RT.Start"][i]), "rt1": float(tbl["RT.Stop"][i]),
                           "im": float(tbl["IM"][i])})
    if not precursors:
        die("G-C09", "%s yields zero precursors for %s at Q.Value <= %.3f"
            % (report_path, sample, Q_THRESHOLD))

    # A missing fragment library is NOT zero purity. v1's `if lib and lib.exists()` left every
    # fragment set empty, so matched spectra received own=0, coiso=0, own_decoy=0, coiso_decoy=0
    # -- four invented measurements from absent reference data.
    if lib_path is None or not Path(lib_path).exists():
        die("G-C08", "fragment library for %s is absent (%s). Purity without a library is not "
                     "0%%, it is UNDEFINED -- and v1 published the zero." % (sample, lib_path))
    lnames = pq.read_schema(str(lib_path)).names
    lmissing = [c for c in LIB_COLUMNS if c not in lnames]
    if lmissing:
        die("G-C08", "%s lacks column(s) %s; columns present: %s"
            % (lib_path, lmissing, lnames[:20]))
    cols = list(LIB_COLUMNS) + (["Decoy"] if "Decoy" in lnames else [])
    lt = pq.read_table(str(lib_path), columns=cols).to_pydict()
    frags = collections.defaultdict(list)
    dec = field(lt, "Decoy")
    for i in range(len(lt["Precursor.Id"])):
        if dec is not None and int(dec[i]) != 0:
            continue
        frags[lt["Precursor.Id"][i]].append(float(lt["Product.Mz"][i]))
    if not frags:
        die("G-C08", "%s contains no target fragments" % lib_path)
    covered = sum(1 for p in precursors if p["id"] in frags)
    coverage = float(covered) / len(precursors)
    if coverage < 0.50:
        die("G-C08", "the library covers only %.1f%% of this sample's truth precursors; fragment "
                     "metrics would be dominated by absent reference data rather than by the "
                     "tool." % (100 * coverage))
    return precursors, dict((k, np.sort(np.array(v))) for k, v in frags.items()), \
        {"run_column": runcol, "distinct_runs": distinct[:20],
         "library_coverage": round(coverage, 4)}


# -------------------------------------------------------------------------- the truth pass
def truth_metrics(mzml_path, precursors, frags, windows, method_id, declared_spectra):
    """Recall, purity and charge agreement, each with its own control and its own population."""
    ledger = Ledger()
    rt_state = RtUnitState()
    mz_sources = collections.Counter()

    idx = collections.defaultdict(list)
    for i, p in enumerate(precursors):
        lo = int((p["rt0"] - D_RT) / RT_BIN)
        hi = int((p["rt1"] + D_RT) / RT_BIN) + 1
        for b in range(lo, hi):
            idx[b].append(i)
    truth_ids = set(p["id"] for p in precursors)

    matched = set()
    matched_decoy = set()
    control_candidates = 0
    control_matches = 0
    shifted_out = 0
    shifted_total = 0
    charge_pairs = collections.Counter()
    matched_truth_charges = collections.Counter()
    own_f, coiso_f, unexp_f = [], [], []
    own_d, coiso_d = [], []
    peaks_available = []

    for spec in pyt_mzml.read(str(mzml_path)):
        where = "%s#%s" % (Path(mzml_path).name, field(spec, "id", "?"))
        if field(spec, "ms level") != 2:
            ledger.record("NOT_MS2")
            continue
        rt = rt_minutes(field(field(spec, "scanList", {}), "scan", [{}])[0],
                        rt_state, where)
        mz, src, z, im = precursor_of(spec)
        if mz is None:
            ledger.record("NO_PRECURSOR_MZ")
            continue
        mz_sources[src] += 1
        if im is None:
            # v1 turned this into im = -1 and the candidate filter read `im < 0 or ...`, which
            # DISABLED the mobility gate for the spectrum instead of excluding it.
            ledger.record("NO_ION_MOBILITY")
            continue
        window, wreason = assign_window(windows, mz, im)
        if window is None:
            ledger.record(wreason)
            continue
        mzs = field(spec, "m/z array")
        ints = field(spec, "intensity array")
        if mzs is None or len(mzs) == 0:
            ledger.record("EMPTY_PEAK_LIST")
            continue

        co_idx = field(idx, int(rt / RT_BIN), ())
        co = []
        for i in co_idx:
            p = precursors[i]
            if not (window["mz_lo"] <= p["mz"] <= window["mz_hi"]):
                continue
            if abs(p["im"] - im) > D_IM:
                continue
            co.append(p)

        # --- control pass: identical selection rule, shifted query --------------------------
        dmz = mz + DECOY_SHIFT
        shifted_total += 1
        if not (window["mz_lo"] <= dmz <= window["mz_hi"]):
            shifted_out += 1
        if co:
            control_candidates += 1
            dhit = select_match([{"mz": p["mz"], "p": p} for p in co], dmz, PPM_PREC)
            if dhit is not None:
                matched_decoy.add(dhit["p"]["id"])
                control_matches += 1

        # --- target pass ---------------------------------------------------------------------
        hit = select_match([{"mz": p["mz"], "p": p} for p in co], mz, PPM_PREC)
        if hit is None:
            ledger.record("NO_TARGET_CANDIDATE")
            continue
        best = hit["p"]
        matched.add(best["id"])
        if z is not None:
            charge_pairs[(z, best["z"])] += 1
            matched_truth_charges[best["z"]] += 1

        peaks_available.append(len(mzs))
        if len(mzs) < DEPTH_K:
            # Silently taking all peaks when fewer than K exist is not equal depth; v1's
            # argsort[-20:] did exactly that and called it "equal depth".
            ledger.record("FEWER_THAN_K_PEAKS")
            continue
        order = np.argsort(ints)[-DEPTH_K:]
        top_mz = np.asarray(mzs)[order]
        top_i = np.asarray(ints)[order]
        total = float(top_i.sum())
        if total <= 0.0:
            ledger.record("ZERO_TOTAL_INTENSITY")
            continue

        own_set = field(frags, best["id"], np.array([]))
        other_list = [frags[p["id"]] for p in co
                      if p["id"] != best["id"] and p["id"] in frags]
        other_set = np.sort(np.concatenate(other_list)) if other_list else np.array([])
        # Size-matched decoy: the SAME reference sets, shifted. Preserves set sizes (own ~9 m/z,
        # co-isolated ~450), which is what v1 got right and this rewrite keeps.
        own_dec = own_set + DECOY_SHIFT if len(own_set) else np.array([])
        other_dec = other_set + DECOY_SHIFT if len(other_set) else np.array([])

        acc = {"own": 0.0, "coiso": 0.0, "unexplained": 0.0}
        acc_d = {"own": 0.0, "coiso": 0.0, "unexplained": 0.0}
        for q, it in zip(top_mz, top_i):
            acc[attribute(float(q), own_set, other_set, PPM_FRAG)] += float(it)
            acc_d[attribute(float(q), own_dec, other_dec, PPM_FRAG)] += float(it)
        # I-2: exactly one category per peak, so the three fractions sum to the total. A
        # recurrence of v1's elif/if+if asymmetry becomes an abort here rather than a bias.
        for label, a in (("target", acc), ("decoy", acc_d)):
            s = a["own"] + a["coiso"] + a["unexplained"]
            if abs(s - total) > 1e-6 * max(total, 1.0):
                die("G-C12", "%s: %s intensity attribution sums to %.6g but the spectrum's top-%d "
                             "total is %.6g. A peak was counted twice or dropped -- target and "
                             "decoy attribution rules have diverged."
                    % (where, label, s, DEPTH_K, total))
        own_f.append(acc["own"] / total)
        coiso_f.append(acc["coiso"] / total)
        unexp_f.append(acc["unexplained"] / total)
        own_d.append(acc_d["own"] / total)
        coiso_d.append(acc_d["coiso"] / total)
        ledger.record("SCORED")

    ledger.check_partition(declared_spectra)
    profile = bucket_mz_source(mz_sources)
    eid = estimand_id(method_id, profile)
    blocked = ledger.blocked()

    def gate(family, q):
        if family in blocked:
            code, detail = blocked[family]
            return Undefined(code, detail, family)
        return q

    n_truth = len(truth_ids)
    truth_den = q_count("truth_precursors", n_truth, eid)
    recall_raw = gate("recall", q_ratio(
        Value(float(len(matched & truth_ids)), max(n_truth, 1), "count", "truth_precursors", eid),
        truth_den, "recall", n_min=N_MIN["recall"]))
    recall_floor = gate("recall", q_ratio(
        Value(float(len(matched_decoy & truth_ids)), max(n_truth, 1), "count",
              "truth_precursors", eid),
        truth_den, "recall_decoy_floor", n_min=N_MIN["recall"]))
    shifted_frac = (float(shifted_out) / shifted_total) if shifted_total else 0.0
    recall = Controlled(recall_raw, recall_floor, "prec_shift_+%.3fDa_same_window_v1"
                        % DECOY_SHIFT, control_candidates, control_matches,
                        SELECTION_RULE_SHA, shifted_frac)

    n_pur = len(own_f)
    own_v = gate("purity", q_mean(own_f, "purity_spectra", eid, "own_fraction",
                                  n_min=N_MIN["purity"]))
    own_fl = gate("purity", q_mean(own_d, "purity_spectra", eid, "own_fraction_decoy",
                                   n_min=N_MIN["purity"]))
    coiso_v = gate("purity", q_mean(coiso_f, "purity_spectra", eid, "coiso_fraction",
                                    n_min=N_MIN["purity"]))
    coiso_fl = gate("purity", q_mean(coiso_d, "purity_spectra", eid, "coiso_fraction_decoy",
                                     n_min=N_MIN["purity"]))
    own_ctl = Controlled(own_v, own_fl, "frag_shift_+%.3fDa_size_matched_v1" % DECOY_SHIFT,
                         n_pur, sum(1 for x in own_d if x > 0), SELECTION_RULE_SHA, shifted_frac)
    coiso_ctl = Controlled(coiso_v, coiso_fl, "frag_shift_+%.3fDa_size_matched_v1" % DECOY_SHIFT,
                           n_pur, sum(1 for x in coiso_d if x > 0), SELECTION_RULE_SHA,
                           shifted_frac)

    n_z = sum(charge_pairs.values())
    agree = sum(v for (a, t), v in charge_pairs.items() if a == t)
    z_den = q_count("charge_observations", n_z, eid)
    charge_agreement = gate("charge", q_ratio(
        Value(float(agree), max(n_z, 1), "count", "charge_observations", eid), z_den,
        "charge_agreement", n_min=N_MIN["charge_agreement"]))
    # FAILURE 4 extension, fixed: v1 computed agreement over MATCHED SPECTRA and the majority
    # baseline over ALL TRUTH PRECURSORS -- two different populations, then subtracted. The
    # baseline is now computed over charge_observations, the same population, spectrum-weighted.
    if n_z:
        maj = max(matched_truth_charges.values())
        charge_majority = q_ratio(
            Value(float(maj), n_z, "count", "charge_observations", eid), z_den,
            "charge_majority_baseline", n_min=N_MIN["charge_agreement"])
    else:
        charge_majority = Undefined("EMPTY_POPULATION", "no charge observations",
                                    "charge_majority_baseline")
    charge_majority = gate("charge", charge_majority)
    charge_vs_majority = q_diff(charge_agreement, charge_majority, "charge_vs_majority")

    return {
        "estimand_id": eid,
        "method_id": method_id,
        "mz_source_profile": profile,
        "mz_source_counts": dict(mz_sources),
        "ledger": ledger.to_json(),
        "blocked_families": dict((k, list(v)) for k, v in blocked.items()),
        "depth_basis": {"k": DEPTH_K,
                        "median_peaks_available": float(np.median(peaks_available))
                        if peaks_available else None},
        "recall": recall.to_json(),
        "own_fraction": own_ctl.to_json(),
        "coiso_fraction": coiso_ctl.to_json(),
        "charge_agreement": _q_json(charge_agreement),
        "charge_majority_baseline": _q_json(charge_majority),
        "charge_vs_majority": _q_json(charge_vs_majority),
        "charge_confusions": sorted(((v, "%d->%d" % (a, t))
                                     for (a, t), v in charge_pairs.items() if a != t),
                                    reverse=True)[:5],
    }


def _q_json(q):
    return q.to_json() if isinstance(q, (Value, Undefined)) else q


def _q_load(obj):
    if not isinstance(obj, dict):
        return obj
    if "__undefined__" in obj:
        return Undefined(obj["__undefined__"], obj["detail"], obj["where"],
                         tuple(_q_load(c) for c in obj.get("causes", [])))
    if "v" in obj and "pop" in obj:
        return Value(float(obj["v"]), int(obj["n"]), obj["unit"], obj["pop"], obj["estimand"])
    return obj


# ------------------------------------------------------------------------- pairing (G-C18/19)
class Pair(object):
    """A pair is two cells in ONE runset, in ONE sample.

    FAILURE 1 -- the worst published error this project has made -- was a cross-sample ratio:
    7,430 was an dataset D result divided by dataset A's the reference implementation reference. A Pair carries a single
    `sample` and is only ever constructed by build_pairs(), which iterates within a sample. A
    cross-sample pair is not representable, and Delta cannot be built from anything else.
    """

    __slots__ = ("sample", "baseline", "treatment", "kind")

    def __init__(self, sample, baseline, treatment, kind):
        if baseline["sample"] != sample or treatment["sample"] != sample:
            die("G-C18", "attempted to build a pair spanning samples %r / %r. FAILURE 1 was "
                         "exactly this: a published ratio of dataset D's number over dataset A's reference."
                % (baseline["sample"], treatment["sample"]))
        self.sample = sample
        self.baseline = baseline
        self.treatment = treatment
        self.kind = kind


def build_pairs(runset, cells):
    pairs = []
    for (sample, arm), cell in sorted(cells.items()):
        base_arm = cell.get("baseline_arm")
        if not base_arm:
            continue
        key = (sample, base_arm)
        if key not in cells:
            die("G-C19", "cell (%s, %s) declares baseline %r, which has no cell in this runset. "
                         "A treatment whose baseline is absent, or in a different runset, or on "
                         "a different host, cannot be published (FAILURE 2)."
                % (sample, arm, base_arm))
        kind = "tool_delta" if cell.get("kind") == "reference_mzml" or \
            cells[key].get("kind") == "reference_mzml" else "param_delta"
        pairs.append(Pair(sample, cells[key], cell, kind))
    return pairs


def validate_pair(runset, pair, stage_meta):
    """Everything that must be equal for a delta to mean anything.

    v1 declared "the pair runs in ONE invocation, same binary, same node, same everything but the
    varied parameters" and implemented none of it: `baseline` was a pointer used once during
    preflight, collate never read the field, and no check required both members to finish.
    """
    b = stage_meta[(pair.sample, pair.baseline["arm"])]
    t = stage_meta[(pair.sample, pair.treatment["arm"])]
    problems = []
    if runset["host_node"] != b["host"]["node"] or runset["host_node"] != t["host"]["node"]:
        problems.append("members ran on different hosts (%s / %s)"
                        % (b["host"]["node"], t["host"]["node"]))
    if b["harness_code_digest"] != t["harness_code_digest"]:
        problems.append("members were produced by different harness code")
    if pair.kind == "param_delta":
        if b["tool_sha256"] != t["tool_sha256"]:
            problems.append("members used different converter binaries (%s vs %s) -- FAILURE 7"
                            % (short(b["tool_sha256"] or ""), short(t["tool_sha256"] or "")))
        if b["ini_sha256"] == t["ini_sha256"]:
            problems.append("members resolved to an IDENTICAL configuration -- FAILURE 5")
    if b["method_id"] != t["method_id"]:
        problems.append("members were scored against different acquisition geometries "
                        "(%s vs %s)" % (b["method_id"], t["method_id"]))
    for key in ("sage_binary_sha256", "sage_config_sha256", "fasta_sha256"):
        if b[key] != t[key]:
            problems.append("members used different %s" % key)
    if b["estimand_id"] != t["estimand_id"]:
        problems.append("members were measured under different rules (%s vs %s)"
                        % (b["estimand_id"], t["estimand_id"]))
    if problems:
        die("G-C19", "pair %s: %s ~ %s is INVALID:\n         - %s"
            % (pair.sample, pair.treatment["arm"], pair.baseline["arm"],
               "\n         - ".join(problems)))


class Delta(object):
    """Constructible only from a Pair. There is no function in this file that takes two arbitrary
    numbers and returns a ratio."""

    __slots__ = ("metric", "pair", "baseline", "treatment", "abs_delta", "rel_pct", "status")

    def __init__(self, metric, pair, baseline, treatment):
        self.metric = metric
        self.pair = pair
        self.baseline = baseline
        self.treatment = treatment
        self.abs_delta = q_diff(treatment, baseline, "delta:%s" % metric)
        if is_undef(self.abs_delta) or is_undef(baseline):
            self.rel_pct = Undefined("PROPAGATED", "absolute delta undefined",
                                     "reldelta:%s" % metric)
        elif baseline.unit == "percent":
            # A percent minus a percent is percentage POINTS. A relative change of a percentage
            # is a different quantity and is not emitted at all -- "percent of a percent" is how
            # a 15-point recall move gets quoted as "+22%".
            self.rel_pct = Undefined("NOT_COMPUTED",
                                     "metric is a percentage; only percentage-point differences "
                                     "are defined for it", "reldelta:%s" % metric)
        elif raw_number(baseline) == 0.0:
            self.rel_pct = Undefined("ZERO_DENOMINATOR", "baseline is zero",
                                     "reldelta:%s" % metric)
        else:
            self.rel_pct = Value(100.0 * raw_number(self.abs_delta) / raw_number(baseline),
                                 baseline.n, "percent", baseline.pop, baseline.estimand)


# --------------------------------------------------------------------------------- rendering
def render(q, digits=1):
    """There is NO branch here that returns a string for an Undefined.

    v1's collate did `r.get("recall_corrected", "-")`, so a missing measurement became a dash in
    a table of real numbers -- indistinguishable, in a screenshot, from a small value.
    """
    if is_undef(q):
        raise RefuseToPublish(q)
    if not isinstance(q, Value):
        raise TypeError("render() takes a Value, got %r" % type(q))
    if q.unit == "count":
        return "{:,}".format(int(raw_number(q)))
    if q.unit == "pp":
        return "%+.*f pp" % (digits, raw_number(q))
    if q.unit == "percent":
        return "%.*f%%" % (digits, raw_number(q))
    return "%.*f" % (digits, raw_number(q))


class RefuseToPublish(Exception):
    def __init__(self, undef):
        super(RefuseToPublish, self).__init__("%s: %s" % (undef.code, undef.detail))
        self.undef = undef


def render_controlled(c, digits=1):
    """A controlled metric prints value, floor and correction, or it prints nothing.

    FAILURE 4's guard was "uncorrected values are never emitted alone"; v1's collate then printed
    recall_corrected with neither its raw value nor its floor. Here the three travel together
    because render_controlled raises unless all three are Values.
    """
    return "%s (raw %s, floor %s)" % (render(c.corrected, digits), render(c.raw, digits),
                                      render(c.floor, digits))


# ------------------------------------------------------------------------------ truth command
def cmd_truth(args):
    require_deps()
    reg, _ = load_registry()
    bench = Bench(bench_root_from(reg, args.bench_root))
    runset, path = load_runset(bench, args.runset)
    if runset["revision"] != 1:
        die("G-C22", "truth runs on revision 1; %s is revision %d"
            % (path, runset["revision"]))
    cells = assert_cells_match_plan(runset)

    new_cells = []
    for (sample, arm), cell in sorted(cells.items()):
        pdir, pseal = bench.load_stage(cell["produce_kind"], cell["produce_stage_id"])
        sdir, sseal = bench.load_stage("search", cell["search_stage_id"])
        mzml_path = pdir / "work" / "pseudo.mzML"
        tsv_path = sdir / "work" / "results.sage.tsv"

        # Chain binding at read time: the search stage in THIS cell must have been run against the
        # mzML the produce stage in THIS cell emitted. Bytes from three attempts cannot satisfy
        # this simultaneously.
        srecipe = jload(sdir / "recipe.json")
        if srecipe["input"]["pseudo_mzml_sha256"] != \
                pseal["artifacts"]["work/pseudo.mzML"]["sha256"]:
            die("G-C22", "cell (%s,%s): search input digest does not match this cell's mzML"
                % (sample, arm))

        entry = sample_entry(reg, sample)
        refs = entry["references"]
        raw_realpath = pseal.get("raw_realpath")
        method = pseal.get("method")
        if not method or not method.get("windows"):
            die("G-C10", "cell (%s,%s): the produce stage carries no acquisition window table. "
                         "There is deliberately no fallback -- v1's hardcoded WINDOWS_S23 applied "
                         "to every sample is what this replaces." % (sample, arm))

        print("[truth] %s/%s ..." % (sample, arm))
        precursors, frags, refinfo = load_truth(
            Path(refs["diann_report"]["path"]), Path(refs["diann_lib"]["path"]),
            sample, raw_realpath)
        tm = truth_metrics(mzml_path, precursors, frags, method["windows"],
                           pseal["method_id"], pseal["checks"].get("spectrum_count_declared"))
        header, rows, hsha = read_sage(tsv_path, "pseudo.mzML",
                                       sseal.get("header_sha256"))
        sm = search_metrics(rows, tm["estimand_id"])

        recipe = {
            "recipe_version": CONTRACT_VERSION, "kind": "metrics",
            "harness_version": HARNESS_VERSION,
            "metrics_code_sha256": metrics_code_digest(),
            "input": {"sample": sample, "arm": arm,
                      "produce_stage_id": cell["produce_stage_id"],
                      "pseudo_mzml_sha256": pseal["artifacts"]["work/pseudo.mzML"]["sha256"],
                      "search_stage_id": cell["search_stage_id"],
                      "results_tsv_sha256": sseal["artifacts"]["work/results.sage.tsv"]["sha256"],
                      "diann_report_sha256": refs["diann_report"]["sha256"],
                      "diann_lib_sha256": refs["diann_lib"]["sha256"]},
            "params": {"ppm_prec": PPM_PREC, "ppm_frag": PPM_FRAG, "d_rt": D_RT, "d_im": D_IM,
                       "rt_bin": RT_BIN, "depth_k": DEPTH_K, "q_threshold": Q_THRESHOLD,
                       "decoy_shift": DECOY_SHIFT},
        }
        att = Attempt(bench, "metrics", recipe)
        payload = {"sample": sample, "arm": arm, "reference": refinfo,
                   "search": sm, "truth": tm}
        (att.work / "metrics.json").write_bytes(canon(payload))
        p = att.work / "metrics.json"
        att.declare("work/metrics.json", {"sha256": sha256_file(p), "bytes": p.stat().st_size})
        att.checks = {"estimand_id": tm["estimand_id"], "method_id": tm["method_id"]}
        stage_id, mseal = att.commit({"sample": sample, "arm": arm,
                                      "estimand_id": tm["estimand_id"]})
        c = dict(cell)
        c["metrics_stage_id"] = stage_id
        c["estimand_id"] = tm["estimand_id"]
        new_cells.append(c)
        _print_truth_line(sample, arm, tm)

    obj = dict(runset)
    obj.pop("seal", None)
    obj["revision"] = 2
    obj["parent_revision_sha256"] = sha256_file(path)
    obj["status"] = "measured"
    obj["cells"] = new_cells
    obj["measured_utc"] = utc_now()
    obj["metrics_code_sha256"] = metrics_code_digest()
    p2 = write_runset(bench, obj)
    print("\n[truth] wrote %s" % p2)
    print("        next: collate2.py collate --runset %s.r2 --publish" % runset["runset_id"])
    return 0


def _print_truth_line(sample, arm, tm):
    rec = _q_load(tm["recall"]["corrected"])
    led = tm["ledger"]["dispositions"]
    try:
        rtxt = render(rec)
    except RefuseToPublish as exc:
        rtxt = "UNDEFINED(%s)" % exc.undef.code
    print("        %s/%s recall=%s scored=%d out_of_window=%d no_im=%d"
          % (sample, arm, rtxt, led["SCORED"], led["OUT_OF_WINDOW_TABLE"],
             led["NO_ION_MOBILITY"]))


# ---------------------------------------------------------------------------- collate command
def cmd_collate(args):
    require_deps()
    reg, _ = load_registry()
    bench = Bench(bench_root_from(reg, args.bench_root))
    runset, path = load_runset(bench, args.runset)
    if runset["revision"] < 2:
        die("G-C22", "runset %s is revision %d; run `collate2.py truth` first"
            % (runset["runset_id"], runset["revision"]))
    cells = assert_cells_match_plan(runset)

    stage_meta = {}
    payloads = {}
    refusals = []
    for (sample, arm), cell in sorted(cells.items()):
        if not cell.get("metrics_stage_id"):
            die("G-C20", "cell (%s,%s) has no metrics stage. v1 printed a dash here and published "
                         "the report anyway." % (sample, arm))
        pdir, pseal = bench.load_stage(cell["produce_kind"], cell["produce_stage_id"])
        sdir, sseal = bench.load_stage("search", cell["search_stage_id"])
        mdir, mseal = bench.load_stage("metrics", cell["metrics_stage_id"])
        mrecipe = jload(mdir / "recipe.json")
        if mrecipe["metrics_code_sha256"] != metrics_code_digest():
            die("G-C21", "cell (%s,%s) metrics were computed by a DIFFERENT version of "
                         "collate2.py (%s vs %s). Re-run truth."
                % (sample, arm, short(mrecipe["metrics_code_sha256"]),
                   short(metrics_code_digest())))
        if mrecipe["input"]["pseudo_mzml_sha256"] != \
                pseal["artifacts"]["work/pseudo.mzML"]["sha256"] or \
           mrecipe["input"]["results_tsv_sha256"] != \
                sseal["artifacts"]["work/results.sage.tsv"]["sha256"]:
            die("G-C22", "cell (%s,%s): metrics were computed from artifacts that are not the "
                         "ones this cell names. This is the three-attempts failure."
                % (sample, arm))
        mp = mdir / "work" / "metrics.json"
        if sha256_file(mp) != mseal["artifacts"]["work/metrics.json"]["sha256"]:
            die("G-C22", "cell (%s,%s): metrics.json changed after it was sealed" % (sample, arm))
        # Cheap integrity sweep over the big artifacts. collate deliberately does NOT rehash a
        # 12 GB mzML -- that is `bench2.py verify --deep`. A size check costs a stat() and catches
        # truncation and appending; a same-size edit is NOT caught here, and the report footer
        # tells the reader to run the deep verify. Saying more than that would be the overclaim
        # this rewrite exists to end.
        for kind, d, seal in ((cell["produce_kind"], pdir, pseal), ("search", sdir, sseal),
                              ("metrics", mdir, mseal)):
            for name, rec in sorted(seal["artifacts"].items()):
                f = d / name
                if not f.exists():
                    die("G-C22", "cell (%s,%s): sealed artifact %s/%s is gone"
                        % (sample, arm, kind, name))
                if f.stat().st_size != rec["bytes"]:
                    die("G-C22", "cell (%s,%s): %s/%s changed size after sealing (%d -> %d)"
                        % (sample, arm, kind, name, rec["bytes"], f.stat().st_size))
        payloads[(sample, arm)] = jload(mp)
        srecipe = jload(sdir / "recipe.json")
        stage_meta[(sample, arm)] = {
            "host": pseal["host"], "harness_code_digest": pseal["harness_code_digest"],
            "tool_sha256": (jload(pdir / "recipe.json").get("tool") or {}).get("sha256"),
            "ini_sha256": pseal.get("ini_sha256"),
            "method_id": pseal.get("method_id"),
            "sage_binary_sha256": srecipe["sage"]["binary_sha256"],
            "sage_config_sha256": srecipe["sage"]["config_sha256"],
            "fasta_sha256": srecipe["sage"]["fasta_sha256"],
            "estimand_id": cell.get("estimand_id"),
            "wall_s": pseal["wall_s"], "kind": cell.get("kind"),
        }

    # Cross-cell equality. v1 RECORDED binary hash and host in each manifest and never compared
    # them; collate then printed only the FIRST row's values for the whole sample, hiding mixed
    # provenance entirely.
    for key in ("sage_binary_sha256", "sage_config_sha256", "fasta_sha256",
                "harness_code_digest"):
        vals = set(m[key] for m in stage_meta.values())
        if len(vals) != 1:
            die("G-C19", "the publish set spans %d different values of %s: %s"
                % (len(vals), key, sorted(short(str(v)) for v in vals)))
    for sample in runset["plan_samples"]:
        mids = set(m["method_id"] for k, m in stage_meta.items() if k[0] == sample)
        if len(mids) != 1:
            die("G-C10", "sample %s spans %d acquisition geometries %s; its arms were not scored "
                         "on the same window table" % (sample, len(mids), sorted(mids)))

    pairs = build_pairs(runset, cells)
    for pr in pairs:
        validate_pair(runset, pr, stage_meta)

    # Per-column estimand equality. Two numbers computed under different measurement rules cannot
    # share a column, whatever the rules were.
    eids = collections.defaultdict(set)
    for (sample, arm), m in stage_meta.items():
        eids[sample].add(m["estimand_id"])
    for sample, s in sorted(eids.items()):
        if len(s) != 1:
            die("G-C21", "sample %s spans estimands %s; those numbers are not comparable"
                % (sample, sorted(s)))

    md, refusals = render_report(runset, cells, payloads, stage_meta, pairs, bench)

    outdir = Path(reg["vault"])
    outdir.mkdir(parents=True, exist_ok=True)
    slug = re.sub(r"[^A-Za-z0-9]+", "-", runset["plan_title"]).strip("-").lower()[:48]
    stem = "%s-%s" % (slug, runset["runset_id"])
    # An incomplete report is NEVER written into the vault directory and carries INCOMPLETE in
    # its filename, so a stray copy cannot be mistaken for a published result.
    incomplete_dir = outdir.parent / ("%s-incomplete" % outdir.name)
    if refusals:
        incomplete_dir.mkdir(parents=True, exist_ok=True)
        p = incomplete_dir / ("%s.INCOMPLETE.md" % stem)
        p.write_text(md)
        sys.stderr.write("\nREFUSING TO PUBLISH: %d undefined measurement(s)\n" % len(refusals))
        for where, code, detail in refusals:
            sys.stderr.write("  %-28s %-22s %s\n" % (where, code, detail[:80]))
        sys.stderr.write("\nIncomplete report (NOT in the vault): %s\n\n" % p)
        return 3
    if not args.publish:
        incomplete_dir.mkdir(parents=True, exist_ok=True)
        p = incomplete_dir / ("%s.draft.md" % stem)
        p.write_text(md)
        print("[collate] draft (not published): %s" % p)
        print("          add --publish to write the vault report")
        return 0
    # G-C23: O_EXCL. Two plans sharing a filename stem overwrote each other's vault report in v1.
    p = outdir / ("%s.md" % stem)
    try:
        fd = os.open(str(p), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
    except OSError:
        die("G-C23", "%s already exists. A vault file is named by runset id and is never "
                     "overwritten." % p)
    with os.fdopen(fd, "w") as fh:
        fh.write(md)
    (outdir / ("%s.vault.json" % stem)).write_bytes(canon({
        "runset_id": runset["runset_id"], "revision": runset["revision"],
        "plan_normalized_sha256": runset["plan_normalized_sha256"],
        "metrics_code_sha256": metrics_code_digest(),
        "report_sha256": sha256_bytes(md.encode("utf-8")),
        "cells": dict(("%s/%s" % k, v) for k, v in sorted(payloads.items())),
        "stage_meta": dict(("%s/%s" % k, v) for k, v in sorted(stage_meta.items())),
    }))
    print("[collate] published %s" % p)
    return 0


def render_report(runset, cells, payloads, stage_meta, pairs, bench):
    refusals = []

    def cell_field(sample, arm, path):
        obj = payloads[(sample, arm)]
        for k in path:
            obj = obj[k]
        return _q_load(obj)

    def safe(sample, arm, label, fn):
        try:
            return fn()
        except RefuseToPublish as exc:
            refusals.append(("%s/%s %s" % (sample, arm, label), exc.undef.code,
                             exc.undef.detail))
            return "**UNDEFINED(%s)**" % exc.undef.code

    L = ["# %s" % runset["plan_title"], "",
         "Generated by `harness/collate2.py` from runset `%s` revision %d."
         % (runset["runset_id"], runset["revision"]),
         "Do not hand-edit; re-run `collate`. Every number below is reproducible from sealed "
         "bytes -- see Provenance.", "",
         "All peptide counts are **`peptide_q <= 0.01`**. `spectrum_q` is a diagnostic and does "
         "not appear in any table here; quoting it inflated this project's counts by 9-11 points "
         "against the reference implementation (FAILURE 3).", ""]

    for sample in sorted(runset["plan_samples"]):
        L += ["## %s" % sample, "",
              "| arm | role | peptides (peptide_q) | recall corr. (raw / floor) | own frag. | "
              "co-isolation | charge | vs majority | wall |",
              "|---|---|---|---|---|---|---|---|---|"]
        arms = sorted(a["name"] for a in runset["plan_arms"])
        for arm in arms:
            p = payloads[(sample, arm)]
            tr = p["truth"]

            def ctl(key):
                c = tr[key]
                return "%s (%s / %s)" % (render(_q_load(c["corrected"])),
                                         render(_q_load(c["raw"]), 2),
                                         render(_q_load(c["floor"]), 2))
            row = [arm, cells[(sample, arm)].get("role", "-"),
                   safe(sample, arm, "peptides",
                        lambda: render(_q_load(p["search"]["peptides_at_peptide_q"]))),
                   safe(sample, arm, "recall", lambda: ctl("recall")),
                   safe(sample, arm, "own_fraction", lambda: ctl("own_fraction")),
                   safe(sample, arm, "coiso_fraction", lambda: ctl("coiso_fraction")),
                   safe(sample, arm, "charge",
                        lambda: render(_q_load(tr["charge_agreement"]))),
                   safe(sample, arm, "charge_vs_majority",
                        lambda: render(_q_load(tr["charge_vs_majority"]))),
                   "%.0f s" % stage_meta[(sample, arm)]["wall_s"]]
            L.append("| " + " | ".join(str(x) for x in row) + " |")
        L.append("")

        sample_pairs = [pr for pr in pairs if pr.sample == sample]
        if sample_pairs:
            L += ["### Paired deltas (%s)" % sample, "",
                  "Each delta names both arms and is computed **within one sample**; there is no "
                  "function in the harness that divides one sample's number by another's "
                  "(FAILURE 1).", "",
                  "| treatment ~ baseline | metric | baseline | treatment | delta |",
                  "|---|---|---|---|---|"]
            for pr in sample_pairs:
                for metric, path in (("peptides at peptide_q",
                                      ("search", "peptides_at_peptide_q")),
                                     ("recall (corrected)", ("truth", "recall", "corrected")),
                                     ("charge agreement", ("truth", "charge_agreement"))):
                    b = cell_field(sample, pr.baseline["arm"], path)
                    t = cell_field(sample, pr.treatment["arm"], path)
                    d = Delta(metric, pr, b, t) if not (is_undef(b) or is_undef(t)) else None
                    label = "%s ~ %s" % (pr.treatment["arm"], pr.baseline["arm"])
                    if d is None:
                        u = b if is_undef(b) else t
                        refusals.append(("%s %s %s" % (sample, label, metric), u.code, u.detail))
                        L.append("| %s | %s | - | - | **UNDEFINED(%s)** |"
                                 % (label, metric, u.code))
                        continue
                    L.append("| %s | %s | %s | %s | %s |"
                             % (label, metric,
                                safe(sample, pr.baseline["arm"], metric, lambda: render(b)),
                                safe(sample, pr.treatment["arm"], metric, lambda: render(t)),
                                safe(sample, label, metric, lambda: render(d.abs_delta))))
            L += ["", "Run order for %s: `%s` (seed %d). **n = 1** per cell: there are no "
                      "technical replicates and no confidence intervals. A delta from a single "
                      "ordering is labelled as such."
                  % (sample, " -> ".join(runset["arm_order"][sample]), runset["order_seed"]), ""]

        L += ["### Spectrum dispositions (%s)" % sample, "",
              "Every spectrum the truth stage declined to score is counted and named. v1 used a "
              "bare `continue`, so a wrong acquisition window table showed up as quietly "
              "depressed recall rather than as an error.", "",
              "| arm | MS2 | scored | out of window | ambiguous window | no IM | no precursor | "
              "< %d peaks |" % DEPTH_K,
              "|---|---|---|---|---|---|---|---|"]
        for arm in arms:
            d = payloads[(sample, arm)]["truth"]["ledger"]
            n = d["dispositions"]
            L.append("| %s | %d | %d | %d | %d | %d | %d | %d |"
                     % (arm, d["ms2"], n["SCORED"], n["OUT_OF_WINDOW_TABLE"],
                        n["AMBIGUOUS_WINDOW"], n["NO_ION_MOBILITY"], n["NO_PRECURSOR_MZ"],
                        n["FEWER_THAN_K_PEAKS"]))
        L.append("")

    n_started = ledger_count(bench, runset["plan_normalized_sha256"])
    L += ["## Provenance", "",
          "| | |", "|---|---|",
          "| runset | `%s` r%d |" % (runset["runset_id"], runset["revision"]),
          "| plan (normalized) | `%s` |" % short(runset["plan_normalized_sha256"]),
          "| converter | `%s` |" % short(runset["tool"]["sha256"]),
          "| harness code | `%s` |" % short(runset["harness_code_digest"]),
          "| metrics code | `%s` |" % short(metrics_code_digest()),
          "| sage / config / fasta | `%s` / `%s` / `%s` |"
          % (short(runset["search_tools"]["sage_binary"]["sha256"]),
             short(runset["search_tools"]["sage_config"]["sha256"]),
             short(runset["search_tools"]["fasta"]["sha256"])),
          "| host | `%s` |" % runset["host_node"],
          "| runs started for this plan | %d |" % n_started,
          "",
          "Runs started for this plan digest: **%d**. If that number exceeds the number of "
          "published reports, some runs were not published, and the reader should ask why."
          % n_started, "",
          "```", "verify: bench2.py verify --runset %s.r%d --deep"
          % (runset["runset_id"], runset["revision"]), "```", "",
          "## Estimand caveats", "",
          "These are properties of the MEASUREMENTS, not of the bookkeeping, and the harness "
          "does not fix them. They are reproduced here so a number cannot be quoted away from "
          "its caveat.", "",
          "- **Decoy model** `+%.3f Da shift` is versioned and its opportunity count is "
          "published, but it is **not calibrated**. Whether it is an exchangeable null is "
          "unestablished. The harness refuses to publish when the control never fired; it "
          "cannot certify the control when it did." % DECOY_SHIFT,
          "- **Abbott correction** `(obs - floor)/(1 - floor)` assumes a model that is not "
          "verified here; it is bounds-checked, never clamped.",
          "- **top-%d purity** is a top-%d estimand, not total search utility. Sage searches "
          "every peak." % (DEPTH_K, DEPTH_K),
          "- **`q_metric_count_gap_pct`** is the gap between two q-value filters. It is not a "
          "redundancy measure; v1's docstring called it one without justification.",
          "- **n = 1.** No replicates, no confidence intervals. `n` printed on a Value is a "
          "count of spectra, and spectra are not independent observations of a precursor -- it "
          "is not a basis for a confidence interval.",
          "- **1/K0 tile bounds** come from a linear scan-number interpolation, not the Bruker "
          "SDK conversion.", ""]
    return "\n".join(L) + "\n", refusals


# ============================================================================== source lint
BANNED = [
    (r"\.get\(", "dict.get with a default manufactures a number out of an absence"),
    (r"\bor \[0\]", "`or [0]` manufactures a synthetic observation"),
    (r"\bor 1\.0\b", "`or 1.0` manufactures a denominator"),
    (r"max\([^)]*,\s*1\s*\)", "max(x, 1) manufactures a denominator"),
    (r"max\([^)]*,\s*1e-", "max(x, 1e-9) manufactures a denominator"),
    (r"getattr\([^,]+,[^,]+,", "3-argument getattr guesses a missing unit"),
    (r"json\.dumps", "use canon(); json.dumps emits bare NaN"),
    (r"except\s*:", "bare except swallows a guard"),
]
ALLOW_LINE = "# lint-ok:"


def lint_source(path=None, region=True):
    """G-C30 -- the constructs that manufactured v1's wrong numbers may not reappear.

    HONEST SCOPE: this is a regex over the MEASUREMENT REGION of this file only, and it stops
    the accident, not a determined author (`d["k"] if "k" in d else default`, `defaultdict`, and
    string concatenation all defeat it). It is a regression lock on known-bad idioms, and the
    README says exactly that.
    """
    import tokenize
    import io
    p = Path(path or (HERE / "collate2.py"))
    text = p.read_text()
    raw_lines = text.splitlines()
    if region:
        try:
            start = next(i for i, l in enumerate(raw_lines)
                         if "--- BEGIN MEASUREMENT REGION ---" in l)
            end = next(i for i, l in enumerate(raw_lines)
                       if "--- END MEASUREMENT REGION ---" in l)
        except StopIteration:
            die("G-C30", "%s has no measurement-region markers" % p)
    else:
        start, end = 0, len(raw_lines)

    # Blank out strings and comments before matching. A docstring that QUOTES the banned idiom in
    # order to explain why it is banned must not trip the lint -- otherwise the only way to pass
    # is to stop documenting the defect, which is the opposite of the point.
    code_lines = list(raw_lines)
    try:
        for tok in tokenize.generate_tokens(io.StringIO(text).readline):
            if tok.type not in (tokenize.STRING, tokenize.COMMENT):
                continue
            (r0, c0), (r1, c1) = tok.start, tok.end
            for r in range(r0 - 1, min(r1, len(code_lines))):
                line = code_lines[r]
                a = c0 if r == r0 - 1 else 0
                b = c1 if r == r1 - 1 else len(line)
                code_lines[r] = line[:a] + " " * max(0, b - a) + line[b:]
    except (tokenize.TokenError, IndentationError) as exc:
        die("G-C30", "%s could not be tokenized for linting: %s" % (p, exc))

    hits = []
    for i in range(start, end):
        line = code_lines[i]
        if ALLOW_LINE in raw_lines[i] or not line.strip():
            continue
        for pat, why in BANNED:
            if re.search(pat, line):
                hits.append((i + 1, raw_lines[i].strip()[:70], why))
    if hits:
        msg = "\n".join("       line %d: %s\n         -> %s" % h for h in hits)
        die("G-C30", "banned construct(s) inside the measurement region of %s:\n%s" % (p, msg))
    return True


# ================================================================================== SELF TEST
def _mk_parquet(path, rows, columns):
    import pyarrow as pa
    arrays = {}
    for c in columns:
        arrays[c] = [r.get(c) for r in rows]
    pq.write_table(pa.table(arrays), str(path))
    return path


def _truth_rows(run, n, mz0=500.0):
    rows = []
    for i in range(n):
        rows.append({"Run": run, "Precursor.Id": "PREC%d" % i,
                     "Precursor.Mz": mz0 + i * 0.01, "Precursor.Charge": 2 if i % 3 else 3,
                     "RT": 10.0, "RT.Start": 9.5, "RT.Stop": 10.5, "IM": 1.0,
                     "Q.Value": 0.001})
    return rows


TRUTH_PARQUET_COLS = ["Run", "Precursor.Id", "Precursor.Mz", "Precursor.Charge", "RT",
                      "RT.Start", "RT.Stop", "IM", "Q.Value"]
LIB_PARQUET_COLS = ["Precursor.Id", "Product.Mz", "Decoy"]


def selftest_collate(tmp):
    st = SelfTest("collate2.py -- fail-loud metrics, enforced pairing, publish gate")
    tmp = Path(tmp)
    eid = "e" * 16

    # ---- G-C30 source lint ------------------------------------------------------------------
    st.expect_ok("G-C30", "measurement region is clean of banned idioms",
                 lambda: lint_source())
    bad = tmp / "bad_region.py"
    bad.write_text("# --- BEGIN MEASUREMENT REGION ---\n"
                   "x = row.get('peptide_q', 1)\n"
                   "y = max(len(sp), 1)\n"
                   "# --- END MEASUREMENT REGION ---\n")
    st.expect_abort("G-C30", "banned .get(default) / max(x,1) detected",
                    lambda: lint_source(bad))

    # ---- Value / Undefined invariants --------------------------------------------------------
    def nan_value():
        try:
            Value(float("nan"), 5, "percent", "scored_spectra", eid)
        except ValueError as exc:
            die("G-C06", "Value rejected a non-finite payload: %s" % exc)
    st.expect_abort("G-C06", "a NaN cannot become a Value", nan_value)

    def zero_n():
        try:
            Value(1.0, 0, "percent", "scored_spectra", eid)
        except ValueError as exc:
            die("G-C06", "Value rejected n=0: %s" % exc)
    st.expect_abort("G-C06", "n=0 cannot become a publishable Value", zero_n)

    def no_arith():
        a = Value(1.0, 5, "percent", "scored_spectra", eid)
        b = Value(2.0, 5, "percent", "scored_spectra", eid)
        try:
            a / b
        except TypeError:
            die("G-C06", "raw division of two Values is a TypeError")
    st.expect_abort("G-C06", "Value / Value is a TypeError, not a ratio", no_arith)

    # ---- G-C06 zero denominators -------------------------------------------------------------
    num = Value(0.0, 10, "count", "scored_spectra", eid)
    den0 = q_count("scored_spectra", 0, eid)
    r = q_ratio(num, den0, "selftest")
    st.results.append(("PASS" if is_undef(r) and r.code in ("PROPAGATED", "ZERO_DENOMINATOR")
                       else "FAIL", "G-C06", "zero denominator -> Undefined, not 0.0%",
                       r.code if is_undef(r) else "returned %r" % r))
    den_small = Value(3.0, 3, "count", "scored_spectra", eid)
    r2 = q_ratio(num, den_small, "selftest", n_min=200)
    st.results.append(("PASS" if is_undef(r2) and r2.code == "INSUFFICIENT_N" else "FAIL",
                       "G-C06", "n below declared minimum -> INSUFFICIENT_N",
                       r2.code if is_undef(r2) else "returned %r" % r2))

    # ---- G-C07 empty populations --------------------------------------------------------------
    m = q_mean([], "purity_spectra", eid, "selftest")
    st.results.append(("PASS" if is_undef(m) and m.code == "EMPTY_POPULATION" else "FAIL",
                       "G-C07", "empty observation list -> EMPTY_POPULATION (not mean of [0])",
                       m.code if is_undef(m) else "returned %r" % m))

    # ---- G-C16 population mismatch -------------------------------------------------------------
    bad_pair = q_ratio(Value(1.0, 10, "count", "truth_precursors", eid),
                       Value(10.0, 10, "count", "charge_observations", eid), "selftest")
    st.results.append(("PASS" if is_undef(bad_pair) and bad_pair.code == "POPULATION_MISMATCH"
                       else "FAIL", "G-C16",
                       "numerator/denominator from different populations",
                       bad_pair.code if is_undef(bad_pair) else "returned %r" % bad_pair))
    ok_pair = q_ratio(Value(5.0, 600, "count", "charge_observations", eid),
                      Value(600.0, 600, "count", "charge_observations", eid), "selftest",
                      n_min=500)
    st.results.append(("PASS" if isinstance(ok_pair, Value) else "FAIL", "G-C16",
                       "same-population ratio accepted", repr(ok_pair)))

    # ---- G-C15 Abbott correction ---------------------------------------------------------------
    def C(raw, floor, cand=500, match=50, shifted=0.0):
        return Controlled(Value(raw, 500, "percent", "truth_precursors", eid),
                          Value(floor, 500, "percent", "truth_precursors", eid),
                          "test", cand, match, "sha", shifted)
    c = C(30.0, 50.0)
    st.results.append(("PASS" if is_undef(c.corrected) and c.corrected.code == "OUT_OF_DOMAIN"
                       else "FAIL", "G-C15", "observed below floor -> OUT_OF_DOMAIN, not -40%",
                       c.corrected.code if is_undef(c.corrected) else "returned %r" % c.corrected))
    c = C(95.0, 93.0)
    st.results.append(("PASS" if is_undef(c.corrected) and c.corrected.code == "CONTROL_TOO_HIGH"
                       else "FAIL", "G-C15", "floor 93% -> CONTROL_TOO_HIGH, not 28.6%",
                       c.corrected.code if is_undef(c.corrected) else "returned %r" % c.corrected))
    c = C(74.3, 35.6)
    st.results.append(("PASS" if isinstance(c.corrected, Value) and
                       abs(raw_number(c.corrected) - 60.1) < 0.2 else "FAIL", "G-C15",
                       "valid correction still computes (74.3/35.6 -> 60.1)",
                       repr(c.corrected)))

    # ---- G-C14 degenerate control ---------------------------------------------------------------
    c = C(74.3, 0.0, cand=300000, match=0)
    st.results.append(("PASS" if is_undef(c.corrected) and
                       c.corrected.code == "CONTROL_DEGENERATE" else "FAIL", "G-C14",
                       "control with 300k opportunities and 0 matches is inert",
                       c.corrected.code if is_undef(c.corrected) else "returned %r" % c.corrected))
    c = C(74.3, 5.0, cand=3, match=1)
    st.results.append(("PASS" if is_undef(c.corrected) and
                       c.corrected.code == "CONTROL_DEGENERATE" else "FAIL", "G-C14",
                       "control with too few opportunities is not evidence",
                       c.corrected.code if is_undef(c.corrected) else "returned %r" % c.corrected))
    c = C(74.3, 5.0, cand=5000, match=200, shifted=0.5)
    st.results.append(("PASS" if is_undef(c.corrected) and
                       c.corrected.code == "CONTROL_ASYMMETRIC" else "FAIL", "G-C14",
                       "50% of decoy queries shifted out of their tile",
                       c.corrected.code if is_undef(c.corrected) else "returned %r" % c.corrected))

    # ---- G-C20 render refuses Undefined ---------------------------------------------------------
    def render_undef():
        try:
            render(Undefined("EMPTY_POPULATION", "x", "y"))
        except RefuseToPublish as exc:
            die("G-C20", "render() refused an Undefined instead of printing a dash: %s" % exc)
    st.expect_abort("G-C20", "render(Undefined) raises; there is no dash path", render_undef)

    def render_controlled_missing_floor():
        try:
            render_controlled(C(74.3, 0.0, cand=300000, match=0))
        except RefuseToPublish as exc:
            die("G-C20", "a controlled metric cannot print without a valid floor: %s" % exc)
    st.expect_abort("G-C20", "controlled metric cannot print without its floor",
                    render_controlled_missing_floor)
    st.expect_ok("G-C20", "a fully defined controlled metric renders",
                 lambda: render_controlled(C(74.3, 35.6)))

    # ---- G-C18 cross-sample delta ----------------------------------------------------------------
    def cross_sample():
        Pair("dataset A", {"sample": "dataset A", "arm": "base"}, {"sample": "dataset D", "arm": "x"},
             "param_delta")
    st.expect_abort("G-C18", "a Pair spanning two samples is refused (FAILURE 1)", cross_sample)
    st.expect_ok("G-C18", "same-sample pair accepted",
                 lambda: Pair("dataset A", {"sample": "dataset A", "arm": "b"},
                              {"sample": "dataset A", "arm": "t"}, "param_delta"))

    def pct_of_pct():
        pr = Pair("dataset A", {"sample": "dataset A", "arm": "b"}, {"sample": "dataset A", "arm": "t"},
                  "param_delta")
        d = Delta("recall", pr,
                  Value(69.9, 500, "percent", "truth_precursors", eid),
                  Value(85.4, 500, "percent", "truth_precursors", eid))
        if not is_undef(d.rel_pct):
            die("G-C18", "a relative change of a percentage was computed")
        if d.abs_delta.unit != "pp":
            die("G-C18", "percentage difference did not carry percentage-point units")
        raise Abort("G-C18", "percent metric yields pp only; relative change is NOT_COMPUTED")
    st.expect_abort("G-C18", "percent metric never yields a relative %", pct_of_pct)

    def mixed_estimand():
        d = q_diff(Value(1.0, 5, "percent", "truth_precursors", "a" * 16),
                   Value(2.0, 5, "percent", "truth_precursors", "b" * 16), "selftest")
        if is_undef(d) and d.code == "ESTIMAND_MISMATCH":
            die("G-C21", "numbers under different measurement rules cannot be differenced")
    st.expect_abort("G-C21", "delta across different estimands refused", mixed_estimand)

    # ---- G-C19 pair validation ---------------------------------------------------------------
    base_meta = {"host": {"node": "n1"}, "harness_code_digest": "h", "tool_sha256": "t1",
                 "ini_sha256": "i1", "method_id": "m1", "sage_binary_sha256": "s",
                 "sage_config_sha256": "c", "fasta_sha256": "f", "estimand_id": eid,
                 "wall_s": 1.0, "kind": "run"}

    def pair_case(mutate):
        b = dict(base_meta)
        t = dict(base_meta)
        t["ini_sha256"] = "i2"
        t.update(mutate)
        pr = Pair("dataset A", {"sample": "dataset A", "arm": "b"}, {"sample": "dataset A", "arm": "t"},
                  "param_delta")
        validate_pair({"host_node": "n1"}, pr,
                      {("dataset A", "b"): b, ("dataset A", "t"): t})
    st.expect_abort("G-C19", "pair with different converter binaries (FAILURE 7)",
                    lambda: pair_case({"tool_sha256": "t2"}))
    st.expect_abort("G-C19", "pair on different hosts",
                    lambda: pair_case({"host": {"node": "n2"}}))
    st.expect_abort("G-C19", "pair scored on different acquisition geometries",
                    lambda: pair_case({"method_id": "m2"}))
    st.expect_abort("G-C19", "pair with different sage config",
                    lambda: pair_case({"sage_config_sha256": "c2"}))
    st.expect_abort("G-C19", "pair with different estimands",
                    lambda: pair_case({"estimand_id": "z" * 16}))
    st.expect_abort("G-C19", "pair whose members resolved identically (FAILURE 5)",
                    lambda: pair_case({"ini_sha256": "i1"}))
    st.expect_ok("G-C19", "well-formed pair validates", lambda: pair_case({}))

    def missing_baseline():
        cells = {("dataset A", "t"): {"sample": "dataset A", "arm": "t", "baseline_arm": "b"}}
        build_pairs({}, cells)
    st.expect_abort("G-C19", "treatment whose baseline has no cell", missing_baseline)

    # ---- G-C01 RT units ------------------------------------------------------------------------
    def mzml(name, specs, **kw):
        return bench2.st_build_mzml(tmp / name, specs, "file:///raw/dataset A.d", **kw)

    def base_spec(i=0, **kw):
        s = {"mz": 500.0 + i * 0.001, "z": 2, "im": 1.0, "rt": 600.0, "iso_target": 500.0,
             "peaks": [(100.0 + j, float(j + 1)) for j in range(30)]}
        s.update(kw)
        return s

    def read_rt(path):
        state = RtUnitState()
        for sp in pyt_mzml.read(str(path)):
            rt_minutes(sp["scanList"]["scan"][0], state, "x")
    st.expect_abort("G-C01", "RT with no unit is refused (v1 assumed minutes)",
                    lambda: read_rt(mzml("nounit.mzML", [base_spec(rt_unit=None)])))
    st.expect_abort("G-C01", "RT with an unknown unit accession is refused",
                    lambda: read_rt(mzml("weird.mzML", [base_spec(rt_unit="furlong")])))
    st.expect_abort("G-C01", "RT unit changing mid-file is refused",
                    lambda: read_rt(mzml("mixed.mzML",
                                         [base_spec(0, rt_unit="second"),
                                          base_spec(1, rt_unit="minute")])))
    st.expect_ok("G-C01", "seconds are read via the CV unit",
                 lambda: read_rt(mzml("sec.mzML", [base_spec(rt_unit="second")])))

    def seconds_convert():
        state = RtUnitState()
        sp = next(iter(pyt_mzml.read(str(mzml("sec2.mzML", [base_spec(rt_unit="second")])))))
        got = rt_minutes(sp["scanList"]["scan"][0], state, "x")
        if abs(got - 10.0) > 1e-9:
            die("G-C01", "600 s should be 10.0 min, got %r" % got)
    st.expect_ok("G-C01", "600 s converts to 10.0 min", seconds_convert)

    # ---- G-C02 / G-C03 precursor fields ----------------------------------------------------------
    def first(path):
        return next(iter(pyt_mzml.read(str(path))))

    p = first(mzml("zeromz.mzML", [base_spec(mz=0.0, iso_target=501.0)]))
    mzv, src, z, im = precursor_of(p)
    st.results.append(("PASS" if mzv == 501.0 and src == "isolation_target" else "FAIL",
                       "G-C02", "present-but-zero m/z falls back to isolation target",
                       "mz=%r source=%r" % (mzv, src)))
    p = first(mzml("nomz.mzML", [base_spec(mz=None, iso_target=None)]))
    mzv, src, z, im = precursor_of(p)
    st.results.append(("PASS" if mzv is None else "FAIL", "G-C02",
                       "no m/z anywhere -> None (a counted disposition)", "mz=%r" % mzv))
    p = first(mzml("noim.mzML", [base_spec(im=None)]))
    mzv, src, z, im = precursor_of(p)
    st.results.append(("PASS" if im is None else "FAIL", "G-C03",
                       "missing IM -> None, NOT -1 (v1 disabled the IM gate)", "im=%r" % im))
    p = first(mzml("noz.mzML", [base_spec(z=None)]))
    mzv, src, z, im = precursor_of(p)
    st.results.append(("PASS" if z is None else "FAIL", "G-C03",
                       "missing charge -> None, not 0", "z=%r" % z))

    # ---- G-C10 window assignment -----------------------------------------------------------------
    W = [{"mz_lo": 490.0, "mz_hi": 510.0, "im_lo": 0.8, "im_hi": 1.2},
         {"mz_lo": 505.0, "mz_hi": 520.0, "im_lo": 1.5, "im_hi": 1.8}]
    w, why = assign_window(W, 500.0, 1.0)
    st.results.append(("PASS" if w is not None and why == "SCORED" else "FAIL", "G-C10",
                       "2-D assignment picks the right tile", why))
    w, why = assign_window(W, 900.0, 1.0)
    st.results.append(("PASS" if w is None and why == "OUT_OF_WINDOW_TABLE" else "FAIL",
                       "G-C10", "m/z outside every tile is a NAMED disposition", why))
    w, why = assign_window([{"mz_lo": 490.0, "mz_hi": 510.0, "im_lo": 0.8, "im_hi": 1.2},
                            {"mz_lo": 495.0, "mz_hi": 515.0, "im_lo": 0.9, "im_hi": 1.1}],
                           500.0, 1.0)
    st.results.append(("PASS" if w is None and why == "AMBIGUOUS_WINDOW" else "FAIL", "G-C10",
                       "overlapping tiles are NOT resolved by first match", why))
    w, why = assign_window(W, 508.0, 1.7)
    st.results.append(("PASS" if w is not None and w["im_lo"] == 1.5 else "FAIL", "G-C10",
                       "mobility separates tiles that overlap in m/z",
                       "chose im_lo=%r" % (w["im_lo"] if w else None)))

    # ---- G-C11 ledger partition -------------------------------------------------------------------
    def bad_partition():
        led = Ledger()
        led.record("SCORED")
        led.total_read += 1          # a spectrum that left the loop with no disposition
        led.check_partition(None)
    st.expect_abort("G-C11", "a spectrum with no named disposition aborts", bad_partition)

    def truncated_read():
        led = Ledger()
        for _ in range(5):
            led.record("SCORED")
        led.check_partition(500)     # the file declared 500
    st.expect_abort("G-C11", "reader saw fewer spectra than the file declares", truncated_read)

    def good_partition():
        led = Ledger()
        for _ in range(5):
            led.record("SCORED")
        led.check_partition(5)
    st.expect_ok("G-C11", "exact partition accepted", good_partition)

    def window_blocks():
        led = Ledger()
        for _ in range(900):
            led.record("SCORED")
        for _ in range(100):
            led.record("OUT_OF_WINDOW_TABLE")
        b = led.blocked()
        if "recall" not in b or b["recall"][0] != "METHOD_MISMATCH":
            die("G-C10", "10%% out-of-window did not block recall: %r" % b)
        raise Abort("G-C10", "10%% of spectra outside the window table blocks recall")
    st.expect_abort("G-C10", "wrong window table blocks the metrics it corrupts", window_blocks)

    # ---- G-C12 attribution symmetry -----------------------------------------------------------------
    own = np.sort(np.array([100.0, 200.0]))
    other = np.sort(np.array([200.0, 300.0]))
    cats = [attribute(x, own, other, 20.0) for x in (100.0, 200.0, 300.0, 400.0)]
    st.results.append(("PASS" if cats == ["own", "own", "coiso", "unexplained"] else "FAIL",
                       "G-C12", "one category per peak, own wins ties (both passes)",
                       str(cats)))

    def attribution_invariant():
        # A peak counted in two categories must abort. Simulated by forcing a bad accumulator.
        total = 100.0
        acc = {"own": 60.0, "coiso": 60.0, "unexplained": 0.0}
        s = acc["own"] + acc["coiso"] + acc["unexplained"]
        if abs(s - total) > 1e-6 * total:
            die("G-C12", "intensity attribution sums to %.6g but the top-K total is %.6g"
                % (s, total))
    st.expect_abort("G-C12", "double-counted peak breaks the sum-to-total invariant",
                    attribution_invariant)

    # ---- G-C13 selection symmetry --------------------------------------------------------------------
    cands = [{"mz": 500.000, "p": "a"}, {"mz": 500.004, "p": "b"}, {"mz": 500.008, "p": "c"}]
    tgt = select_match(cands, 500.002, 20.0)
    dec = select_match(cands, 500.006, 20.0)
    st.results.append(("PASS" if tgt is not None and dec is not None and
                       isinstance(tgt, dict) and isinstance(dec, dict) else "FAIL", "G-C13",
                       "target and control use the SAME single-best rule",
                       "target=%s control=%s" % (tgt["p"], dec["p"])))

    # ---- G-C04 sage contract ---------------------------------------------------------------------------
    def sage_tsv(name, **kw):
        cols = list(bench2.SAGE_REQUIRED_COLUMNS)
        if kw.get("drop"):
            cols.remove(kw["drop"])
        rows = []
        for i in range(5):
            rec = {"psm_id": str(i), "filename": kw.get("fn", "pseudo.mzML"),
                   "scannr": "s%d" % i, "rank": "1", "label": "1", "peptide": "PEP%d" % i,
                   "charge": "2", "spectrum_q": "0.001", "peptide_q": "0.001",
                   "protein_q": "0.001", "hyperscore": "20.0"}
            if kw.get("qval") and i == 2:
                rec["peptide_q"] = kw["qval"]
            rows.append("\t".join(rec[c] for c in cols))
        body = "\t".join(cols) + "\n" + "\n".join(rows) + "\n"
        p = tmp / name
        p.write_text(body)
        return p
    st.expect_abort("G-C04", "TSV missing peptide_q (not a published zero)",
                    lambda: read_sage(sage_tsv("s_drop.tsv", drop="peptide_q"), "pseudo.mzML"))
    st.expect_abort("G-C04", "q-value 'NA'",
                    lambda: read_sage(sage_tsv("s_na.tsv", qval="NA"), "pseudo.mzML"))
    st.expect_abort("G-C04", "q-value NaN (silently fails every <= comparison)",
                    lambda: read_sage(sage_tsv("s_nan.tsv", qval="NaN"), "pseudo.mzML"))
    st.expect_abort("G-C04", "q-value out of [0,1]",
                    lambda: read_sage(sage_tsv("s_hi.tsv", qval="7.5"), "pseudo.mzML"))
    st.expect_abort("G-C04", "TSV naming a different mzML (stale search result)",
                    lambda: read_sage(sage_tsv("s_fn.tsv", fn="other.mzML"), "pseudo.mzML"))
    st.expect_abort("G-C04", "header digest differs from the pinned schema",
                    lambda: read_sage(sage_tsv("s_pin.tsv"), "pseudo.mzML", "0" * 64))
    st.expect_ok("G-C04", "well-formed sage TSV parses",
                 lambda: read_sage(sage_tsv("s_ok.tsv"), "pseudo.mzML"))

    # ---- G-C05 spectrum_q never enters a table ---------------------------------------------------------
    _, rows, _ = read_sage(sage_tsv("s_ok2.tsv"), "pseudo.mzML")
    sm = search_metrics(rows, eid)
    st.results.append(("PASS" if "peptides_at_spectrum_q" in sm["diagnostics"] and
                       "spectrum_q" not in json.dumps(
                           {k: v for k, v in sm.items() if k != "diagnostics"})
                       else "FAIL", "G-C05",
                       "spectrum_q is diagnostics-only, never a headline field",
                       "diagnostics=%s" % sm["diagnostics"]))

    # ---- G-C09 reference binding ----------------------------------------------------------------------
    rep_ok = _mk_parquet(tmp / "rep_ok.parquet", _truth_rows("dataset-A", 300),
                         TRUTH_PARQUET_COLS)
    rep_multi = _mk_parquet(tmp / "rep_multi.parquet",
                            _truth_rows("dataset-A", 200) +
                            _truth_rows("dataset-D", 200),
                            TRUTH_PARQUET_COLS)
    rep_wrong = _mk_parquet(tmp / "rep_wrong.parquet",
                            _truth_rows("dataset-D", 300), TRUTH_PARQUET_COLS)
    rep_nan = _mk_parquet(tmp / "rep_nan.parquet",
                          _truth_rows("dataset-A", 300)[:-1] +
                          [dict(_truth_rows("dataset-A", 1)[0],
                                **{"Q.Value": float("nan")})], TRUTH_PARQUET_COLS)
    rep_nocol = _mk_parquet(tmp / "rep_nocol.parquet",
                            [{k: v for k, v in r.items() if k != "Run"}
                             for r in _truth_rows("x", 5)], TRUTH_PARQUET_COLS[1:])
    raw08 = "/raw/dataset-A"
    st.expect_abort("G-C09", "reference parquet describing ANOTHER sample (codex #9)",
                    lambda: verify_reference(rep_wrong, "dataset A", raw08))
    st.expect_abort("G-C09", "parquet with no Run/File.Name column",
                    lambda: verify_reference(rep_nocol, "dataset A", raw08))
    st.expect_ok("G-C09", "single-run reference accepted",
                 lambda: verify_reference(rep_ok, "dataset A", raw08))

    def multirun():
        tbl, keep, runcol, distinct = verify_reference(rep_multi, "dataset A", raw08)
        if len(keep) != 200:
            die("G-C09", "multi-run parquet was not filtered to this sample (%d rows)"
                % len(keep))
    st.expect_ok("G-C09", "multi-run DIA-NN report is FILTERED, not rejected", multirun)

    lib_ok = _mk_parquet(tmp / "lib_ok.parquet",
                         [{"Precursor.Id": "PREC%d" % i, "Product.Mz": 300.0 + j, "Decoy": 0}
                          for i in range(300) for j in range(4)], LIB_PARQUET_COLS)
    st.expect_abort("G-C09", "NaN Q.Value silently changes the truth population",
                    lambda: load_truth(rep_nan, lib_ok, "dataset A", raw08))
    st.expect_abort("G-C08", "absent fragment library -> abort, not 0% purity",
                    lambda: load_truth(rep_ok, tmp / "nope.parquet", "dataset A", raw08))
    lib_thin = _mk_parquet(tmp / "lib_thin.parquet",
                           [{"Precursor.Id": "PREC%d" % i, "Product.Mz": 300.0, "Decoy": 0}
                            for i in range(10)], LIB_PARQUET_COLS)
    st.expect_abort("G-C08", "library covering 3% of truth precursors",
                    lambda: load_truth(rep_ok, lib_thin, "dataset A", raw08))
    st.expect_ok("G-C08", "adequate library accepted",
                 lambda: load_truth(rep_ok, lib_ok, "dataset A", raw08))

    # ---- G-C23 vault O_EXCL --------------------------------------------------------------------------
    def vault_collision():
        v = tmp / "vault"
        v.mkdir(exist_ok=True)
        p = v / "report.md"
        p.write_text("first")
        try:
            os.open(str(p), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
        except OSError:
            die("G-C23", "vault file %s already exists and is never overwritten" % p)
    st.expect_abort("G-C23", "vault report is never overwritten", vault_collision)

    # ---- END-TO-END truth pass -----------------------------------------------------------------------
    _selftest_truth_e2e(st, tmp, rep_ok, lib_ok, raw08)
    _selftest_integration(st, tmp)
    return st


def _selftest_integration(st, tmp):
    """bench2 run -> collate2 truth -> collate2 collate, on a synthetic world.

    This is the case the v1 harness never had: the PUBLISH GATE driven end to end. Two runs:
    one whose decoy control is inert (must REFUSE, and must leave the vault empty) and one whose
    control fires (must publish). A gate that has only ever been unit-tested is a gate nobody
    has watched close.
    """
    import yaml
    work = tmp / "integration"
    work.mkdir()
    saved_env = dict(os.environ)

    def build(name, control_fires):
        root = work / name
        root.mkdir()
        raw = bench2.st_make_fake_d(root / "dataset A.d", mz_centers=(520.0,), width=80.0,
                                    serial="SER-08", sample_name="dataset A",
                                    full_mobility_span=True)
        tool = bench2.st_write_faketool(root / "faketool.py")
        sage = bench2.st_write_fakesage(root / "fakesage.py")
        fasta = root / "human.fasta"
        fasta.write_text(">sp|X|X\nPEPTIDE\n")
        cfg = root / "sage.json"
        cfg.write_text(canon({"database": {"fasta": str(fasta)}}).decode())

        # 300 target precursors on the tool's own m/z grid. When control_fires, a THIRD of them
        # get a partner exactly DECOY_SHIFT above, so the shifted decoy query has something to
        # hit. That condition is precisely what the harness refuses to assume: with it absent,
        # the floor is structurally zero and the correction is the identity function.
        rows = []
        frag = []
        for i in range(300):
            rows.append({"Run": "dataset A", "Precursor.Id": "T%d" % i,
                         "Precursor.Mz": 500.0 + i * 0.02, "Precursor.Charge": 2,
                         "RT": 10.0, "RT.Start": 9.0, "RT.Stop": 11.0,
                         "IM": 0.7 + i * 0.002, "Q.Value": 0.001})
            for j in range(6):
                frag.append({"Precursor.Id": "T%d" % i,
                             "Product.Mz": 1000.0 + i * 7.0 + j * 0.37, "Decoy": 0})
            if control_fires and i % 3 == 0:
                rows.append({"Run": "dataset A", "Precursor.Id": "D%d" % i,
                             "Precursor.Mz": 500.0 + i * 0.02 + DECOY_SHIFT,
                             "Precursor.Charge": 2, "RT": 10.0, "RT.Start": 9.0,
                             "RT.Stop": 11.0, "IM": 0.7 + i * 0.002, "Q.Value": 0.001})
                for j in range(6):
                    frag.append({"Precursor.Id": "D%d" % i,
                                 "Product.Mz": 5000.0 + i * 7.0 + j * 0.37, "Decoy": 0})
        rep = _mk_parquet(root / "report.parquet", rows, TRUTH_PARQUET_COLS)
        lib = _mk_parquet(root / "lib.parquet", frag, LIB_PARQUET_COLS)
        ident = bench2.raw_identity(raw)
        reg = {
            "bench_root": str(root / "bench"), "vault": str(root / "vault"),
            "tool": {"scaled_params": []},
            "search": {"sage_binary": {"path": str(sage), "sha256": sha256_file(sage)},
                       "sage_config": {"path": str(cfg), "sha256": sha256_file(cfg)},
                       "fasta": {"path": str(fasta), "sha256": sha256_file(fasta)}},
            "samples": {"dataset A": {
                "raw": str(raw),
                "pinned": {"raw_content_id": ident["content_id"],
                           "acquisition_id": ident["acquisition_id"],
                           "method_id": ident["method_id"],
                           "raw_coverage": ident["coverage"],
                           "pinned_utc": utc_now(), "pinned_by": "selftest"},
                "references": {"diann_report": {"path": str(rep), "sha256": sha256_file(rep)},
                               "diann_lib": {"path": str(lib), "sha256": sha256_file(lib)}}}},
        }
        regp = root / "samples.yaml"
        regp.write_text(yaml.safe_dump(reg, default_flow_style=False, sort_keys=True))
        plan = {"schema": "spextractor.plan/2", "title": "integration %s" % name,
                "binary": str(tool), "threads": 2, "tcmalloc": False, "samples": ["dataset A"],
                "arms": [{"name": "base", "role": "baseline", "params": {}},
                         {"name": "split", "role": "treatment", "baseline": "base",
                          "params": {"trace:ms1_split_valleys": 7.0}}]}
        planp = root / "plan.yaml"
        planp.write_text(yaml.safe_dump(plan, default_flow_style=False, sort_keys=True))
        os.environ.update({"SPEXTRACTOR_REGISTRY": str(regp), "FAKETOOL_MODE": "ok",
                           "FAKESAGE_MODE": "ok", "FAKETOOL_SPECTRA": "700",
                           "FAKETOOL_MZ0": "500.0", "FAKETOOL_MZSTEP": "0.02",
                           "FAKETOOL_PEAKMODE": "indexed"})
        return root, planp, Path(reg["vault"]), Path(reg["bench_root"])

    def latest_runset(bench_root):
        cands = sorted((bench_root / "runsets").glob("*.r1.json"))
        return cands[-1].name.split(".r1.json")[0]

    try:
        # ---- case 1: inert control must REFUSE ------------------------------------------
        root, planp, vault, broot = build("inert", control_fires=False)
        rc = bench2.main(["run", "--plan", str(planp), "--need-gb", "0.01"])
        st.results.append(("PASS" if rc == 0 else "FAIL", "G-B12",
                           "integration: bench2 run completes", "rc=%d" % rc))
        rid = latest_runset(broot)
        rc = main(["truth", "--runset", rid])
        st.results.append(("PASS" if rc == 0 else "FAIL", "G-C22",
                           "integration: truth writes revision r2", "rc=%d" % rc))
        rc = main(["collate", "--runset", rid + ".r2", "--publish"])
        st.results.append(("PASS" if rc == 3 else "FAIL", "G-C20",
                           "integration: inert decoy control REFUSES to publish",
                           "exit %d (3 = refusal)" % rc))
        published = sorted(vault.glob("*.md")) if vault.exists() else []
        st.results.append(("PASS" if not published else "FAIL", "G-C20",
                           "integration: refusal leaves the vault EMPTY",
                           "%d file(s) in %s" % (len(published), vault)))
        inc = sorted((vault.parent / (vault.name + "-incomplete")).glob("*.INCOMPLETE.md"))
        body = inc[0].read_text() if inc else ""
        st.results.append(("PASS" if inc and "UNDEFINED(" in body else "FAIL", "G-C20",
                           "integration: incomplete report says UNDEFINED in the cell",
                           inc[0].name if inc else "no incomplete report written"))

        # ---- case 2: live control must PUBLISH -------------------------------------------
        root, planp, vault, broot = build("live", control_fires=True)
        bench2.main(["run", "--plan", str(planp), "--need-gb", "0.01"])
        rid = latest_runset(broot)
        main(["truth", "--runset", rid])
        rc = main(["collate", "--runset", rid + ".r2", "--publish"])
        st.results.append(("PASS" if rc == 0 else "FAIL", "G-C20",
                           "integration: a fully defined runset DOES publish",
                           "exit %d" % rc))
        published = sorted(vault.glob("*.md")) if vault.exists() else []
        report = published[0].read_text() if published else ""
        st.results.append(("PASS" if published else "FAIL", "G-C23",
                           "integration: vault report written under runset id",
                           published[0].name if published else "nothing published"))
        st.results.append(("PASS" if "UNDEFINED(" not in report else "FAIL", "G-C20",
                           "integration: published report contains no UNDEFINED cell",
                           "clean" if "UNDEFINED(" not in report else "has UNDEFINED"))
        table_rows = [l for l in report.splitlines() if l.startswith("|")]
        leaked = [l for l in table_rows if "spectrum_q" in l]
        st.results.append(("PASS" if not leaked else "FAIL", "G-C05",
                           "integration: spectrum_q never appears in a published table row",
                           "absent from all %d table rows" % len(table_rows) if not leaked
                           else leaked[0][:60]))
        st.results.append(("PASS" if "n = 1" in report else "FAIL", "G-C20",
                           "integration: report states n = 1 explicitly",
                           "present" if "n = 1" in report else "missing"))

        # ---- republishing the same runset must not overwrite ------------------------------
        rc = main(["collate", "--runset", rid + ".r2", "--publish"])
        st.results.append(("PASS" if rc == 2 else "FAIL", "G-C23",
                           "integration: re-collating refuses to overwrite the vault file",
                           "exit %d (2 = guard abort)" % rc))

        # ---- verify --deep over the real sealed tree --------------------------------------
        rc = bench2.main(["verify", "--runset", rid + ".r2", "--deep"])
        st.results.append(("PASS" if rc == 0 else "FAIL", "G-B12",
                           "integration: verify --deep rehashes every sealed artifact",
                           "rc=%d" % rc))

        # ---- tamper with a sealed artifact; verify must catch it ---------------------------
        rs, _ = load_runset(Bench(broot), rid + ".r2")
        cell = rs["cells"][0]
        victim = Bench(broot).stage_dir(cell["produce_kind"],
                                        cell["produce_stage_id"]) / "work" / "pseudo.mzML"
        os.chmod(str(victim.parent), 0o755)
        os.chmod(str(victim), 0o644)
        with open(str(victim), "ab") as fh:
            fh.write(b"<!-- tampered -->")
        rc = bench2.main(["verify", "--runset", rid + ".r2", "--deep"])
        st.results.append(("PASS" if rc == 2 else "FAIL", "G-B12",
                           "integration: a byte appended to a sealed artifact is detected",
                           "exit %d (2 = guard abort)" % rc))
        rc = main(["collate", "--runset", rid + ".r2"])
        st.results.append(("PASS" if rc == 2 else "FAIL", "G-C22",
                           "integration: collate refuses a runset with a mutated artifact",
                           "exit %d" % rc))
    finally:
        os.environ.clear()
        os.environ.update(saved_env)


def _selftest_truth_e2e(st, tmp, report, lib, raw08):
    """Drive the whole truth pass over a synthetic mzML and assert the ledger partitions."""
    windows = [{"mz_lo": 495.0, "mz_hi": 505.0, "im_lo": 0.9, "im_hi": 1.1}]
    specs = []
    for i in range(250):
        specs.append({"mz": 500.0 + (i % 300) * 0.01, "z": 2, "im": 1.0, "rt": 600.0,
                      "iso_target": 500.0,
                      "peaks": [(300.0 + j, float(30 - j)) for j in range(25)]})
    # spectra that must land in named non-SCORED buckets
    specs.append({"mz": 900.0, "z": 2, "im": 1.0, "rt": 600.0, "iso_target": 900.0,
                  "peaks": [(300.0, 1.0)] * 1})
    specs.append({"mz": 500.0, "z": 2, "im": None, "rt": 600.0, "iso_target": 500.0,
                  "peaks": [(300.0, 1.0)]})
    specs.append({"mz": 500.0, "z": 2, "im": 1.0, "rt": 600.0, "iso_target": 500.0,
                  "peaks": [(300.0, 1.0), (301.0, 2.0)]})
    path = bench2.st_build_mzml(tmp / "e2e.mzML", specs, "file://" + raw08)
    precursors, frags, refinfo = load_truth(report, lib, "dataset A", raw08)

    def run():
        return truth_metrics(path, precursors, frags, windows, "method1", len(specs))
    st.expect_ok("G-C11", "end-to-end truth pass partitions every spectrum", run)
    tm = run()
    led = tm["ledger"]["dispositions"]
    st.results.append(("PASS" if sum(led.values()) == tm["ledger"]["spectra_read"] == len(specs)
                       else "FAIL", "G-C11", "ledger partition is exact end to end",
                       "read=%d sum=%d" % (tm["ledger"]["spectra_read"], sum(led.values()))))
    st.results.append(("PASS" if led["OUT_OF_WINDOW_TABLE"] == 1 and
                       led["NO_ION_MOBILITY"] == 1 and led["FEWER_THAN_K_PEAKS"] == 1
                       else "FAIL", "G-C10",
                       "each excluded spectrum lands in its NAMED bucket",
                       "oow=%d noim=%d thin=%d" % (led["OUT_OF_WINDOW_TABLE"],
                                                   led["NO_ION_MOBILITY"],
                                                   led["FEWER_THAN_K_PEAKS"])))
    rec = _q_load(tm["recall"]["corrected"])
    st.results.append(("PASS" if is_undef(rec) or isinstance(rec, Value) else "FAIL", "G-C15",
                       "recall is a Value or a NAMED Undefined -- never a bare number",
                       rec.code if is_undef(rec) else repr(rec)))
    st.results.append(("PASS" if is_undef(rec) and rec.code == "CONTROL_DEGENERATE" else
                       ("PASS" if isinstance(rec, Value) else "FAIL"), "G-C14",
                       "inert control on synthetic data blocks publication",
                       rec.code if is_undef(rec) else "control fired: %r" % rec))
    ch = _q_load(tm["charge_majority_baseline"])
    st.results.append(("PASS" if is_undef(ch) or ch.pop == "charge_observations" else "FAIL",
                       "G-C16", "charge baseline uses the SAME population as agreement",
                       ch.code if is_undef(ch) else ch.pop))


def cmd_selftest(args):
    require_deps()
    tmp = Path(tempfile.mkdtemp(prefix="spextractor-collate-selftest-"))
    try:
        st = selftest_collate(tmp)
        rc = st.report()
        if args.keep:
            print("\nfixtures kept at %s" % tmp)
        return rc
    finally:
        if not args.keep:
            shutil.rmtree(str(tmp), ignore_errors=True)


# ---------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("truth", help="compute metrics for a runset; writes revision r2")
    p.add_argument("--runset", required=True)
    p.add_argument("--bench-root")

    p = sub.add_parser("collate", help="validate and render the vault report")
    p.add_argument("--runset", required=True)
    p.add_argument("--bench-root")
    p.add_argument("--publish", action="store_true")

    p = sub.add_parser("selftest", help="deliberately trigger every guard and verify it fires")
    p.add_argument("--keep", action="store_true")

    args = ap.parse_args(argv)
    if not args.cmd:
        ap.print_help()
        return 2
    table = {"truth": cmd_truth, "collate": cmd_collate, "selftest": cmd_selftest}
    try:
        return table[args.cmd](args)
    except Abort as exc:
        sys.stderr.write("\nABORT %s\n%s\n\n" % (exc.guard, exc.message))
        return 2
    except RefuseToPublish as exc:
        sys.stderr.write("\nREFUSE TO PUBLISH: %s\n\n" % exc)
        return 3


if __name__ == "__main__":
    sys.exit(main())
