# The third axis: indexing ion mobility, and what it is worth

**Question (user, 2026-09-03):** the RT axis is a frame index and the m/z axis is the flight-time
bin -- can ion mobility be indexed the same way, and what does it buy?

**Short answer: yes, and the mobility index matters less for its own sake than for what it unlocks.
Indexed on its own it saves ~6 bytes per trace and no measurable time. Combined with the RT index it
turns the scoring gate from a scan into a lookup, and the gate is currently streaming an estimated
83 TB of memory per dataset D run to discard 99.4% of what it reads.**

## 1. The mobility axis is discrete, and for the same reason

A TIMS device ramps a field and reads out in SCANS: a frame has ~1,000 of them, and 1/K0 is a
calibrated function of (scan, frame). So mobility is an integer at the source, and the tool has been
storing a quantised version of the derived quantity -- `uint16` over [0.4, 1.8], the same shape of
approximation the flight-time work removed for m/z.

**The axis can be derived without the vendor calibration.** Within one frame the distinct 1/K0
values ARE the scans, so ranking them recovers the scan index exactly: sort the distinct values of a
frame (about a thousand), and a peak's rank is its scan. This mirrors how the RT axis is built and
keeps the tool independent of the TIMS calibration table, which is a different table from the m/z
one and is not needed to say that two peaks came from the same scan.

**A per-frame index, not a global one.** The TIMS calibration is frame-dependent, so the same scan
maps to slightly different 1/K0 in different frames -- which is exactly why comparing SCANS is
better than comparing calibrated mobilities: the gate stops spending part of its tolerance on
calibration drift. As with m/z, the calibrated value is needed only when writing.

## 2. What indexing mobility alone is worth

| | now | indexed | saving |
|---|---|---|---|
| per peak | `uint16` quantised 1/K0 | `uint16` scan | 0 |
| per trace | `double im` (8 B) | `uint16` scan (2 B) | 6 B |
| ~15 M traces in flight | | | **~90 MB** |
| gate comparison | `fabs(f.im - pc.im) > delta_im` on doubles | integer difference | negligible |

**On its own this is not worth doing.** 90 MB against a 100-150 GB peak, and an integer compare
replacing a double compare in a loop that is memory-bound rather than ALU-bound.

## 3. What it is worth combined with the RT index: the 158x gate

Measured on dataset D (`score gate` instrumentation, 2026-09-03):

    710,046 precursors visited 867,449,463,431 fragments  (1,221,680 per precursor)
    5,480,062,634 passed the RT gate                       (158.3x scanned per survivor)

Fragments are sorted by mobility only, so the mobility window is a binary search but the RT gate is
a LINEAR SCAN of everything inside it. 99.4% of the visits are discarded.

The cost is memory bandwidth, not arithmetic. `Trace` is 96 bytes, and the scan strides one per
visit:

| | traffic per dataset D run |
|---|---|
| scanning `Trace` records at 96 B | **~83 TB** |
| scanning an 8 B parallel array instead | ~6.9 TB |
| touching only survivors, via a 2-D index | **~0.04 TB** |

That reframes the score stage. It is 35% of window-loop time and its dominant term is not the dot
products, it is walking 867 billion records to find 5.5 billion.

## 4. The two changes, in order of value per unit of risk

**(a) Parallel arrays for the gate fields. One line, no algorithm change.** `frag_im` already exists
as a parallel `vector<double>` purely so the mobility bound can be bisected. Adding the same for the
apex frame index means the RT scan streams 8 B per element instead of 96. **Expected: ~12x less
traffic in the scan; if the scan is bandwidth-bound as the numbers suggest, most of the gate's share
of the score stage.** Output-neutral by construction: the same fragments are visited in the same
order. Falsifier: the score stage does not move, which would mean the scan was not bandwidth-bound
and the dot products dominate after all.

**(b) A 2-D bucket index over (scan, frame).** With both axes integer, the fragment list can be
bucketed by mobility scan, and within a bucket sorted by apex frame, so the gate becomes: three
bucket lookups, then a bisect on the frame range. **Expected: the 158x collapses to ~1x, removing
~99.4% of the visits.** Output-neutral: the fragments that pass are identical; only the ones that
are never examined change. Falsifier: the same.

Neither depends on the mobility axis being the vendor's scan number rather than the current
quantised 1/K0 -- bucketing works on either. **The scan index is the right representation, but the
gate fix is what pays**, and it should be measured first so the two are not confounded.

## 5. Estimated effect, stated as a range

The score stage is 3,853 s of 10,960 s of window-elapsed time on dataset D (35%), and the window loop is
~80% of wall. If the gate scan is 60-85% of the score stage -- which the bandwidth arithmetic
supports but has not been isolated -- then:

| change | window-loop time | expected dataset D wall |
|---|---|---|
| now | -- | 12:53 |
| (a) parallel arrays | −15 to −25% of the score stage | 12:00 - 12:30 |
| (a) + (b) index | −60 to −85% of the score stage | **9:30 - 11:00** |
| mobility indexed as scans | ~0 | unchanged |

Memory: (a) adds one 8 B array per window (~70 MB), (b) adds bucket offsets of similar size, and
the mobility index removes ~90 MB from the traces. Net effect on the 100-150 GB peak: none worth
predicting.

**How to falsify the whole plan cheaply:** the instrumentation that produced the 158x is already in
the binary. Run (a) alone and read the score-stage timer. If it does not move, the scan was never
the cost, and (b) will not help either -- in which case the remaining time is genuinely in the
correlation arithmetic, and the only lever left is reducing how many (precursor, fragment) pairs are
correlated at all, which is a scientific decision about fragment sharing rather than an engineering
one.
