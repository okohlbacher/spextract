# The trace struct has one degree of freedom per point, and stores three

**Origin.** User, 2026-09-04: a trace is contiguous by definition, so it needs only its entry point
on the RT axis and its intensities -- not a coordinate per point. Same for the flight-time bin.

## 1. What is stored today, and why it is wrong

```cpp
struct Trace {                       // src/spextractor.cpp:243
  double mz, rt, im, intensity;      //  32 B  per trace, all DERIVABLE from what follows
  uint32_t tof; double b;            //  16 B  per trace (b is the apex FRAME's factor: derivable)
  vector<uint32_t> xi;               //   4 B  per point: frame index      <- redundant
  vector<float>    xv;               //   4 B  per point: intensity        <- the only payload
  vector<uint32_t> xt;               //   4 B  per point: flight-time bin  <- redundant
};                                   // 120 B struct + 3 malloc headers + 12 B/point
```

Measured: 35.6 GB peak, the largest structure in the integer path. A trace point has ONE degree of
freedom -- its intensity -- because the frame is implied by position and the bin is the trace's.
`xi` and `xt` encode information that is either constant along the trace (the bin, within
tolerance) or arithmetic on the position (the frame). The doubles are caches of things a few
instructions recompute.

## 2. "Contiguous" needs one refinement, and it decides the design

Two facts from the data (dataset D):

* **33,553 frames, 17,343 distinct RTs, 1,343 MS1 frames, 24 windows.** A window's CONSECUTIVE
  frames are ~13 apart on the GLOBAL RT axis, because the MS1 frame and the other windows' frames
  interleave in every cycle. So a trace is NOT contiguous on the global axis. It IS contiguous in
  its own window's frame sequence (`PeakSlab::rt_index`, window-local index 0..F-1), and for MS1
  traces in the MS1 map's sequence. **The entry point must be window-local**, and `rt_index[]`
  maps it to the global axis when a global RT is needed.
* **Traces have gaps.** The detector, following OpenMS, tolerates up to 5 consecutive missing
  frames and accepts a trace if points/frames >= 0.5 over its span. So "contiguous" means "a span
  with holes", and the holes must be representable. **Dense over the span, zero = missing.** By
  the sample-rate rule at most half the entries are padding, so the worst case is 8 B per real
  point of intensity storage against 12 B per point today, and the typical case (few gaps) is ~4.
  Every consumer already skips zero intensities (the scorer's `overlap` counts only non-zero
  precursor slots; `mean`/`invnorm` sums are unaffected by zeros), so the semantics are unchanged.

## 3. The redesigned trace

```cpp
struct Trace
{
  uint32_t frame0;     // first frame, WINDOW-LOCAL (MS1: index in the MS1 map)
  uint32_t off;        // offset of this trace's intensities in the window's arena
  uint32_t tof;        // the trace's flight-time bin (apex peak's, or centroid, per estimator)
  uint16_t len;        // frames spanned (traces are < 65k frames; the gradient is 1,400)
  uint16_t apex;       // offset of the apex within the span
  uint16_t imq;        // quantised mobility centroid
};                     // 18 B per trace, no allocation per trace
// per window:  vector<float> arena;   // all traces' intensities, 4 B per spanned frame
```

Everything the old struct cached is an accessor over the slab and the arena:

| old field | now | cost |
|---|---|---|
| `rt` | `rtAxis()[sl.rt_index[frame0 + apex]]` | two loads |
| `mz` | `tofAxis().mzOf(tof, sl.b[frame0 + apex])` | one calibration call, only at gate/export time |
| `im` | `IM_LO + imq / IM_Q` | one multiply |
| `intensity` | `arena[off + apex]` | one load |
| `b` | `sl.b[frame0 + apex]` | one load |
| `xi[i]` | `frame0 + i` | arithmetic |
| `xt[i]` | not needed: the bin is per trace; the apex's measured bin IS `tof` | -- |
| `xv[i]` | `arena[off + i]` | one load |

**Bytes.** ~7.5e6 traces per window, ~10.7 points per trace measured, call it 14 spanned frames
with gaps: today 120 + 48 + 12x10.7 = ~296 B/trace; redesigned 18 + 4x14 = 74 B/trace. **4x
smaller, and 3 mallocs per trace become 0** (the arena grows by append inside one band task, then
is concatenated per window). The measured 35.6 GB becomes ~9 GB.

## 4. What it does to the consumers -- the part that makes this more than a compaction

* **`FragGrid` disappears.** Its job was to resample every fragment onto a common per-window RT
  grid as CSR `(slot, intensity)` pairs. With dense window-local traces the grid IS the window's
  frame sequence and the CSR IS the arena: `flat[k]` = `arena[off + k]`, slot = `frame0 + k`.
  No union, no `slot_of_frame`, no copy: 1.7 GB and one pass over every point gone. `mean` and
  `invnorm` stay as two per-trace doubles computed once. The precursor-to-fragment nearest-frame
  table (`nearest_slot`) stays, now indexed by window-local frame -- ~1,400 entries.
* **The scorer's inner loop becomes a contiguous walk** `for k in [0,len)` reading
  `pdense[frame0 + k]` and `arena[off + k]`: two sequential streams, no gather through a slot
  index. Skipping `arena == 0` preserves `overlap` exactly.
* **`xicCorr` (isotope partners)** is an overlap of two `[frame0, frame0+len)` intervals and an
  aligned dot product -- no merge of index lists.
* **`atrousLevels` / FWHM** read `rtAxis()[rt_index[frame0 + i]]` where they read `rtAt(i)`.
* **Valley splitting** (`splitIntegerTraces`) builds its `MassTrace` from `(rt_index, arena)` and
  the split traces come back as spans of the parent: a child is `(frame0 + a, len b-a, same tof)`
  -- it can be represented WITHOUT copying intensities, as a sub-span of the parent's arena range.
  This removes the post-split re-conversion copy.
* **The OpenMS detector path** produces `MassTrace` objects whose peaks are the map's spectra;
  `toTrace()` finds `frame0` as the index of the first peak's RT in the map's spectrum sequence
  (the map is RT-sorted, so a binary search over ~1,400 RTs) and writes intensities into the
  arena with zeros for gaps. Both detectors then produce the SAME structure.
* **Precursor inference and assembly** read `mz`/`rt`/`im`/`intensity` through accessors that take
  the slab; the hot fragment gate already uses SoA arrays (`frag_rt`, `frag_im`), which become
  per-window arrays filled once from the accessors.

## 5. What could go wrong, and how to know

| risk | symptom | check |
|---|---|---|
| `len` overflows uint16 | a trace spanning > 65,535 frames: impossible on a 1,400-frame gradient, but assert it | `assert(len < 65536)` at construction, abort not truncate |
| a consumer relied on per-point bins (`xt`) | export m/z of split traces reverts to a centroid | the apex bin is `tof`; the ONLY consumer of per-point bins was `splitIntegerTraces` giving EPD per-peak m/z, and EPD ignores per-peak m/z (it splits on RT/intensity; SNR filter off) -- so `toTrace()` after a split must take the apex bin from the PARENT's slab lookup, not from the MassTrace |
| zero-padding changes a statistic | `mean`/`invnorm`/`overlap`/sample-rate differ | zeros contribute nothing to sums; `overlap` tests non-zero; the sample-rate test is done in the detector before the span is built -- verify by digest against the current integer output on dataset D (the integer path is not gated on bit-identity, but THIS change should be, because it is meant to be pure representation) |
| MS1 traces need the MS1 map's frame sequence, not a window's | wrong `rt` for precursors | one `PeakSlab`-like frame table for the MS1 map (`rt_index` only) |
| arena growth reallocates while tasks hold offsets | offsets stay valid (they are indices, not pointers) but concurrent appends race | each BAND task appends to its OWN arena; the window concatenates them after the band tasks join, rebasing `off` -- one add per trace |

**Falsifier for the whole redesign:** the spectrum digest of the integer path on dataset D must be
IDENTICAL before and after (this is pure representation; any difference is a bug), AND the
`traces` ledger line must fall from 35.6 GB to under 12 GB. A digest match with no memory drop
means the arithmetic in section 3 is wrong; a memory drop with a digest change means a consumer
was not as zero-tolerant as section 4 claims.

---

# Revision 2 (2026-09-04) — after adversarial review

Reviewed by codex and vibe (kimi was silent on the first pass and is being retried); full texts in
`docs/reviews/trace-struct-2026-09-04/`. Both reviewers reached the same verdict independently:
**the representation is viable and most of the memory win is real, but the design as written is
NOT pure representation.** Sections 2-5 above are corrected as follows.

## R1. The zero sentinel is a PRESENCE mask, and every consumer must honour it

The claim "every consumer already skips zeros" was false. A zero entry is semantically visible in:
`overlap` (counts every stored fragment point whose precursor slot is non-zero -- a padding zero
at such a slot would be counted), `touched` (pushes a slot whenever the precursor value was zero
before adding, so a zero point still pushes), `atrousLevels` (the median inter-point gap collapses
to the frame interval), both FWHM routines (a padding zero becomes an artificial half-max
crossing), the wavelet smoother (runs over the packed observed sequence and needs uniform sampling
of REAL points), `xicCorr` (intersects real points), and ElutionPeakDetection (must be handed real
points only: it derives scan time and smooths over the sequence it is given).

The zero sentinel is valid as presence because every observed intensity is strictly greater than
`noise_threshold_int >= 0`. **The rule is: no consumer walks the span; every consumer walks the REAL
points in chronological order, skipping zeros.** With that rule the floating-point sums are visited
in the same order as today, which is what makes bit-identity achievable (codex §7).

`npts` (real point count, `uint16`) is stored next to `len`, because the `< 5` / `< 3` tests and
`G` need it and scanning for it would be an O(len) hidden cost.

## R2. `G` and the MS1-to-MS2 alignment survive as small tables, not as the frame sequence

"The grid IS the window's frame sequence" was the riskiest line (codex). Two things must stay:

* **`G` = the number of DISTINCT frames that carry at least one real fragment point**, not the
  window's frame count. It is the Pearson denominator and changing it changes every correlation.
  Kept as a bitmask over the window's frames (~1,400 bits) plus its popcount.
* **The precursor-to-fragment nearest-frame table** stays, computed per window from the two RT
  tables: `nearest_frag_local[ms1_local]`, ties to the LATER fragment frame, exactly as today's
  `nearest_slot`. An MS1-local index and a window-local index with the same value are different
  times; the mapping goes through `rt_index[]` on both sides.

`FragGrid.flat` and `FragGrid.off` are deleted (the arena and `(off, len)` replace them);
`mean`/`invnorm` per fragment stay, computed once over real points.

## R3. Traces carry a frame-table CONTEXT, and some scalars are not derivable

* A bare `Trace` cannot say whether `frame0` indexes the MS1 map or a window. Accessors take an
  explicit `const FrameTable&` (`rt_index` + `b`); there are two tables per window in flight, the
  MS1 map's and the window's. No type split -- one struct, an explicit argument.
* The OpenMS detector's `mz` under the `median`/`mean` estimators is a function of per-point m/z
  the redesign does not keep, and its `rt` is the CENTROID RT, not the apex RT. Neither is
  derivable. **Step 1 keeps `mz`, `rt`, `im` as stored doubles** (24 B/trace). The points are the
  bulk; the scalars can be revisited once the point encoding has landed and been gated.
* A split child's exported m/z today comes from the child's own measured apex peak. With a single
  per-trace bin the child would inherit the parent's, which codex lists as an output change. Step 1
  therefore keeps a per-frame bin array in a SECOND arena (`uint32`, 4 B per spanned frame) whose
  lifetime ends when valley splitting has assigned each child its apex bin. It exists only between
  detection and splitting.

## R4. Arena sequencing (codex §4, adopted verbatim as the implementation order)

1. each band task appends to a PRIVATE arena; its `off` values are band-local;
2. the detection taskloop joins (the existing taskloop has the implicit taskgroup barrier);
3. prefix-sum the band arena sizes; allocate the window arena once; move each band arena into
   its range and add its base to every `off`; concatenate records in band order;
4. chunked EPD tasks run against the now IMMUTABLE arena; each child stores absolute offsets;
5. EPD joins; concatenate children; canonical sort of the records (safe: offsets are absolute);
6. build the used-frame mask, `G`, `mean`/`invnorm`, and the SoA gate arrays;
7. scorer tasks start only now; the arena and BOTH frame tables (`rt_index`, `b`) live until the
   scoring taskgroup joins. **Today `wslab.b` is freed before scoring** (~line 3202); that is a bug
   under the new accessors and is fixed as part of the change.
8. `assert(arena.size() < 2^32)` before any `off` is written; `assert(len < 65536)`.

## R5. Corrected arithmetic and bound

`sizeof` of the record is 20 B with natural alignment, not 18. The padding factor is bounded by
`1/min_sample_rate`, not 2: the sample-rate test is against frames visited minus trailing misses,
which is at least the span, so `span <= floor(npts / r)`. At the default `r = 0.5` the worst case
is 2x; at `r = 0.3` it is 3.3x; at `r = 0` there is no bound and the span can be the whole window.
Realistic estimate at 10.7 real points and ~75% occupancy: 20 + 4 x 14.3 = 77 B/trace against
248 B/trace tracked today (the 35.6 GB ledger excludes malloc headers), i.e. **35.6 GB -> ~11 GB
typical, ~15 GB at the r = 0.5 worst case.** The "~9 GB" above was unsupported.

## R6. What is deleted, and what is deliberately kept

Deleted: `xi`, `xt` (the persistent one), the three per-trace vectors and their three mallocs,
`FragGrid.flat`/`off`/`slot_of_frame`. Kept, with the reviewers' reasons: `off` (a prefix sum only
works if records are never re-sorted, and they are); `apex` (2 B against a 10-20 float scan at
every sort/gate/emit); `imq` for MS2 (already the compact form of a quantity that is needed);
`len` (dropping it needs an offsets array and costs random access).

## R7. The falsifier, staged

The spectrum-list SHA-256 of the integer path on dataset D must be identical before and after -- this IS
a pure representation change once R1-R3 hold, and codex confirms the arena concatenation does not
by itself alter floating-point order. Because the existing `[det]` trace digest hashes only
`mz`/`rt`/`intensity` and is order-insensitive, three staged checks precede it so a mismatch can be
localised rather than debugged from the output: (1) exact trace count and the exact ordered
`(local frame, intensity bits)` tuples before valley splitting; (2) exact child boundaries and
scalar bits after it; (3) exact used-frame mask, `G`, `mean` and `invnorm`. A digest match with no
drop in the `traces` ledger line means R5's arithmetic is wrong; a memory drop with a digest change
means a consumer still walks the span.

## Implementation status (2026-09-04, 06:15)

Implemented in commit `3a1d91e` on `coverage-analysis-2026-07-24`, exactly as revision 2 specifies:
`Trace` = scalars + `(st, frame0, off, len, npts, apex)`; `TraceStore` per window and for the MS1
map (`rt_index`, `b`, `frame_rt`, `inten` arena, `bins` arena until valley splitting);
`makeSpan()`/`packReal()`; `FragStats` replaces `FragGrid`; band and chunk tasks append to private
stores absorbed after each join with offsets rebased; the MS1 release compacts the arena. Compiles
clean. **The gate is running:** the 8-check suite, then the dataset D spectrum digest of BOTH detectors
against the pre-refactor outputs, then the memory ledger's `traces` line.


## Gate result (2026-09-04, 06:41) — PASSED on both detectors

| dataset D | digest vs pre-span | ledger simultaneous peak | visited | traces | process RSS |
|---|---|---|---|---|---|
| OpenMS detector | **IDENTICAL** (656,254 spectra) | 68.6 -> 32.4 GB | -- | 35.6 -> 7.3 GB | 164 -> 147 GB |
| integer detector | **IDENTICAL** (656,256 spectra) | 68.6 -> 27.9 GB | 4.87 -> 0.61 GB | 35.6 -> 6.9 GB | 105 -> 88 GB |

So the span representation IS pure representation, as claimed, once one bug introduced by the
refactor was removed: `splitIntegerTraces` had gained a walk over the parent traces to look a
child's flight-time bin up by frame, but parents within a chunk are ordered by mobility and m/z,
not frame, and their frame ranges overlap freely. The pre-span code never looked anything up -- the
bin travels inside the MassTrace as each peak's m/z -- and restoring that fixed it. The OpenMS arm
being digest-identical WHILE the integer arm lost 3.3% of peptides is what localised it, since that
arm exercises every shared consumer.

Landed together with the ledger corrections it made measurable: `visited` as a bitset (4.87 -> 0.61
GB), the slab and prep released before valley splitting, the trace bytes guarded rather than
add-only (the old 35.6 GB line was cumulative over the run, not a peak), and the phantom 6.0 GB
`list<Peak2D>` line replaced by the two `MassTrace` payloads that really coexist.

Seed fraction, measured: **100% of above-noise peaks are seeds** (`min_apex = snr x noise` with
snr 1.0 equals the membership threshold), so `order` is P-sized, not P/4 -- the exact-reserve item
is worth its full 2.4 GB, and `member` is redundant with the seed list as well as with `inten`.
