# Where the integer detector's memory goes

**Method.** A byte ledger over the window loop's data structures (`SPEXTRACTOR_MEM_LEDGER=1`,
`MemGuard` RAII around each allocation), recording PEAK CONCURRENT bytes per structure rather than a
sum over time -- windows and their band tasks run at once, so simultaneity is the whole question.
Measured on dataset D, 100 threads, 128-core node, `trace:detector=integer`. Process peak RSS 105 GB.

The first ledger run over-reported traces at 63.5 GB because `buildFragGrid` consumes every trace
profile and nothing subtracted it, making that one counter cumulative. The tell was that it exceeded
the total point-bytes in the file. Corrected below.

## 1. Per-structure inventory (measured)

| structure | peak concurrent | bytes/element | what it is |
|---|---|---|---|
| **traces (`xi`/`xv`/`xt` + `Trace`)** | **35.6 GB** | 12 B/point + 120 B/trace | detected traces, from detection until the grid consumes them |
| **prep (`frame_of`/`member`/`order`)** | **14.3 GB** | 4 + 1 B/peak, 4 B/seed | per-window preparation shared by the band tasks |
| **slab (`tof`/`inten`/`imq`)** | **12.0 GB** | 10 B/peak | the peaks themselves |
| EPD `list<Peak2D>` | 6.0 GB | ~48 B/point | valley-splitting round trip |
| emitted spectra | 5.0 GB | 16 B/peak | accumulate across ALL windows until the write |
| `visited` flags | 4.9 GB | 1 B/peak **per band task** | consume-and-mark state |
| EPD `MassTrace` + split | 3.6 GB | ~500 B/trace, doubled | `mts` and `split` coexist |
| FragGrid | 1.7 GB | 8 B/point | correlation grid |
| pdense scratch | 0.02 MB | 4 B/grid slot | per-thread |
| **sum of per-structure peaks** | 83.1 GB | | (not simultaneous) |
| **simultaneous peak of the total** | **68.6 GB** | | ledger's own peak |
| process peak RSS | 105 GB | | the 36 GB gap is allocator overhead and untracked OpenMS internals |

For comparison the OpenMS path's ledger accounts only 27 GB of its 162 GB RSS, because
`MassTraceDetection`'s band sub-maps and its own `MassTrace` objects are inside OpenMS and not
instrumented -- its materialised `PeakMap` alone is 24.0 GB against the slab's 12.0 GB for the same
peaks.

## 2. Duplication, with lifetimes

| pair | simultaneous? | cost | verdict |
|---|---|---|---|
| `TracePrep.member[P]` vs `slab.inten[P]` | yes, whole window | 1.2 GB | **pure duplication**: `member[k]` IS `inten[k] > noise`, an O(1) recompute |
| `TracePrep.frame_of[P]` vs `slab.frame_off[F+1]` | yes, whole window | 4.8 GB | **derivable**: a peak's frame is a binary search over 1,400 boundaries, and every hot-path use already knows the frame |
| `Trace.xt` vs `slab.tof` | yes -- the slab outlives the traces | ~10 GB | **duplication**: `xt` copies bins the slab still holds; storing the PEAK INDEX would give bin, frame and intensity |
| `Trace.xv` vs `slab.inten` | yes | ~10 GB | same |
| `Trace.{xi,xv,xt}` vs EPD `list<Peak2D>` | brief, per chunk | 12 vs 48 B/point | the profiles are freed as each list is built, so only one chunk's worth overlaps |
| EPD `mts` vs `split` | yes, per chunk | 3.6 GB | `splitValleys` returns a new vector while the input is alive |
| slab vs compact store | no | -- | `toSlab` consumes the frames as it goes |

## 3. Repeated passes over the same data

Per peak, from slab to spectrum: (1) `toSlab` conversion; (2) `prepareTracing` frame/member pass;
(3) `prepareTracing` order-fill pass; (4) the intensity sort; (5) N candidate visits during
extension. **Passes 2 and 3 are two sweeps that could be one.** Per trace point: written once by the
detector, read once into a `list<Peak2D>`, read once back out of the split `MassTrace`, read once by
`buildFragGrid` -- four touches of data that is 12 bytes wide, three of them only because valley
splitting round-trips through OpenMS objects.

## 4. `visited`: the clearest single defect

Each band task allocates `vector<char> visited(P)` at the WHOLE window's peak count, although it
seeds only within its own twelfth of the flight-time range. Measured 4.9 GB concurrent; the
reviewers' static estimate is 14.4 GB at full concurrency. As a bitset it is 8x smaller. **And a
per-band `visited` is not merely wasteful, it is wrong**: two bands can each claim the same
out-of-core peak, because extension is unrestricted while the flag is not shared. One shared
per-window bitset is both 12x smaller and closer to the single-threaded semantics it reimplements --
at the cost of atomics on the mark.

## 5. Proposed changes, ranked by measured bytes per unit of risk

| # | change | saving | output | risk |
|---|---|---|---|---|
| 1 | `member` -> inline `inten[k] > noise` | 1.2 GB | none | none |
| 2 | drop `frame_of`: carry the frame in `got` (the loop already knows it) and binary-search the two places that do not | 4.8 GB | none | none |
| 3 | `visited` -> shared per-window bitset with atomic mark | 4.3 GB, and removes the cross-band double-claim | changes output (fixes a real defect) | low, needs both engines |
| 4 | `Trace` stores the PEAK INDEX instead of `xi`+`xt` | ~10 GB | none | low: the slab must outlive the traces, which it already does |
| 5 | one arena allocation per trace instead of three vectors | ~11 GB of allocator headers, and 3e7 small allocations become 1e7 | none | medium |
| 6 | valley splitting without the `list<Peak2D>` round trip | up to 9.6 GB | none if the split is identical | high: it is OpenMS's own step |

Items 1, 2 and 4 are output-neutral by construction and together are ~16 GB of the 68.6 GB ledger
peak. Item 3 is the one worth doing for correctness rather than for bytes.

**What this does not explain.** 105 GB RSS against a 68.6 GB ledger peak. Some is allocator
fragmentation from ~3e7 small vector allocations per window (item 5 attacks it directly); some is
OpenMS internals inside `ElutionPeakDetection`. Before spending effort on item 5 the right
measurement is `MALLOC_ARENA_MAX` or a different allocator, which is a one-run experiment and was
already falsified once for speed on this codebase -- but never tested for peak RSS.
