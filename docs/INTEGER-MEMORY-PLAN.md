The plan below is written against the **current** file (`<repo>/src/spextract.cpp`, 3853 lines), not against the file the eight analyses read. That distinction is load-bearing and is the first thing in the document.

---

# SpeXtract integer detector: memory and lifetime plan

**Reading note that changes everything below.** The eight analyses and both skeptic passes were run against the pre-redesign file (3754 lines). Since then four commits landed — `1ed875b` (design), `a3ce034` (revision 2), `3a1d91e` (implementation), `3bdcb60` (comment) — all *after* `9adcfcf`, the commit that produced the measured ledger. **The 68.6 GB / 105 GB ledger of record does not describe the code in the tree**, and no benchmark has been run since the redesign. Two of the eight dimensions are now about structures that no longer exist:

* `Trace::xi/xv/xt` (12 B/point, 3 mallocs/trace, the 35.6 GB line) is gone. A trace is now a span `[off, off+len)` into a per-store `vector<float>` arena (`struct Trace` :232-268, `struct TraceStore` :274-307, `makeSpan` :322-342).
* `FragGrid.flat` / `FragGrid.off` (the 1.7 GB line, and the CSR copy of every point) is gone. `buildFragStats` (:1386-1422) reads the arena in place; only `mean`, `invnorm`, `nearest_local`, `used` remain.

Everything else in the ledger — slab, prep, visited, EPD list, emitted spectra — is byte-for-byte the same code. So the ledger is **half stale and half valid**, and the first experiment is not an optimisation, it is re-taking the instrument.

---

## 1. The point encoding

**A trace point does not cost 12 bytes today.** It costs **4 B** in `TraceStore::inten`, one float per *spanned* frame with zero meaning missing (`Trace::xv` :308, `real()` :264). During the integer detector and valley splitting only, it costs **8 B**, because `TraceStore::bins` (:280) carries a parallel `uint32_t` per spanned frame; that arena is freed at :3298, immediately after `splitIntegerTraces` returns.

What is left is not the points. It is the struct:

```
struct Trace  (:232-268)          bytes
  mz, rt, im, intensity            32     four doubles
  tof (u32) + 4 B padding           8
  b (double)                        8
  st (const TraceStore*)            8
  frame0, off (u32 each)            8
  len, npts, apex (u16 each) + pad  8
                                   ---
                                    72
```

At the measured ~10.7 real points and ~14 spanned frames per trace (`docs/TRACE-STRUCT-REDESIGN.md` §3), a post-split trace costs **72 B struct + 56 B arena = 128 B**, of which **56 % is the struct**. The design projected 18 + 56 = 74 B. The implemented version is **1.7× its own projection**, and the entire excess is the six cached scalars that revision 2 kept so the OpenMS path could share the type.

Two of those six are provably redundant on the integer path:

| field | equals | where set | verdict |
|---|---|---|---|
| `intensity` (8 B) | `(double)xv(apex)` | :1016 `sl.inten[kap]`; `toTrace` :395 `mt.getMaxIntensity(false)` | derivable, both paths |
| `b` (8 B) | `st->b[frame0 + apex]` | :1012 `sl.b[frame_of[kap]]`; :1104 `wst.b[cf]` | derivable on the integer path only — `store.b` is **empty** on the OpenMS path (:3309 `setFrames(ri, nullptr)`) and `b` stays 0.0 there, which `exportMz_` relies on |
| `mz`, `rt`, `im` | — | EPD returns *centroids* (`toTrace` :390-394), not apex values | **not** derivable after valley splitting |

Sketch (the `intensity` half only — it is unconditional):

```cpp
// struct Trace: delete `double intensity;`
double intens() const { return (double)xv(apex); }   // st must be set (:3292/:3297)
```
…and rewrite the ~8 read sites (`sort` comparator :3335, `frag_*` build :3338, `weighted_`, assembly). `b` needs its own accessor with the empty-store fallback and its own check that the OpenMS path's 0.0 semantics survive; do not bundle them.

**Saving: 8 B/trace unconditionally, 16 B/trace if `b` also survives review — 1.6 GB / 3.2 GB at ~2.0e8 concurrent post-split traces.** That trace count is itself derived from the pre-redesign ledger and must be re-measured (Experiment 0). This is a low-priority item precisely because the number rests on an unmeasured population; it is listed so nobody re-derives it from scratch.

**What should NOT be attempted:** removing `bins`. It is 4 B/spanned frame and duplicates `sl.tof`, but recovering a trace's bin per frame needs the *peak index*, which is also 4 B. Same width, more code. And the slab is alive during split today anyway (see §3), so `bins` is a pure lifetime artefact that the lifetime fix removes for free.

---

## 2. Compaction

Anchor for every row: the slab line is the only internally consistent one in the measured ledger. 12.0 GB at exactly 10 B/peak (`toSlab` :772 reserves exactly; tof 4 + inten 4 + imq 2) gives **Pc ≈ 1.20e9 concurrent peaks** against 1.26e9 in the file — i.e. essentially every window's peaks are resident at once, and the admission gate is not binding at that instant.

**The prep line is impossible and must not be used as a basis.** `prepareTracing` (:855-877) allocates `frame_of` 4 B/peak + `member` 1 B/peak + `frame_live` F bytes + `order` at 4 B/seed whose capacity can never exceed P (`reserve(P/4)` :867 doubles P/4 → P/2 → P; size ≤ P by construction). Ceiling **9 B/peak = 10.8 GB at Pc**. The reported 14.3 GB is above the ceiling. Cause: `MemGuard g_prep` is declared at :3195 in the window body while `tprep` is declared at :3263 inside the `if (integer_detector)` block that closes at :3304 — the guard books freed bytes through split, grid, scoring and emission. Re-take it before quoting it.

| structure | current | proposed | saving | output effect | risk |
|---|---|---|---|---|---|
| **`visited`** (:894) `vector<char>(P)` per band task | 1 B/peak × live band tasks — **4.90 GB** | `vector<uint64_t>((P+63)/64)` + 2 inline helpers | **4.29 GB** (both skeptics agree exactly) | none by construction — set-identical, digest-gated | low |
| **`frame_of`** (:850) | 4 B/peak — **4.8 GB nominal** | delete; frame comes from the loop variable at every hot site, one `upper_bound`−1 over `frame_off` per *seed* | **1.6–4.8 GB.** Skeptics split: 4.8 GB is the nominal ceiling; 1.6 GB is the floor if only windows actually inside :3263-3292 count. Take the ledger delta after the guard is re-scoped. | none by construction (four value-identical substitutions) | low |
| **`member`** (:851) | 1 B/peak — **1.2 GB nominal** | delete; recompute `(double)sl.inten[k] > noise` at :928 (pass `noise` in) | **0.4–1.2 GB**, same scope disagreement | none by construction — but line 869 must keep **both** conjuncts (`> noise && > min_apex`); `ms2_chrom_peak_snr` may be < 1 | low |
| **`order`** (:853) | 4 B/seed, capacity ∈ {P/4, P/2, P} | count seeds, reserve exactly | **0 – 2.4 GB**, expected ~1.2 GB. Entirely contingent on the unmeasured seed fraction *m*. | none by construction | low, but **measure `m` first** |
| **slab** (:750-760) | 10 B/peak — **12.0 GB** | drop sub-noise peaks in `toSlab`; carry raw `tof` min/max and an F-sized raw occupancy count on `PeakSlab` | **(1−m) × (12.0 + frame_of + member)**. At m = 0.6 → ~7 GB; at m = 0.9 → ~1.8 GB. | claimed none-by-construction, but **only if** band edges come from the RAW tof range (:3260-3261 → :3277-3278) and `SPEXTRACT_INT_EMPTY_RAW` (:971) reads a stored raw count. Get either wrong and every band's seed partition moves silently. | **medium** |

**Not worth its complexity:** band-scoped `visited` (per-band-per-frame offset tables, ~20 lines, needs a reach histogram to pick the margin) buys 0.55 GB *on top of* the bitset's 4.29 GB. A shared-across-bands `visited` buys 0.56 GB and changes output (a peak claimed by a neighbouring band mid-flight becomes `++d.missed` at :985, which can terminate or invalidate a trace). Both are refused on ratio, not on correctness.

**Ponytail note:** rows 2, 3 and 4 all *delete* — a struct field, a P-sized allocation, an `assign()` first-touch pass, and (with `got` carrying pairs) the per-trace `pts` vector at :1005. Row 1 adds four lines. Row 5 adds a counting pre-pass and two `PeakSlab` fields.

---

## 3. Lifetimes and pooling

### 3a. The slab and prep are held across the whole EPD stage for nothing — the headline item

`splitIntegerTraces` (:1052) takes `(vector<Trace>&, double, TraceStore&, int)`. **It has no `PeakSlab` parameter and reads no slab field.** Its only inputs are `wst.bins` (:1077), `wst.b` (:1078, :1104) and the trace spans. `tprep` is last read inside `detectTracesInteger_`, and the band `taskloop` at :3274 carries an implicit taskgroup, so every band task has joined by :3282.

So at line **3292** — before `splitIntegerTraces` at :3296 — the following are provably dead:

| structure | bytes | freed today at |
|---|---|---|
| `wslab.tof / inten / imq / b / frame_off / rt_index` | 12.0 GB | :3300-3303, **after** split |
| `tprep.frame_of / member / frame_live / order` | ≤ 10.8 GB | :3304 (end of block), **after** split |

They are held across exactly the stage where `mts` + `split` + the `MassTrace` payloads peak. The fix is to move :3300-3303 above :3296 and wrap :3260-3292 in a block so `tprep` and a locally-scoped `g_prep` die there. **Net lines added: zero (one brace pair, four lines moved).** Output-neutral by construction — the type system enforces it, because `splitIntegerTraces` cannot name either object.

The ledger will report *no change* unless `MemGuard` gains a release, so this must ship with §3d.

### 3b. The valley-splitting list

`std::list<Peak2D> pk` (:1072) is declared **inside** the per-trace loop and destroyed each iteration at :1080; `MassTrace mt(pk)` (:1081) copies it into a `std::vector<Peak2D>`. At most one list per running thread exists — **1–5 MB**, not 6.0 GB.

`MemGuard g_list(MEM_EPD_LIST, list_bytes)` is constructed at :1086, *after* the loop, with `list_bytes` summed over the whole chunk (:1071). **The 6.0 GB "EPD list<Peak2D>" line is a phantom.** Its numeric coincidence with reality is an accident: 48 B/pt happens to approximate the 32 B/pt that `mts` really holds plus the 32 B/pt that `split` holds beside it.

**A pool is the wrong answer** and should be recorded as refused. `MassTrace` (OpenMS `MassTrace.h`) exposes only two constructors — no setter, no `reserve`, no `clear` for `trace_peaks_` — so a recycled `MassTrace` allocates a fresh internal vector anyway; a pool saves the header and nothing else. And the only poolable object, the point buffer, is just a variable that outlives the per-trace loop.

### 3c. The valley-splitting objects — streaming is the right answer

What must be alive: **one input `MassTrace` and its immediate split output, per thread.** What is alive: the whole chunk's `mts` (24 B/pt `trace_peaks_` + 8 B/pt `smoothed_intensities_`, written in place by `ElutionPeakDetection::smoothData`) *and* the whole chunk's `split` (another 32 B/pt), both through the `toTrace` loop at :1095-1107.

Two independent derivations of the concurrent cost agree:
* from the (phantom-anchored) ledger: 1.25e8 points × 64 B ≈ **8.0 GB** payload + ~3.6 GB headers ≈ 11.6 GB;
* from the chunking directly: chunk = n_in/48 ≈ 156k traces × ~10.7 pts × 64 B ≈ 107 MB per live chunk × ~100 live chunks ≈ **10.7 GB**.

**Take 10–11.6 GB as the range.** Note that one skeptic's 1.0 GB "give-back" objection — that streaming loses today's per-trace input free — **no longer applies**: `Trace::freeProfile()` is a no-op on the arena (:267) and the parents' arena is released wholesale at :1110-1111 after the taskloop.

The change: replace the two-phase chunk body with a per-trace loop calling the single-trace overload `epd.detectPeaks(mt, one)`, converting each element of `one` with `toTrace` straight into `parts[ci]`, using **one reused `vector<Peak2D>` scratch declared *inside* the taskloop body** (the pragma at :1060 is `default(shared)`; a scratch above it is a data race). `splitValleys` (:1025) then has no caller and should be deleted. **Net: a function removed.**

### 3d. What makes a spectrum FINAL, and what enforces order

**Content** is final the moment `assembleOne_` returns (:3441) — peaks, precursor m/z, charge, both isolation offsets, drift time, `spx_guessed`, `spx_n_isotopes`. After that only `merge:rt_window` (:3642) and `consolidate:delta_rt` (:3771) can mutate it; both default to 0 = OFF and both are documented FALSIFIED on dataset D in their own option text.

**Identity is NOT final until the total order is known.** `assembleFromList_` never calls `setNativeID`, so OpenMS's `MzMLHandler` sets `renew_native_ids = true` and writes `id="spectrum=<index>"` — **every spectrum is named by its rank in the canonical sort**.

Order is enforced by exactly two things:
1. every result is written to an **index-addressed slot** — `win_out[wi]` (:3352), `pslot[pi-plo]` (:3443), `per[b]` (:3279) — never appended in completion order (the comment at :3156-3157 states this as an invariant);
2. the **canonical total-order sort** at :3595-3610 over (RT, precursor m/z, charge, size, fragment sequence), which requires every spectrum resident simultaneously.

RT is a **global** axis and every window emits across the whole gradient, so completion-order writing cannot reproduce canonical order. That is the hard constraint on everything in §4.

**Emitted-spectra lifetime.** `MEM_SPECTRA` is `memAdd`-ed once per window at :3450 and **never** subtracted, so its 5.0 GB is the *end-of-loop* value, not a peak-concurrent one — and the `memAdd` sits inside the share-all branch only, so the apportion (:3400) and rp_max (:3511) branches are unaccounted. Nobody knows what the spectra buffer held at the RSS peak. The one mechanism by which held spectra genuinely cost something is the admission gate at :3178, which reads live free RAM: dead output in `win_out` directly reduces the budget lent to later windows.

The two free fixes here: `out.reserve(frags.size())` before the emplace loop (**0.15–0.2 GB**; the 0.4 GB claim was refuted on tcmalloc size-class grounds — `max_fragments` = 500 → 8000 B and 8192 B land in the same 8192 class), and `out_exp.setSpectra(std::move(all_out))` replacing the loop at :3819 (**0 GB at peak; 1.27 GiB steady / 1.90 GiB transient during WRITE**, ~0.5 s, one line replacing two).

---

## 4. Writing during computation

**Neither, and it should not be built.**

**Memory: no.** The saving is `MEM_SPECTRA` at the RSS-peak instant, which is unmeasured and structurally likely to be small — the peak sits mid-loop where slab, prep, traces and EPD coexist, while spectra accumulate monotonically toward the end. Blocks can only be released in ascending key order, and the master admits windows **heaviest-first** (`stable_sort` by compact bytes), which is uncorrelated with the isolation-window key: expected freed prefix at the peak is ~0.1–1.0 GB. Those are 2–8 KB allocations, below the mmap threshold, so they return to arena free lists, not to the OS.

**Runtime: no.** WRITE is one phase of a 413 s run; overlapping it with a 343 s window loop can at best hide it, and it cannot be hidden because of the ordering constraint below.

**Three hard blockers, any one fatal:**
1. `<spectrumList count="N">` is written from `spectra_expected_` at the **first** `consumeSpectrum`. N is unknown until the last window finishes. The harness enforces this (`hash_and_validate_mzml` fails with G-B11 when the count attribute disagrees with the spectrum count).
2. `renew_native_ids` is hardcoded `false` in `MSDataWritingConsumer`, and our spectra carry no nativeID → `id=""` on every spectrum. Same class as the recorded MSFragger `Scans=0` blocker.
3. `merge:rt_window` and `consolidate:delta_rt` both random-access a fully materialised `all_out`; streaming silently turns them into no-ops. They default off, so the dataset D gate would not catch it.

**The ordering constraint, stated exactly:** a spectrum's file position is `rank(RT, precursor m/z, charge, size, fragment sequence)` over the *whole run*, and its nativeID is that rank. No per-window blocking scheme can compute it, because every window emits across the entire RT range. Streaming would require either abandoning canonical order (which changes every nativeID and is an output change requiring the both-engines gate) or a spill-and-merge design — 150 lines of temp-file machinery, a serialisation format and a k-way merge — to relocate ~1 GB in a 105 GB process.

---

## 5. The occupancy gap

**Leading hypothesis: there is no gap.** Occupancy is CPU-seconds ÷ wall-seconds by definition. On dataset A: integer 31,033 CPU-s at 66×, OpenMS 41,530 CPU-s at 84×. Work ratio 31,033/41,530 = **0.747**; occupancy ratio 66/84 = **0.786**. They agree to within 5 %. A path that does 25 % less work in comparable wall time *must* show ~25 % lower occupancy. The unexplained residual is ~5 %, not 20 %.

The direction also argues against every bandwidth explanation offered: a thread stalled on a cache miss is still in `R` and still counts as occupied, so memory pressure *inflates* CPU-seconds. The integer path has fewer.

Three candidate mechanisms were sized and all fail by two to four orders of magnitude against the ~7,850 thread-seconds the residual would represent: `visited` first-touch (15.1 GB → 3.69e6 minor faults ≈ 5–10 CPU-s), `visited` mmap/munmap churn (288 events for the whole run), and the `b_of_rt`-class per-point lookups (~6–20 CPU-s). The one mechanism with a *measured* per-call cost — `SavitzkyGolayFilter` coefficient reconstruction inside `ElutionPeakDetection::smoothData`, measured on a Release OpenMS build at 16.8 µs against 18.0 µs for a whole 12-point `detectPeaks` — is real but is paid by **both** detector paths, so it cannot be the differential either.

**The one measurement that settles it costs nothing and needs no rebuild.** Attach a 1 Hz thread-state sampler to the next scheduled dataset D run of each arm:

```
while :; do date +%s; awk '{print $3}' /proc/$PID/task/*/stat | sort | uniq -c; sleep 1; done
```

* integer shows ~30 threads in `S` where OpenMS shows ~15 → threads are parked; look at the 20 ms admission sleep (:3186) and libgomp taskgroup idling.
* integer shows extra `D` → page-fault / mmap path; the bitset change is then also a throughput fix.
* both show ~100 `R` while occupancy still reads 68× → the metric is arithmetic, the gap is not real, and this line of inquiry closes.

**Do not** spend a run on `peak_visited_bytes / P_window`. It is already computable: 24 × (4.9e9 / 1.26e9) ≈ 93 of 100 threads, which is at the arithmetic ceiling (band tasks hold no scheduling point while holding `visited`) and therefore carries no information.

---

## 6. Ranked plan

Ordered by (measured or well-bounded saving) ÷ (risk × code added). **B** = output-neutral by construction, verifiable by the per-window `[det] frag traces … digest=` line at :3317. **E** = requires the both-engines rule (Sage + MSFragger + entrapment).

| # | change | GB saved | runtime | output | experiment |
|---|---|---|---|---|---|
| 0 | **Verify the landed trace redesign** (no code) | re-takes the ledger | — | **B** | dataset A, HEAD vs `4bd81c2`, `SPEXTRACT_DET=1` + `SPEXTRACT_MEM_LEDGER=1`. Falsifier: any per-window digest differs, **or** `traces` does not fall 35.6 → < 12 GB. |
| 1 | **Free slab + prep before `splitIntegerTraces`** — move :3300-3303 above :3296; scope `tprep`+`g_prep` to a block ending :3292 | ceiling 22.8; realistic **3–8** (× fraction of windows in EPD at the peak) | neutral / mildly positive (gate reads real free RAM) | **B** — `splitIntegerTraces` cannot name either object | ledger delta on `MEM_SLAB` / `MEM_PREP` + process peak RSS. Falsifier: RSS and the simultaneous-peak line both unmoved ⇒ no window is ever in EPD while another peaks. |
| 2 | **`visited` → bitset** (:894, 4 lines) | **4.29** (both skeptics exact) | neutral | **B** | same run as #1 (different ledger line). Falsifier: any digest difference = indexing bug. |
| 3 | **Ledger fix** (~10 lines) | 0 | neutral | **B** | prerequisite for every byte claim below — see §7 E1. |
| 4 | **Stream valley splitting** + reused `vector` scratch + delete `splitValleys` | **10 – 11.6** | faster or neutral | **B** (per-trace call order equals the patched nested serial loop) | isolated build. Falsifier: per-window digest differs ⇒ `detectElutionPeaks_` has cross-trace state. |
| 5 | **Drop `frame_of`** + fold `pts` into `got` as pairs | **1.6 – 4.8** | neutral to faster | **B** | `f0` **must** be `upper_bound(...)−1`, never `lower_bound` (empty frames); :1001 becomes `visited[p.second]`. |
| 6 | **Drop `member`** | **0.4 – 1.2** | ambiguous — adds a 4 B stream to `best()`'s walk | **B** | separate run from #5 or the timer is uninterpretable. Line 869 keeps both conjuncts. |
| 7 | **Exact-reserve `order`** + log seeds/P | **0 – 2.4** (exp. ~1.2) | neutral | **B** | one line + one counter; rides in E1. |
| 8 | **Compact slab to members** | **(1−m) × ~18 GB**; m = 0.6 → ~7, m = 0.9 → ~1.8 | faster (shorter candidate walk) | **B** *only if* band edges come from the raw tof range and raw occupancy is stored; otherwise **E** | run **only** if E1 reports m < 0.8. Assert raw `tlo`/`thi` unchanged before/after. |
| 9 | **SG coefficient cache** (OpenMS patch, ~6 lines) | 0 | measured 87–95 % of per-trace `detectPeaks`; end-to-end share unmeasured | **B** (pure function of frame_length, order) | needs a split-stage timer first. Falsifier: split stage < 3 % of loop CPU. |
| 10 | **`setSpectra(std::move(all_out))`** (:3819) | 0 at peak; 1.27 GiB during WRITE | ~0.5 s | **B** | read the `[perf]` table: ASSEMBLE vanishes, WRITE `rss_end_mb` drops. |
| 11 | **`out.reserve(frags.size())`** | **0.15 – 0.2** | negligible | **B** | ledger only, once #3 reports capacity as well as size. |
| 12 | **Derive `Trace::intensity`** (and possibly `b`) | **1.6** (3.2 with `b`) | neutral | **B** for `intensity`; `b` needs its own check (OpenMS path stores 0.0) | only after E0 gives the real trace count. |

Nothing in this table requires the both-engines rule as written. **That is the point of the ordering:** every item is either a lifetime move, a representation change, or a deletion, and each is falsifiable by a digest comparison that costs one dataset A run rather than a Sage + MSFragger + entrapment cycle. Item 8 is the single item that can slip into **E** if its raw-edge invariant is not held, which is why it is gated behind a cheap measurement of *m*.

---

## 7. Experiment schedule

Runs are dataset A unless stated. Builds that share a run are ones whose effects land on **different ledger lines**, so the ledger attributes them even though process RSS does not.

**E0 — the baseline nobody has (isolated, must be first).**
Two runs of the *unmodified* HEAD plus one of `4bd81c2`, with `SPEXTRACT_DET=1` and `SPEXTRACT_MEM_LEDGER=1`. The two HEAD runs establish whether the per-window `[det]` digest is reproducible at fixed threads — the project record says the tool is **not** digest-reproducible end-to-end (~1e-4 factor; the standing gate is peptide-set overlap), and the `[det]` line is emitted upstream of that. If it is stable, digests gate everything below; if not, every **B** item falls back to peptide sets and the plan gets three times more expensive. Attach the §5 thread-state sampler; it is free.
*Falsifier:* any window digest differs between HEAD and `4bd81c2` ⇒ the redesign changed output and must be fixed before anything else. Or `traces` does not fall below 12 GB ⇒ the redesign's own arithmetic is wrong.

**E1 — instrument + two lifetime/representation fixes + one counter (one build, one run).**
Ledger fix: `MemGuard::release()`; `memSub(MEM_TRACES, …)` at window end (today :3299/:3312 add and nothing subtracts, so that line accumulates over the run and is not a peak); charge the real `MassTrace` payload (24 + 8 B/pt) once per chunk instead of the phantom list; hoist the `MEM_SPECTRA` add out of the share-all branch; sum `capacity` as well as `size` for spectra; rename the stale row labels at :718-720. Plus #1 (free slab+prep), #2 (bitset) and #7's seed counter.
*Attribution:* #1 → `MEM_SLAB` + `MEM_PREP`; #2 → `MEM_VISITED`; the counter → a log line. Disjoint.
*Falsifier:* `MEM_EPD_LIST` does not collapse to single-digit MB ⇒ the loop-scoping analysis of `pk` is wrong. `MEM_PREP` does not land in [0, 10.8] GB ⇒ a second accounting defect remains and no byte claim below is safe.

**E2 — stream valley splitting (isolated).**
It moves `MEM_EPD_LIST` and `MEM_EPD_MT` wholesale, so it cannot share a run with anything that touches those lines, and it needs the corrected ledger from E1 to be readable at all.
*Falsifier:* any per-window digest difference; or process peak RSS falls by materially less than the corrected `MEM_EPD_*` delta, meaning the bytes moved rather than disappeared.

**E3 — drop `frame_of` (isolated).** *Falsifier:* digest difference (⇒ the `upper_bound`−1 or the `got`-pair transcription is wrong), or `MEM_PREP` does not fall by 4 B/peak.

**E4 — drop `member` (isolated, because its runtime direction is genuinely unknown).** *Falsifier:* digest difference; or the split-stage/prefix timer shows `best()` got slower by more than the removed array is worth.

**E5 — compact slab to members (isolated; run only if E1's counter gives m < 0.8).** *Falsifier:* raw `tlo`/`thi` differ before/after, or the digest moves. Either means the band partition shifted and the change is an output change, not a compaction.

**E6 — SG coefficient cache (isolated; run only if a split-stage timer says the split stage exceeds ~3 % of window-loop CPU).** Note the missing instrument: there is **no** timer around `splitIntegerTraces` — :3296 falls inside `t_trace` (:3226 `_t2` → :3343 `_t3`), pooled with `prepareTracing`, all band tasks and the split, and summed over overlapping window tasks. Adding that timer is a prerequisite, not part of the experiment.
*Falsifier:* the split stage is under 3 % of loop CPU ⇒ the whole item buys < 3 % of wall and is not worth a carried OpenMS patch.

**E7 — hygiene bundle (one build): `setSpectra` move + `out.reserve`.** Read the `[perf]` phase table only; no search engine.
*Falsifier:* the ASSEMBLE row does not vanish, or WRITE's `rss_end_mb` does not drop by ~1.3 GiB scaled to dataset A.

**Final gate.** One dataset D run of the accumulated build against the dataset D baseline, on **both** engines plus entrapment — not because any single item requires it, but because the project rule is that the baseline is refreshed with both engines at the end of a cycle, and because the combined build has never been measured on dataset D.

---

## 8. Rejected, and why

**Trace encoding (all obsoleted by the landed redesign):** `p2-arena-index`, `p3-arena-csr`, `both-arenas`, `p3-csr-inplace-grid`, `drop-xt-only`, `soa-trace` — the structures they act on (`xi/xv/xt`, `FragGrid.flat`) no longer exist; each additionally carried a use-after-free on the `split_valleys = 0` early return, a semantic break (`xi` is a *global* RT-axis index, out-of-bounds as a slab peak index), or a peak-RSS regression from pinning the arena past `buildFragStats`.

**`no-arena`** — recorded as decided: after the redesign there are ~1e4 live allocations where there were ~3.6e8, so a bump allocator has ~0.2 MB of per-allocation overhead left to remove, and `Trace` crosses three task boundaries so the arena is not scope-bounded.

**`pack32`** — 10-bit per-window `imq` dictionary cannot hold: a single frame already has ~1,000 scans and TIMS calibration is frame-dependent, so distinct `imq` per window is 1e3–1e5, and the fail-closed check would abort the benchmark.

**`inten-exact` / `inten-log16`** — slab intensities are `PeakPickerIM` cluster *sums*, so u16 is not reachable; log16 collapses distinct integer intensities above ~3,560 to one code, which permutes the greedy seed order (`visited` claim-first) and the strict apex test — a discrete ownership change, not a 0.03 % perturbation.

**`slab-in-place`** — 0 GB by its own arithmetic; pins each window's peaks in thousands of small arena chunks through the trace stage instead of one mmap'd block, plausibly *raising* fragmentation, and re-indexes the hottest loop.

**`exact-reserve-order` as a 5 GB claim** — kept as item #7 at 0–2.4 GB; the 5 GB version was derived from the impossible 14.3 GB prep line.

**`band-bucket-order` / `seed-bucket` / `quantile-bands`** — all three derive a band index by division while the taskloop tests truncated edges (`tlo + (TofIdx)(b*wdt)`, :3277-3278); the two are not inverses, so boundary seeds change bands, change `visited` ownership, and change output. Quantile edges additionally place every boundary in the *densest* region, multiplying the known 0.1 %-of-peptides cross-band duplicate term.

**`fuse-sweeps`** — removes ~262 MB/window of prefetched streaming reads (~30 ms wall) while leaving the seed sort, which does ~3e8 cache-missing indirect comparisons in the same function.

**`bitset-shared` / `reconcile-at-end` / `collide-count`** — sharing `visited` across bands turns a taken peak into `++d.missed` (:985), which can terminate a direction or fail `min_sample_rate`, i.e. it *deletes* traces near boundaries; reconciliation cannot undo it because `visited` steers growth, not just membership; the collision counter measures a first-order rate that does not bound the cascade, and the counterfactual is already available for free via `-perf:trace_bands 1`.

**`vector-not-list` / `pool-reject` / `vector-scratch` as standalone items** — subsumed into the streaming change (#4); alone they are ~0.1–0.5 % of runtime and 0 GB, and a `MassTrace` pool saves nothing because `trace_peaks_` has no refill path.

**`finer-chunks` (nchunk 48 → 192)** — superseded by streaming, which removes the chunk-lifetime residency entirely rather than dividing it.

**`native-split`** — reimplements Savitzky-Golay + `findLocalExtrema` + the sub-trace loop (~250 lines, forked from upstream) for ≤ 0 GB beyond streaming, and gives up the currently-verified cross-node bit-determinism.

**`compact-emission-buffer` / `spill-and-merge` / `stream-window-blocks`** — see §4: all three break `<spectrumList count>`, native-ID renewal, and merge/consolidate, for a saving whose value at the RSS peak is unmeasured and structurally small.

**`bofrt-array`** — `b_of_rt` no longer exists; `splitIntegerTraces` reads `wst.b` directly.

**`tof-minmax-in-toslab`** — correct and free, worth ~0.5–1.3 CPU-s of 31,033; fold it in if someone is already editing `toSlab` for item #8, but attach no runtime claim and do not spend a run.

**`admit-in-master`** — `need` computed in the master while `Release` lives in the task makes `need` a shared reference to an overwritten stack slot; `inflight` then drifts and the `inflight == 0` escape hangs the master.

**`measure-first` (E1/E2/E3 as proposed)** — right instinct, wrong instruments: E1's counters mirror the unbalanced `MEM_TRACES`; E2's allocator A/B is confounded because `availableBytes_()` returns MemAvailable **plus our own RSS**, so a lower-overhead allocator admits more windows and RSS returns to the same ceiling; E3's rt_index collision is already impossible because `rtAxis()` is `sort`+`unique`d.

**`allocator-of-record` / `malloc-info`** — the harness LD_PRELOADs `tcmalloc_minimal` by default, the glibc↔tcmalloc peak A/B is already on record (101.8 → 94.3 GB), and `malloc_info` at exit reports an arena state uncorrelated with the mid-loop peak.

**`ledger-confound` / `stage-timers` / `rss-timeline` / `band-phase-concurrency` / `mem-peak-snapshot` / `ledger-shrink` as standalone** — the useful half of each (a `MemGuard` release, a split-stage timer, a per-window RSS log line in the existing `#pragma omp critical`) is folded into E1; the standalone versions each cost a run to produce a number that is already derivable, or introduce a race on a static snapshot array.