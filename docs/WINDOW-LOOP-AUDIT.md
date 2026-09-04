# Window-loop audit — where the 3,805 s goes

Scope: `src/spextract.cpp:1526-1720` (the `#pragma omp parallel for` over the 24 diaPASEF
isolation tiles) and everything it calls, including OpenMS
(`<OpenMS>/src/openms/`).

Every line reference below was checked against the working tree at `787441a`
(`src/spextract.cpp`, 1,932 lines) and against the OpenMS source on this machine. Claims
are labelled **MEASURED** (a number that exists in a log) or **INFERRED** (derived from code
structure plus arithmetic). The distinction is the point of this document: almost everything
previously asserted about the inside of the window loop is inferred, and one widely repeated
"measurement" is not a measurement at all.

Quality constraint: the dataset D baseline is **10,072 peptides** (`docs/dataset D-BASELINE.md:11`), noise
floor 0.06% same-node, "differences below ~0.2% mean nothing". A speedup that costs peptides is
not a speedup.

---

## 1. WHERE THE LOOP'S TIME GOES

### 1.1 MEASURED

dataset D (9.8 GB, smallest sample), 120 threads, tcmalloc already `LD_PRELOAD`ed:

| phase | time | share | source |
|---|---:|---:|---|
| load + pick + split | 32 s | 0.8% | `clk_()` at `:1257`, `:1314` |
| MS1 tracing + inference | 330 s | 7.9% | `clk_()` at `:1314` → `:1458` |
| **window loop** | **3,805 s** | **91.3%** | `clk_()` at `:1458` → `:1736` |
| total | 4,167 s | | |

| resource | value |
|---|---:|
| wall | 4,428 s (01:13:48) |
| CPU total | 64,521 s (17:55:21) |
| user | 36,456 s (10:07:36) |
| **system** | **28,065 s (07:47:45) = 43.5% of CPU** |
| effective cores (CPU/wall) | **14.57** of 120 |

the reference implementation does the same job in 457 s → we are **9.1×** slower.

### 1.2 The window loop is ONE undivided number — this is the central measurement gap

`phase_clock_()` is defined at `:77`, wrapped by `clk_()` at `:82`, and called at exactly five
sites: `:1222`, `:1257`, `:1314`, `:1458`, `:1736`. **There is no timer anywhere between `:1526`
and `:1720`.** "WINDOW LOOP 3,805 s" is therefore a single opaque interval covering nine distinct
stages. Every attribution below it is inferred.

### 1.3 The loop's structure — MEASURED from code, not timed

- `n_threads = omp_get_max_threads()` = 120 (`:1483`; `TOPPBase.cpp:86` applies `-threads 120`,
  `harness/bench.py:110`).
- `n_bands = ceil(120 / 24) = 5` (`:1489-1492`), `perf:trace_bands` default 0 → auto.
- `n_conc = min(24, 120/5) = 24` (`:1497`) — i.e. **every window is admitted at once**.
- `#pragma omp parallel for schedule(dynamic) num_threads(n_conc)` over 24 iterations with 24
  threads: **`schedule(dynamic)` is inert.** One window per thread, no work stealing, no
  rebalancing. The loop's wall time is the *slowest single window*, not the mean.
- Admission (`:1535-1546`) is a `for(;;)` poll with `sleep_for(200 ms)` inside
  `#pragma omp critical(admit)`, re-reading `/proc/meminfo` each pass. 24 × ~5.8 GB ≈ 139 GB
  against the ~150 GB observed resident ⇒ all 24 were admitted; throttling is not the limiter.

### 1.4 What runs at what width inside one window — INFERRED from code

| lines | stage | threads |
|---|---|---|
| `:1549-1551` | `materializeWindow` + `sortSpectra` (~56 M peaks, ~1.1 GB) | 1 |
| `:633-649` | band partition, `bands` full passes over every peak, both output vectors unreserved | 1 |
| `:652-666` | `mtd.run` per band (`num_threads(bands)`) | **5** |
| `:684-718` | ElutionPeakDetection valley split, 32 batches of 250 k | **nested, see §2** |
| `:721-737` | `[merged-trace]` RT-span instrumentation, 64 MB alloc + 8 M-element sort | 1 |
| `:738-741` | `toTrace` × ~8 M (each with its own `sort`) | 1 |
| `:1568-1573` | canonical sort of 8 M × 56 B = 448 MB | 1 |
| `:1579` | `buildFragGrid` (sorts ~56 M doubles ≈ 450 MB to recover ~1,343 uniques) | 1 |
| `:1591-1712` | scoring, ~40 k precursors | 1 |

### 1.5 Two arithmetic facts that constrain any attribution

**(a) Band parallelism contributes almost nothing to the time integral.** If `mtd.run` at `:652`
dominated the loop, the run would show ≈ 24 × 5 = 120 effective cores. It shows **14.57**. Even
crediting *all* measured CPU to the loop gives 64,521 / 3,805 = **17.0 concurrent** — below 24,
i.e. below one busy thread per window. The wall sits in code running at ~one thread per window.

**(b) Scoring cannot be the bottleneck.** `scoreCandidates_` is one call per precursor on the
default path (`:1648` → `:1168`); the apportion (`:1616`/`:1627`) and `rp_max` (`:1663`) paths are
default-off (`:496-497`) and both are on the falsified list (`docs/dataset D-BASELINE.md:35-36`). At
~2.2 ms/precursor × ~40 k precursors that is **~90-170 s of a ~3,805 s per-window wall = 2-4%**.
Independently: the claimed ~22 TB of gate-scan traffic over 3,805 s is 5.8 GB/s aggregate on a
node that sustains ≥100 GB/s, capping the whole scan at ~12% even in the worst case.

**Conclusion.** The time is in the single-threaded stages of the window body and in the one stage
that is neither single-threaded nor sanely parallel — the nested EPD region. §2 argues the latter.

### 1.6 The instrumentation that must exist before any further ranking

Add `phase_clock_()` deltas around `:1549`, `:1555` (split into MTD / EPD), `:1568`, `:1579`,
`:1591`, plus per-window wall, reduced in the `#pragma omp critical` that already exists at
`:1557`. Six calls. Until they exist, everything in §1.4 is a hypothesis with a plausible
mechanism, not a budget.

---

## 2. THE 43.5% KERNEL TIME

### 2.1 The single best-supported explanation: a 24×-oversubscribed nested OpenMP region

Six facts, each verified at a cited line, compose into one mechanism:

1. **MS2 valley splitting is ON.** `trace:ms2_split_valleys` registers a default of `0.0`
   (`:514`), but `bench/run_open.sh:42` and `docs/dataset D-BASELINE.md:21-22` pass **7.0**. So
   `:684` (`split_valleys > 0.0 && !mts.empty()`) is true on every window and
   `ElutionPeakDetection::detectPeaks` runs *inside* the window body.

2. **Nested parallelism is enabled.** `:1502` executes `omp_set_max_active_levels(2)` whenever
   `n_bands > 1`, which is the default configuration. Level-2 regions are therefore **active**,
   not serialized.

3. **The EPD call site is at nesting level 2.** The band region closes at `:666`;
   `epd.detectPeaks(chunk, part)` is at `:708`, inside the window loop but outside the band loop.

4. **EPD's parallel region has no thread cap.**
   `<OpenMS>/src/openms/source/FEATUREFINDER/ElutionPeakDetection.cpp:326-329`:
   ```cpp
   #ifdef _OPENMP
   #pragma omp parallel for
   #endif
       for (SignedSize i = 0; i < (SignedSize) mt_vec.size(); ++i)
   ```
   No `num_threads` clause ⇒ team size = `nthreads-var` = **120**.

5. **⇒ up to 24 × 120 = 2,880 live OS threads on 120 cores.** Entered 32× per window
   (`BATCH = 250000` at `:696` over ~8 M traces) = **768 fork/join cycles of a 120-thread team**,
   each with a full barrier, while 23 other windows are doing the same thing.

6. **Every result leaves through a program-global lock.**
   `ElutionPeakDetection.cpp:459-463` and `:544-548` are
   `#pragma omp critical (OPENMS_ElutionPeakDetection_mtraces)` — a *named* critical, which in
   OpenMP is one lock **per program**, not per region. Roughly 8 M emits/window × 24 ≈ **192 M
   contended acquisitions**, each doing a heap-allocating deep copy (`push_back` of a `MassTrace`
   that has no move constructor) plus periodic vector reallocation **inside the lock**.

**24× oversubscription is the term that generates system time.** 2,880 runnable threads on 120
cores means every OpenMP barrier and every `critical` exhausts its spin budget before the holder
is even scheduled, and falls through to `futex(FUTEX_WAIT)`/`sched_yield`. That is charged as
kernel time and shows up as involuntary context switches, which is exactly the observed 43.5%
and the 14.57/120 efficiency.

### 2.2 Two aggravating factors at the same call site

- **`ElutionPeakDetection.cpp:334-336`**: `#pragma omp atomic ++progress` executes **per trace**,
  outside the `IF_MASTERTHREAD` guard — one contended RMW on a single shared cache line by a
  120-thread team, ~8 M times per window.
- **`ElutionPeakDetection`'s constructor calls `setLogType(CMD)`
  (`ElutionPeakDetection.cpp:43`)** and `spextract.cpp:686` never silences it — unlike
  `mtd.setLogType(ProgressLogger::NONE)` at `:584` ("quiet + thread-safe (called from the parallel
  window loop)") and `m2.setLogType` at `:655`. `startProgress`/`endProgress` mutate
  `ProgressLogger::recursion_depth_`, a **non-atomic `static int`**
  (`ProgressLogger.cpp:125`, `++` at `:238`, `--` at `:268`), from 24 concurrent window threads.
  That is a genuine data race, independent of how much time the logger costs.

### 2.3 What this explanation does NOT claim

**The named critical alone does not account for 43.5%.** Measured aggregate lock throughput with
a `MassTrace`-sized payload is ~16.6 M emits/s at 16 threads. 192 M emits ⇒ ~12 s of fully
serialized wall and ~170-550 CPU-s of contention — **0.3-0.9% of 64,521 s, not 43.5%**. The lock
is the *funnel*; the 2,880 threads queuing at it are the cost.

Page faults do not account for it either. The three largest fault sources (band partition at
`:633-649` ≈ 550 k/window, `materializeWindow` ≈ 280 k, `buildFragGrid` ≈ 165-230 k) total under
1 M minor faults per window ≈ 24 M for the run ≈ **~24 s** at 1 µs each — three orders of magnitude
short of 6.3 cores permanently in the kernel. THP on `madvise` with no `madvise` call anywhere in
the tool is a real observation and worth a `MADV_HUGEPAGE` experiment, but it is a hypothesis, not
the explanation.

### 2.4 The falsifying measurements — run these first

1. **`grep Threads /proc/$(pgrep -f spextract)/status` during the window loop.**
   ~120 ⇒ this whole section is wrong. ~2,900 ⇒ confirmed. One command, decides everything.
   Note the log line at `:1504-1507` is **not** evidence either way: it prints
   `String(n_conc * n_bands) + " live threads"`, the tool's own arithmetic over the band loop,
   which does not know EPD exists.
2. **`OMP_NUM_THREADS=120,1`** — the correct zero-code A/B. The list form sets `nthreads-var`
   per level; both explicit `num_threads` clauses (`n_conc` at `:1526`, `bands` at `:652`)
   override it, so *only* EPD's bare `parallel for` is capped.
   **`OMP_MAX_ACTIVE_LEVELS=1` does nothing** — `:1502` overrides the environment ICV at runtime.
   Verified: `OMP_MAX_ACTIVE_LEVELS=1` leaves the EPD team at 16/16 and total live threads
   unchanged; `OMP_NUM_THREADS=16,1` drops it to 1.
3. `perf stat -e context-switches,minor-faults` and `perf top -g -K` / `-U` to split kernel time
   between futex and fault handling.

---

## 3. OUTPUT-IDENTICAL FIXES, ranked by (expected speedup)/(risk)

"Output-identical" here means: the emitted mzML is byte-identical, **or** identical to within the
run-to-run nondeterminism the current code already has. That caveat is necessary and is stated per
fix: the canonical comparator at `:1568-1573` sorts on `(im, mz, rt, intensity)` with no index
tie-break, `std::sort` is not stable, and today's trace insertion order is already
thread-schedule-dependent through the EPD critical — so the 10,072-peptide baseline is **not
bit-reproducible today**.

### OI-1 — Cap the nested EPD team. `src/spextract.cpp:686-716`. **← IMPLEMENTED**

**Change.** Immediately before the batch loop that calls `epd.detectPeaks`, and only when already
inside a parallel region, set the *calling task's* `nthreads-var` to the same inner budget the
band loop was given (`bands`); restore it after. `omp_set_num_threads` writes only the calling
task's ICV, so it cannot disturb the explicit `num_threads(bands)` clause at `:652`, and the
`omp_get_level() > 0` guard leaves the **MS1** call (`:1345`, level 0, `ms1_split_valleys` default
7.0) at the full 120 threads where that parallelism is legitimate. Also silence EPD's
`ProgressLogger`, matching `mtd` two lines up.

**Expected gain.** Live threads during the EPD stage drop from ≤ 2,880 to `n_conc × bands` ≈ 120,
i.e. exactly the thread budget. If §2 is right this removes the bulk of the 28,065 s of system
time. A conservative reading — kernel time falling from 43.5% to the ~10% a well-behaved OpenMP
program shows — implies the loop's 3,805 s falls to roughly 1,600-2,400 s, a **1.6-2.4×** on the
loop and ~1.5-2.1× overall. If §2 is wrong the change is a no-op, not a regression.

**Why it cannot change output.** It changes only how many threads execute a loop whose iterations
are independent (`detectElutionPeaks_(mt_vec[i], single_mtraces)`, each reading a distinct
`mt_vec[i]`). The *set* of emitted sub-traces is thread-count-invariant. The *order* in which they
land in `single_mtraces` is today decided by which thread wins the critical — i.e. already
nondeterministic — and becomes closer to deterministic index order under the cap. Downstream, the
MS2 traces are canonically re-sorted at `:1568-1573` before anything reads them. Net:
**set-identical, and strictly more deterministic than today.**

**Falsifying measurement.**
- `grep Threads /proc/<pid>/status` mid-loop still shows ~2,900 ⇒ the cap did not take (libgomp
  and libomp both honour the per-task ICV, but verify).
- Window-loop wall does not move and kernel share stays ≥ 40% ⇒ §2's mechanism is wrong; the time
  is elsewhere and the phase timers of §1.6 become mandatory.
- dataset D peptides move by more than the 0.06% noise floor ⇒ the order argument above is wrong and the
  change must be reverted.

### OI-2 — Give `MassTrace` a move constructor. `OpenMS/.../KERNEL/MassTrace.h:66-72`

**Change.** Two lines: `MassTrace(MassTrace&&) = default;` and
`MassTrace& operator=(MassTrace&&) = default;`. `MassTrace.h:66,69,72` user-declare the destructor,
copy constructor and copy assignment (`= default` still counts as user-declared), which by
`[class.copy.ctor]/8` **suppresses the implicit move members**. `MSSpectrum.h:173/182` and
`MSExperiment.h:107/113` do declare theirs — `MassTrace` is the outlier, so this is an OpenMS
defect, not a design choice.

**Expected gain.** Measured: `std::move(mt)` into a `vector<MassTrace>` is currently a deep copy
(`src.getSize()` stays 7 after the move); `push_back` costs 5.277 allocations/trace unreserved,
2.000 reserved. With the patch, 0.001/trace. Every `std::move` in `detectTraces_` (`:664`, `:669`,
`:703`, `:710`) plus `MassTraceDetection.cpp:675` and `ElutionPeakDetection.cpp:463/548` — the last
two **inside** the global critical — become real moves. ~6-10 deep copies of every trace per
window eliminated, and the lock hold time in OI-1's funnel becomes constant instead of
occasionally O(n).

It also repairs three things silently:
- `:706` `mts[i] = MassTrace()` currently **frees nothing** (measured over 20,000 traces:
  `new=0 free=0`) because the rvalue binds to `const MassTrace&` and `vector::operator=` retains
  capacity. With the patch: `new=0 free=40000`. The documented memory optimisation starts working
  *as written*, with trace order preserved — which is why this is the right fix and reverse-iterate
  + `pop_back` is not.
- `:699` `out.reserve(n_in)` is a **guaranteed** under-reserve (the adjacent comment says "split
  only ever grows the count"), so `out` reallocates and deep-copies millions of traces.
- `:708`'s `part` is unreserved and filled *inside* the program-global critical.

**Why it cannot change output.** A defaulted move and a copy produce identical values.
`OpenMS::String` is nothrow-move-constructible so the defaulted move ctor is `noexcept`, meaning
`vector` growth uses it. No numeric path changes. **Bit-identical by construction** — this is the
only fix on the list with no ordering exposure at all.

**Falsifying measurement.** Rebuild libOpenMS; any change in dataset D peptide count at all (not just
beyond noise) falsifies the "identical values" claim. Memory: per-window peak RSS should drop
visibly from ~5.8 GB.

**Blocker.** Requires rebuilding `libOpenMS`, not just the TOPP tool, so it is out of scope for
the single-file build path used here.

### OI-3 — Gate the `[merged-trace]` instrumentation. `src/spextract.cpp:721-737`

**Change.** The block is guarded by `if (!mts.empty())` — **only the string assembly at `:734` is
guarded by `span_log`**. So a full scan of every peak of every one of ~8 M traces, a 64 MB
`vector<double>`, and an 8 M-element `sort` run **unconditionally, per window, serial**. The result
is printed only when `win_done == 0` (`:1559`). **23 of 24 windows compute it and throw it away.**
Add a `bool want_spans` parameter and gate the whole block.

**Expected gain.** ~0.6 s/window of pure diagnostics on the critical path, plus 64 MB of dead
allocation and ~8 M page touches per window. Small in wall terms, free in risk terms.

**Why it cannot change output.** `spans` is never read by anything but the log string. Note a
latent bug it also fixes: `:714` appends the valley-split summary to `*span_log` and `:734`
**assigns** over it, so the valley-split line is currently discarded whenever splitting ran.
Gating restores it — a *log* change, not an mzML change.

**Falsifying measurement.** Any diff in the emitted mzML. There should be none.

### OI-4 — SoA the RT gate. `src/spextract.cpp:1074-1075`, build at `:1574`

**Change.** `:1074-1075` binds `const Trace& f` and reads **only `f.rt`** — 8 bytes at offset 8 of
a 56-byte object — to reject ~99.7% of the IM band (`gate:delta_rt` 3.0 at `:488`, against a
1,857 s gradient). 24 of those 56 bytes are the `xic` vector header, permanently
`{nullptr,nullptr,nullptr}` after the explicit swap at `:1585`. Build a `vector<double> frag_rt`
alongside the existing `frag_im` at `:1574-1575` and gate on that.

**Expected gain.** 7.0× less traffic on the gate scan: ~15.7 MB → ~2.24 MB per precursor
(`N_band ≈ 8 M × 0.035 ≈ 280 k`, using the project's own measured ±0.01 1/K0 slab occupancy at
`docs/REVIEW-2026-07-22.md:672-674`). This is ~3-5% of the loop *and* it relieves a shared
resource: 24 concurrent streams at 7 GB/s each is 168 GB/s of demand on a node that sustains
~120 GB/s, so the per-thread rate is below 7 GB/s today and cutting traffic helps every thread.

**Why it cannot change output.** Same `double`, same `fabs`, same comparison, same iteration
order — a pure reordering of the same memory. **Placement is load-bearing**: `frag_rt` must be
built after the `remove_if` at `:1563` *and* after the sort at `:1568`, because `fi` indexes the
post-sort array. That is exactly where `frag_im` is already built, so the pattern is proven
in-tree.

**Cost the fix owes.** +64 MB/window × 24 = **1.5 GB RSS** on a memory-bounded run. The strictly
better version packs the four scalars into a 32-byte SoA and drops `frag_traces` entirely after
`:1585`, which *saves* 192 MB/window.

**Falsifying measurement.** mzML diff (should be empty); loop wall should move by ≤ 5%, so if it
moves by 30% the model of where time goes is wrong.

### OI-5 — Release `g.rt`'s capacity. `src/spextract.cpp:440-442`

**Change.** `g.rt` accumulates **every** XIC RT of every fragment — ~56 M doubles, unreserved
(`total` is computed at `:447-448`, one line *after* the loop that needed it) — then `sort` +
`unique` + `erase` reduces it to ≤ ~1,343. `erase` does not release capacity. Measured at dataset D
shape: `size=1,341  capacity=67,108,864` = **537 MB still held**, × 24 windows in flight, for the
whole rest of the window body. One line: `vector<double>(g.rt).swap(g.rt);` after the erase.

**Expected gain.** ~12.9 GB of resident memory returned. **No CPU win** — the 56 M sort is only
~0.55 s/window, 0.35% of the loop, so "avoid the sort" is not a speedup. This buys headroom
against the admission controller and against OOM on a shared node, nothing more.

**Why it cannot change output.** `G = fg.rt.size()` (`:443`) enters `g.mean`/`g.invnorm` at
`:461-462` and the Pearson at `:1066`/`:1115`. Shrinking *capacity* does not touch *size* or
contents, so `G` is untouched. **The bitmap variant that would avoid the sort is blocked**:
`map.clear(true)` at `:650`/`:674` empties `wmap` inside `detectTraces_`, before `buildFragGrid`
runs at `:1579`, so the frame RTs would have to be snapshotted before `:1555`.

**Falsifying measurement.** Peak RSS does not fall ⇒ the capacity was not the term.

### OI-6 — `std::move(tmp_spec)`. `OpenMS/.../FEATUREFINDER/MassTraceDetection.cpp:258`

`:256-258` builds `PeakMap::SpectrumType tmp_spec(it)` (full copy including float data arrays),
`select`s it down, then `work_exp.addSpectrum(tmp_spec)` — an **lvalue**, so it copies again. One
word removes the second copy. Runs per band, i.e. 5× per window. Output-identical trivially.
Requires a libOpenMS rebuild.

### OI-7 — Reserve `frags` / `frag_scores`. `src/spextract.cpp:1166-1167`

Genuinely unreserved and they grow past the 500 cap on 46.7-75.7% of spectra
(`docs/merge-design.md:82-84`). `touched` at `:1043-1044` is *already* reserved — do not "fix" it.
Note hoisting `frags` out of the loop is partially defeated by `frags.swap(keep)` at `:1141`,
which hands the hoisted buffer the small capacity; fix `assembleFromList_` to write into a
caller-owned buffer instead of swapping. Minor.

---

## 4. BEHAVIOUR-CHANGING FIXES — these need the dataset D quality gate

**None of these may ship without a full dataset D benchmark** (`harness/bench.py`, Sage screen then
MSFragger confirm, per `docs/search-engine-policy.md`). Decision rule: progress = peptides
≥ 10,072 and PSMs/peptide not worse; regression = peptides fall ≥ 1%.

### BC-1 — Band MS1 tracing. `src/spextract.cpp:1344-1346`

**Change.** `const int ms1_bands = std::max(1, getIntOption_("perf:trace_bands"));` — and
`perf:trace_bands` defaults to **0** (`:545`), so this is `max(1,0) = 1`. **MS1 tracing is
unbanded under every shipped configuration.** The auto-expansion `n_bands = ceil(threads/windows)`
lives at `:1489-1492`, **145 lines after** this call and local to the window loop.
`MassTraceDetection.cpp` contains **zero** OpenMP directives (verified: every "omp" hit is the word
"computational"/"computed" in doc comments), so MS1 MTD is genuinely one thread on one core.

The comment at `:1337-1343` asserting the fix was "passing n_bands", and citing `:1537` for the
MS2 call (the actual call is `:1555`), **is false**. That comment should be corrected regardless
of whether the fix ships.

**Expected gain.** 330 s = 7.9% of the run, almost entirely serial. Even a 4× here is ~250 s, and
by Amdahl this phase caps *any* total speedup at 12.7×.

**Why it is behaviour-changing.** Banding partitions peaks by m/z with a halo; a truncated MS1
trace changes the neutral mass fed to `inferPrecursors_` at `:1355`, which changes charge and
monoisotope assignment. This is the *riskiest* place in the tool to perturb trace boundaries.

**Falsifier / gate.** dataset D peptides at 10,072 ± noise with `-perf:trace_bands N` for N ∈ {4, 8}.
Note this option **cannot be tested in isolation** (see BC-3).

### BC-2 — Fix the halo formula. `src/spextract.cpp:639`

```cpp
const double h = edge[b] * ppm * 1e-6 * 20.0;
const double lo = edge[b] - h, hi = edge[b + 1] + h;
```
`h` is derived from the **lower** edge and applied to the **upper** boundary, and `:628` forces
`edge.front() = 0.0`, so **band 0's upper halo is exactly zero**. Traces whose centroid sits just
below `edge[1]` are truncated, and the core rule at `:664` then keeps the truncated band-0 copy and
discards band 1's complete one. The help text at `:545` and the comment at `:613-615` ("banding
with a halo is exact") are wrong.

Scope correction: the truncation only bites where `20 × edge[b] < edge[b+1]`, which over a fragment
m/z range holds **only for b = 0**. So this is one ~0.005 Da sliver at `edge[1]`, and it does *not*
grow with band count.

**Behaviour-changing** — it changes the trace set at that boundary. It is a correctness fix, so it
should be gated but is expected to be neutral-to-positive.

### BC-3 — Raise `perf:trace_bands`. `src/spextract.cpp:545`, `:1344`, `:1489`

**Not a single-axis experiment.** The same option feeds `:1344`, so `-perf:trace_bands 30`
simultaneously flips MS1 tracing from 1 band to 30. Any measurement of it is BC-1 and BC-3 fired
together. Split the option first.

Further, `:652` is `num_threads(bands)` over `bands` iterations — **inner threads always equal
inner iterations**, so raising the band count reproduces the same zero-headroom static mapping at
a larger size. The defensible version is to decouple them: keep `n_conc = 24` and give `:652` a
band count *greater* than its thread count so `schedule(dynamic)` can actually balance.

Two further costs: halo double-consumption (a halo peak can be consumed by a kept trace in band
*b* **and** a different kept trace in band *b+1*, because MTD's visited-set is per-`run()`) scales
with the number of boundaries — roughly 0.08% of peaks at b=5 and **~0.6% at b=30**, above the
0.06% noise floor. And `MassTraceDetection.cpp:264-268` **throws** `InvalidValue` on a map with
< 3 spectra; `:626-627` can emit duplicate quantile edges on tied m/z, producing an empty band.
The throw is captured at `:1529` and rethrown, killing the run. Probability rises with band count.

### BC-4 — Rebalance the window loop. `src/spextract.cpp:1497`, `:1526`

24 iterations on 24 threads with `schedule(dynamic)` is a static one-to-one map with zero
balancing headroom; `n_bands` is frozen at 5 before the loop; nothing rebalances. Window wall time
varies severalfold (`docs/REVIEW-2026-07-22.md:124` puts the density peak at m/z 468-790, windows
2-15; `harness/collate.py:171-175` shows window 0 is 141.3 m/z wide against 32.8 for the next). The
loop's wall is the slowest window and its tail runs on **5 threads of 120**. Any fix here (adaptive
band counts, splitting the widest tiles) changes trace boundaries ⇒ gate it.

Note: the frequently-quoted "window 0 holds 34% of the work" (`s = 0.34`) is **unsourced**. Nothing
in the repo supports it, and every ceiling derived from it — the "1.8× hard ceiling", "8.2 effective
cores", "87 cores of kernel time in the first 322 s" — is downstream of that one number plus a
parameter `f` that was *solved* to reproduce the measurement it is then said to predict. Do not
plan against those figures.

### BC-5 — Add an index tie-break to the canonical trace sort. `src/spextract.cpp:1568-1573`

The comparator is `(im, mz, rt, intensity)` with no final tie-break, and `std::sort` is unstable.
It is therefore **not a total order**, and the run is not bit-reproducible today. Adding a stable
tie-break would make the baseline reproducible — which is behaviour-changing exactly once, and
worth doing before any further optimisation work, because it converts every OI fix's caveat into a
provable "no diff". Same argument applies to the canonical spectrum sort at `:1748-1764`, which
compares RT, precursor m/z, charge, spectrum size and every fragment but **not** `pc.im`, though
`pc.im` is written to the output at `:1153`.

---

## 5. REJECTED — attractive, but not supported by the code or the profile

| proposal | why it is rejected |
|---|---|
| **`OMP_MAX_ACTIVE_LEVELS=1` as the zero-code A/B** | **No-op.** `:1502` calls `omp_set_max_active_levels(2)` at runtime, overriding the environment ICV. Measured: team sizes and total live threads unchanged. Had it worked it would also have serialized the *deliberate* banded-MTD region at `:652`, confounding the result. Use `OMP_NUM_THREADS=120,1`. |
| **"The named critical explains the 43.5% kernel time"** | Off by ~40×. Measured lock throughput ~16.6 M emits/s aggregate ⇒ 192 M emits ≈ 170-550 CPU-s ≈ 0.3-0.9% of 64,521 s. The lock is the funnel; the oversubscription is the cost. |
| **Memory fixes will "raise `n_conc`"** | Twice refuted. `project()` at `:1518-1522` returns `compact_bytes*10 + 256 MB` — a pure function of *input* size, blind to actual allocations, so freeing memory cannot move it. And `n_conc = min(24, 120/5) = 24 = window_list.size()`: **there is no concurrency left to buy.** |
| **Intensity prefilter before scoring (the `key = c·I`, `c ≤ 1` theorem)** | The theorem is sound (`:1170-1171`, `:1135-1137` is a strict total order so the top-k set is unique). It is *pointless*: `scoreCandidates_` is 2-4% of the window wall. Also three implementation traps: `S` is undefined until 500 *emitted* survivors exist (`:1081`, `:1116` can cull the top-500-by-intensity); `g.flat` stores intensities as **float** (`:457`) while `mean`/`invnorm` are computed from **doubles** (`:458-462`), so `c` can exceed 1 by O(1e-7) systematically; and `toTrace` never dedupes `xic` (`:181`), so duplicate grid indices would make `invnorm` too large and `c > 1` by an unbounded amount. |
| **Cache-tiling / blocking / sparsifying `pdense`** | `pdense` is `vector<float>(fg.rt.size())` at `:1586` ≈ 5,372 B = 84 cache lines. L1-resident six times over. Nothing to tile. |
| **Sorting or reordering the CSR traversal** | Already strictly ascending. `frag_traces` is IM-major sorted (`:1568`), `fg` is built in that index order (`:450`) so `g.off` is monotone in `fi`, `lo`/`hi` come from `lower_bound`/`upper_bound` on sorted `frag_im`, and within a fragment `f.xic` is sorted (`:181`) with `gi` from `lower_bound` on sorted `g.rt` (`:456`). |
| **False sharing on `pdense`** | Declared *inside* the loop body at `:1586` — per-iteration, thread-private. The only genuinely shared adjacent write is `win_out` (`:1472`): 24 × 24 B = 9 cache lines, ~31 k header writes over 3,805 s. |
| **`nth_element` instead of `partial_sort` at `:1134`** | ~18 k comparisons per emitted spectrum against a ~400 k × 56 B scan for the same precursor — four orders of magnitude too small to list as a win. And not strictly safe: `MSSpectrum::sortByPosition` uses `std::stable_sort` when no data arrays are present, and `out` at `:1146` has none, so push order *is* preserved into the emitted spectrum. |
| **Drop fragments with < `min_corr_pts` support once per window** | **Behaviour-changing, mislabelled.** `g.rt` at `:440-442` is the union over **all** fragments including those. Dropping them removes grid points ⇒ `G` changes ⇒ `pmean = psum/G` at `:1066` changes ⇒ **every** Pearson at `:1115` changes and the precursor resampling at `:1047-1055` snaps to a different grid. |
| **Rebuild `g.rt` from `wmap`'s frame RTs** | `wmap` is **destroyed** before it could be read: `map.clear(true)` at `:650` (banded) and `:674` (unbanded), inside `detectTraces_`, while `buildFragGrid` runs at `:1579`. Salvageable only by snapshotting the RTs before `:1555`. |
| **Reverse-iterate + `pop_back` to make `:706` free memory** | Changes trace order for no reason. OI-2 makes the existing line work as written, order preserved. |
| **IM-sort the precursor visit order** | After OI-4 the band is 2.24 MB/precursor ≈ 0.32 ms × 40 k = 12.8 s/window = **0.34%** — below the project's own 0.06%/0.2% noise floor. And it costs a permanent `apportion == 0 && rp_max <= 0` gate that silently changes behaviour if those defaults ever move (`:1618` is order-dependent float accumulation; `:1672-1679` tie-breaks on `(long)pi`). |
| **Hoist `G * pmean` at `:1115`** | `*` is left-associative, so it already evaluates as `(G * pmean) * fg.mean[fi]`. Hoisting is plain CSE, bit-identical — and worthless, GCC already does it. |
| **Parallelise the per-precursor scoring loop** | Correct in principle (indexed output preserves `pi`-ascending order ⇒ output-identical), but it targets 2-4% of the window wall. Do it after the phase timers exist, if at all. |
| **THP / `MADV_HUGEPAGE`** | A real observation (THP is `madvise`, the tool never calls `madvise`, ~150 GB on 4 KB pages) but page faults account for ~24 s of the 28,065 s system time. Retest *after* OI-1; if kernel time is still high it becomes the next hypothesis. |
| **`TCMALLOC_RELEASE_RATE=0`** | tcmalloc was already `LD_PRELOAD`ed in the measured run. The allocator lever is spent; the remaining allocator cost is addressed by OI-2, not by tuning. |

---

## 6. Order of work

1. `grep Threads /proc/<pid>/status` mid-window-loop. One command.
2. **OI-1** (implemented) — measure the window-loop wall and the kernel share.
3. The six `phase_clock_()` calls of §1.6. Everything after this point should be planned against
   measurements, not against §1.4.
4. **OI-2** + **OI-6** (libOpenMS rebuild), then **OI-3**, **OI-5**, **OI-4**.
5. Only then consider §4, each behind a full dataset D gate, one axis at a time — and split
   `perf:trace_bands` into MS1 and MS2 options first (BC-3).
