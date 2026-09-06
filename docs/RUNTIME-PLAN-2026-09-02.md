# Runtime + memory plan (2026-09-02) — measured, not guessed

## The profile that replaces every hypothesis
Per-phase timers added 2026-09-02 (`[perf]` table; CPU from getrusage so the parallel factor is
real). dataset D, release binary, 100 threads, 128-core/995 GB node, warm page cache:

| phase | wall s | share | parallel | note |
|---|---|---|---|---|
| LOAD (streaming .d) | ~1,046 | 53% | **1.0x** | serial: zstd + TOF->m/z + scan->IM + MSSpectrum build, on one thread |
| WINDOW_LOOP | 769 | 39% | **63x** | already well parallelised |
| MS1_TRACE | 83 | 4% | 5.0x | |
| PRECURSOR_INFER | 63 | 3% | 1.0x | serial, small |
| SORT + ASSEMBLE + WRITE | 24 | 1.2% | 1.0x | **not** a lever |

Every prior guess about "the unattributed 53%" (canonical sort, consolidation, the 9.3 GB mzML write)
totals 1.2%. The whole runtime problem is one serial phase, and the 27x whole-run parallel factor
is that phase dragging a 63x loop down. NOTE: the clock epoch was set on first call, which on the
streaming path is the "load done" line -- so the load was invisible to all earlier `[t=..]`
markers. Fixed (epoch anchored at main_).

## Root causes in the loader (verified in source, file:line in the archaeology reports)
1. **Every MS2 frame is decoded once per isolation window of its group** (window-outside-frame loop,
   BrukerTimsFile.cpp:1589-1610; `save_to_buffs` does not cache, it `close()`s). dataset D: 2 windows per
   group -> **2.0x re-decode of the entire 6.2 GB run**, from SQL: 16,105 physical frames, 32,210 decodes.
2. **The loop is strictly serial** (`for (fid..) consumer.consumeSpectrum(spec)`), and our consumer
   peak-picks + compacts 1.26e9 peaks inline on that thread.
3. I/O is NOT the cost: tdf_bin is 6.18 GB, contiguous append-only (0 out-of-order TimsId pairs of
   17,448), mmap'ed -- ~3 s sequential on NVMe. It is CPU work: zstd, per-peak conversions, allocation.
4. opentims documents the thread-safe decode recipe (per-thread ZSTD_DCtx + buffer: decompress(buf,ctx)
   -> save_to_buffs(..,ctx) -> close()) and uses it itself (opentims.cpp:599-634). No lock, no upstream change.

## Measures, ranked by (win) / (risk x effort). Status as of writing.
| # | measure | expected | byte-identical? | status |
|---|---|---|---|---|
| 1 | **Parallel batched pick+compact in our consumer** (64-frame batches, ordered append) | load 1,046 -> 637 s; total 33:05 -> **25:28 (-23%)** | intensities differ at 1e-5..1e-4 rel. (chaotic amplification of last-ulp), m/z identical, **Sage 12,200 = identical peptide count**; self-determinism run pending | **DONE, measured** |
| 2 | **Decode-once, frame-major, parallel decode in the loader** (per-thread ZSTD_DCtx; batch via SPEXTRACTOR_LOAD_BATCH) | removes the 2x re-decode AND parallelises the remaining ~637 s serial path; predicted load -> **~50-100 s** | expected identical to #1's output (per-window frame order preserved) | **built; sweep 64/256/1024 running** |
| 3 | Cull traces with < min_correlation_points support before EPD/grid (codex #2; "R2") | 8-18% wall, 15-35 GB RSS | yes | next |
| 4 | Link tcmalloc/jemalloc in CMake (kimi: July LD_PRELOAD A/B measured -18.4% RSS, -10.7% wall); free A/B `MALLOC_ARENA_MAX=4` first | -20..-60 GB RSS | yes | queued (needs quiet node) |
| 5 | Band construction: two lower_bounds instead of O(bands x peaks) scan ("R3") | 3-8% wall | yes | next |
| 6 | Fixed-bands concurrency cap sweep (n_conc 8/4/2 with trace_bands=12 held) | 25-40 GB RSS, 0-15% wall cost | yes | queued; the old "more windows is slower" sweep was CONFOUNDED (it also changed bands 12->4) |
| 7 | MS1 loop (1,343 frames) same treatment as #2 | ~small (4% of frames) | yes | after #2 |
| 8 | Cross-window task graph ("R5") | 15-30% of what remains | no | only after 1-6; its 2x claim ignored load |
| - | Parallel/compressed mzML write | write is 22 s | - | **dropped** (kimi's #1, refuted by measurement) |
| - | Naive dynamic thread apportionment, spectrum merging/collapse | - | - | falsified earlier; do not re-propose |

Ring buffer / producer thread (the question asked): with #2 the design is batch-synchronous
(decode a batch in parallel, hand off in order, pick in parallel). Decode-agent sizing formula:
N ~ k*(2B + D*t_p/t_d), k in [1.5,2]; with decode parallel, decode is the bottleneck, so a depth
of 2-3 batches suffices and extra depth is wasted RAM -- spend threads on decoders, not on depth.
Memory per buffered frame ~1.5 MB (37k peaks), so 256 frames ~ 400 MB, negligible vs the 12 GB
compact store. The batch sweep measures this instead of trusting the formula. A separate reader
PROCESS is unnecessary: mmap + per-thread contexts give the same overlap in-process.

## Memory (120 GB peak, measured inside the window loop)
Peak RSS is a temporal high-water mark inside WINDOW_LOOP (RSS 27 GB after load, 87 GB at loop
end, 124 GB peak). July measured ~87% of a similar peak as glibc arena retention (tcmalloc -18.4%).
Live data: 12 GB compact store + per-window working sets + ~10 GB of buffered output spectra.
Order of attack: #4 (free), #3 (cuts working sets), #6 (caps concurrency), then streaming the
output instead of buffering 924k spectra.

## Pre-registered predictions for the decode-once loader sweep (recorded before the numbers)
Kimi's first-principles decomposition (plan review, 2026-09-02): post-#1 load 637 s = ~619 s serial
decode+build + ~18 s pick wall; decode-once removes ~75-135 core-s; what remains (~485-545 core-s
decode/slice/build + ~430 core-s pick) is fully parallel -> at 35-70% efficiency on 100 threads
**LOAD wall 15-30 s (central ~20-25 s), whole run ~16 min.** My own number in the table above was
50-100 s, which -- as kimi points out -- silently assumes ~10-15% parallel efficiency. Reading key:
<12 s suspect the timers; 15-30 s kimi confirmed; 50-100 s = contention (allocator / NUMA first-touch),
not "confirmed"; >150 s or parallel factor ~1.0 = the pragma did not take. Batch 64/256/1024 should be
within +-20% of each other; a large win for 1024 means per-batch overhead nobody has seen in the code.
Sweep is gated to start only after the concurrent MSFragger run has exited (it was not, initially).

## Plan review verdicts (kimi + codex, 2026-09-02) and what changed
**Bugs found and fixed (before the sweep binary is built):**
- **Chunking capped the parallelism** (codex #6/#7). `schedule(dynamic, 8)` over a 64-frame batch =
  8 chunks = 8 workers; the "parallel" pick ran 8-way, so 1046 -> 637 s is the 8-way number, not the
  100-way one (implied serial pick ~467 s, residual ~58 s). The loader's `(dynamic, 4)` over 256
  frames capped the decoder at 64 workers (16 at batch 64). Both are now chunk 1, pick batch 256
  (`SPEXTRACTOR_PICK_BATCH`), and the running sweep will build this binary first.
- **Exceptions in both OMP regions** (codex #11) would have aborted the process; both loops now
  capture and rethrow serially like the full-load path.
- **Vendor-SDK converter thread-safety is unvouched** (codex #9, kimi): the loader runs serially
  when the SDK converter is active unless `SPEXTRACTOR_SDK_PARALLEL` is set. The TDF-table converter
  is read-only after construction. Cluster runs use the table (provenance `tdf_table_modeltype1`).
- **Sweep timing gated** on the concurrent MSFragger job exiting (it was not).

**Acceptance criterion, replacing "Sage count identical":** (a) two fixed-config runs digest-identical
(`bench/semantic_digest.py`, pickdet running); (b) baseline vs optimized: same spectrum count, keys,
RT, precursors, quantized m/z/IM, and a pre-registered intensity envelope with the cause explained;
(c) both engines: accepted-peptide SET overlap, not count. **Measured for (c) now:** the v0.2.0
serial-pick dataset D run (`final_t100`, 12,200) and the batched-pick run (12,200) accept **exactly the same
12,200 peptides (symmetric difference 0)**. For calibration, the same-code run at another thread count
(`tbl_s30`, 11,976) differs from it by 2,730 peptides (20% of the union) -- the cross-config chaos is
a 20% membership churn at a 2% count level, so count equality alone is indeed a weak test and set
equality is a strong one. The 1e-5..1e-4 intensity differences are tens-to-hundreds of float ULPs
(codex #2), cause not yet shown; `SPEXTRACTOR_PICK_SERIAL=1` on the new build is the direct A/B.

**Pre-registered load predictions for the sweep:** kimi 15-30 s (fully parallel, 35-70% efficiency);
codex 105 / 75 / 70 s for batch 64 / 256 / 1024 -- with the OLD chunking; with chunk 1 codex's
worker cap no longer applies, so codex's numbers become an upper bound. Mine (50-100 s) is withdrawn.

**Re-ranked measures (after the sweep lands; wall and RSS ranked separately, codex #15/kimi):**
1. determinism gate (#1, #2) -> 2. allocator env A/B `MALLOC_ARENA_MAX=4` (free, byte-identical to
verify) -> 3. exact `lower_bound` band endpoints (#5; codex #13: same half-open interval, provable)
-> 4. WINDOW_LOOP parallel efficiency 63x -> ?, PRECURSOR_INFER (63 s serial) and MS1_TRACE (83 s
at 5x) -- three rows the plan lacked (kimi) -> 5. short-trace skip ONLY at the dot-product-only
location (#3 demoted: pre-EPD/grid cull is behaviour-changing, codex #12) -> 6. tcmalloc/jemalloc
as separate A/Bs -> 7. concurrency cap for RSS (#6; no memory target stated yet: 124 GB on a 995 GB
node) -> 8. output streaming, cross-window task graph (#8) re-estimated against the post-sweep
profile (codex #20: reprofile decode/convert/slice/sort/pick/compact before assigning percentages).
**Loose ends kept:** 33:05 - 409 s = 26:16, not 25:28 -- 48 s moved elsewhere, per-phase tables to be
compared (both reviewers); the memory section's 400 MB/batch undercounts (codex #18): the sweep's
RSS column is the measurement; the "ring buffer" prose over-sold a pipeline that is batch-synchronous
(codex #19) -- with decode and pick both parallel there is no producer/consumer overlap to size, only
the per-batch barrier, which chunk 1 minimises.

## Self-determinism: FAILED, and what the failure looks like (2026-09-02)
Two same-config runs of the batched-pick binary (100 threads, dataset D) are not digest-identical. Per-spectrum
diff of the first 85,144 spectra: **2,001+ spectra differ; in each, EVERY peak's intensity differs by a
common ~1e-4 relative factor, m/z untouched** (e.g. 18.34720 vs 18.34847). That is not a picker race
(picking is per-spectrum pure: `compactify` and `ensureIMArrayName` hold no shared state, MassTraceDetection
has no OpenMP inside, band edges are an exact data-derived partition, band results are concatenated by
index) -- it is a per-precursor quantity moving by 1e-4 and propagating into every fragment's
`corr^corr_power` weight. The accepted-peptide SET is nonetheless identical (12,200 = 12,200, symdiff 0),
so it is quality-neutral, but the README's fixed-thread reproducibility claim is false as stated and the
serial-pick v0.2.0 streaming binary was never self-tested. Running now (`bench_det2.sh`): S1/S2 (serial
pick) and P1/P2 (parallel pick) on the fixed build -- if S1==S2 the batched pick is the entry point; if
S1!=S2 it is downstream. Chained (`bench_det3.sh`): `SPEXTRACTOR_DET=1` logs an order-insensitive digest of
MS1 traces, precursors, and each window's fragment traces, so two runs pinpoint the first stage that differs.

**The 48 s "discrepancy" (both reviewers):** per-phase tables of the three runs resolve it as run-to-run
variance on a shared node, not arithmetic -- MS1_TRACE 83 / 59 / 81 s, PRECURSOR_INFER 63 / 32 / 56 s,
WINDOW_LOOP 769 / 770 / 848 s across instr_s30 / parload_s30 / parload_s30b (the last ran beside an
MSFragger job). Serial phases move by +-30 s between runs; quote LOAD deltas, never totals, and never
time a sweep beside another job (now enforced in the scripts). The batched-pick LOAD parallel factor was
**3.2x** (CPU 2,028 s / wall 637 s), which is the 8-worker chunking cap made visible.
**det2 phase 1, first half (10:35):** S1 == S2 (serial pick, chunk-1 loader): digest
`17850198e0ce…` both runs, so the serial-pick path is self-deterministic on this build and the
nondeterminism enters at or after the parallel pick. LOAD with the serial pick is 801 s at 1.6x --
the pick, not the decode, is the serial floor; P1's LOAD wall is the number the sweep was for.
**Nondeterminism, code-reading hypothesis (10:55):** `PeakPickerIM::pickIMCluster` is `const` and
per-spectrum pure, but seeds clusters from `std::sort(intensity desc)` -- an UNSTABLE sort on raw TDF
counts, where ties are everywhere. A tie broken the other way moves one boundary point between two
clusters: intensity moves by that point's count (1 on 1e4 = 1e-4), m/z by ~1e-6 Da (invisible after
1e-5 Da quantisation) -- exactly the observed signature. Same input must give the same tie order in the
same binary, so the open question is why the parallel path's INPUT differs: per-thread history
(firstprivate picker state) or a race/uninitialised read. det4 (chained behind det3, same build):
`schedule(static)` pair T1/T2 -- T1==T2!=dynamic -> history; T1!=T2 -> race; T1==dynamic -> neither.
**det2 phase 1 complete (11:36): S1 == S2 == P1 == P2**, one digest `17850198e0ce…` for all four
runs (serial pick x2, parallel pick x2, chunk-1 loader). The fixed build is byte-identical at fixed
config AND the parallel pick is byte-identical to the serial pick. The earlier "1e-5..1e-4
intensity drift in ~2% of spectra" was the chunking-bug build (8-way pick / 64-capped loader), not
a property of the parallel pick. Consequence: the stage-digest localisation runs (det3/det4) are
unnecessary and were cancelled; the diagnostic code was stashed, not merged.
Measured: LOAD 801 s (serial pick) -> 432 s (parallel pick, 4.1x); total 28:06 -> 22:01.
Remaining sweep: batch 64, batch 1024, MALLOC_ARENA_MAX=4.

## mzPeak reconversion (mzpeak-convert v0.9.2, ims-compact integer-TOF) -- 2026-09-02
| file | vendor `tdf_bin` | old .mzpeak peaks table (Jul) | **v0.9.2 archive peaks table** | container incl. vendor side-files |
|---|---|---|---|---|
| dataset D | 6.18 GB | 5.68 GB | **6.46 GB** (zstd 3) · 6.44 (zstd 9) · 6.39 (zstd 15) | 9.98 GB (+3.52 GB vendor) |
| dataset A | 8.95 GB | 9.45 GB | **9.45 GB** | 13.69 GB (+4.24 GB vendor) |
- Layout: `point` table {spectrum_index u64 delta, intensity i32 BYTE_STREAM_SPLIT, 1/K0 f64 dictionary,
  tof i32 BYTE_STREAM_SPLIT}, all ZSTD; 2,471,762,904 rows for dataset D = exactly the .d peak count;
  17,448 spectra (1,343 MS1 / 16,105 MS2). Conversion 13 s (dataset D) / 21 s (dataset A) on the M-series Mac.
- The archive is 4.5-5.6% larger than the vendor binary; zstd-15 buys 1% for 2.3x the time;
  `--no-tims-recalibration` changes nothing (6.45 GB). The 1/K0 column (ratio 0.57) and tof (0.52)
  are the bulk; the old dataset D file being smaller (5.68) is not reproduced by any v0.9.2 setting.
- **Converter bug (v0.9.2): `--ims-chunked` aborts on TDF** -- "layout family mismatch between spectrum
  facets: spectra_data is 'point' but spectra_peaks would be 'chunk'". Both dataset D and dataset A.
- The archive still stores raw TOF + the two-point transform (`mzpeak:transform_params`); a reader
  must convert, and the exact ModelType-1 calibration is NOT in the archive (only in the embedded
  `vendor/analysis.tdf.gz`). m/z accuracy from mzPeak therefore = the two-point model unless the
  reader re-derives the table model from the embedded tdf.
- Reader: OpenMS' in-tree `MzPeakFile` (2,571 lines, direct Arrow) has no ims-compact/tof or
  per-peak-IM support and returned 0 spectra on these files; the mzPeak C++ library fork
  (okohlbacher/mzpeak-openms trunk 139cc64, 2026-09-01; 164 ahead of / 1 behind OpenMS/mzpeak)
  reads ims-compact, per-peak IM and per-window mobility bands. SpeXtractor now carries
  `src/MzPeakStreamLoad.h`: a streaming .mzpeak loader on that library with the same (MS1 first,
  then per-(frame, window)) hand-off contract as the .d loader; build-optional via `-DMZPEAK_ROOT`.
  Cluster build: library against contrib Arrow 23 (pin relaxed from >=24), Boost 1.89 from the
  OpenDIAlyzer env, libzip 1.11.4 from contrib. Benchmark queued: dataset D archive vs .d, same flags.
**Measure #4 (allocator), arena arm, 13:02:** `MALLOC_ARENA_MAX=4` at 100 threads is HARMFUL: the window
loop ran at ~636% CPU instead of ~6300% (allocator serialised across threads), 16/24 windows after 56 min
vs 22 min for the whole run; aborted. The July -10.7% was not this setting at this thread count. The
remaining allocator arm is a real replacement allocator (tcmalloc/jemalloc via LD_PRELOAD), not arena capping.
**mzPeak reader wiring (13:20):** library built on the cluster (fork trunk, GCC 13 + Arrow 23 with two
compat patches in `patches/`, recipe `scripts/build_mzpeak_lib.sh`). Reading the v0.9.2 dataset D archive:
17,448 spectra, 32,210 MS2 windows (2 per frame), metadata sweep 0.43 s, 24-frame decode 0.07 s. Two
archive/reader defects worked around in `MzPeakStreamLoad.h`: the per-window 1/K0 band is stored as
name-only CV params (typed limits empty), and `precursor_index` is NULL so the reader attaches both ions
of a frame to the first precursor -> ions are flattened and matched to windows by selected_ion_mz.
Without this every window would have carried the whole frame (2x peaks, garbage IDs).

## LOAD split measured (P256det, 13:26) -- the serial floor is the MS1 path, not the MS2 loader
| inside LOAD(stream) = 444.7 s | wall |
|---|---|
| MS2 batched decode (72 parallel regions) | **4.4 s** |
| MS2 serial hand-off incl. the consumer's parallel pick flush | 11.4 s (pick 11.3 s) |
| **everything else = the MS1 path** (1,343 frames, 927 M peaks: serial decode, serial FrameCentroider, serial `picker_.pickIMCluster` in `consumeMS1Spectrum_`) | **~425 s** |
Measure #7 ("MS1 loop, ~small: 4% of frames") was mis-ranked: MS1 is 4% of frames but 37% of peaks and
100% serial, i.e. ~95% of LOAD wall and ~32% of the whole run. Format-independent -> mzPeak cannot fix
it. Next lever: batch MS1 frames, decode + centroid in parallel (per-thread FrameCentroider), pick MS1 in
the consumer's parallel flush (buffer MS1 like MS2, keep frame order). Expected LOAD 445 -> ~30 s, total
22:22 -> ~15:30.

## mzPeak input vs .d input (dataset D, same binary, same flags, 100 threads) -- 13:55
| | `.d` (P256det) | `.mzpeak` v0.9.2 archive (MZP) |
|---|---|---|
| LOAD(stream) wall | 444.7 s | **395.8 s** (-11%) |
| LOAD CPU / parallel factor | 1,785 s / 4.0x | 7,296 s / 18.4x (system time 1 h 05 min vs 10 min) |
| WINDOW_LOOP | 774.3 s | 765.4 s |
| total | 22:22 | **21:24** (-4%) |
| pseudo-MS2 written | 922,902 | 925,595 |
| peak RSS | 107 GB | 112.5 GB |
- Verdict: **reading from mzPeak does not resolve the performance issue.** The serial floor inside LOAD
  is the MS1 pick in the consumer (~380 s), which is format-independent; both paths pay it. mzPeak
  removes only the small MS2/MS1 decode share and spends 4x the CPU doing so (a reader Index was opened
  per 12-frame work item; now thread-local, expected to cut the system time).
- Spectrum content differs from the .d path by construction (two-point m/z transform vs exact table
  model; converter-recalibrated 1/K0 vs our rational TIMS calibration): 50/50 stage digests differ,
  +0.3% spectra. Sage accepted-peptide set overlap and the ppm shift are the meaningful comparison
  (pending in this chain).
- The lever is the MS1 path, for BOTH inputs: parallel MS1 decode in the loader (`[SpeXtractor ms1-par]`)
  and parallel MS1 pick in the consumer (`flushMS1_`), queued as the next timed run.
**mzPeak input, identification (13:58):** Sage @1%: `.d` 12,217 peptides / 97,793 PSMs; mzPeak 10,785 / 89,012
(**-11.7% peptides**; 9,987 shared, 2,230 only-.d, 798 only-mzPeak). Mass errors of accepted PSMs:
`.d` precursor median +4.63 ppm, fragment +3.39; mzPeak precursor +8.60, **fragment +10.05** [p10 +7.96,
p90 +12.36]. The archive's two-point TOF transform sits ~7 ppm off the exact ModelType-1 model and pushes
fragments toward the 20 ppm search tolerance -- the same bow measured on 2026-09-01. mzPeak input is
therefore slower per CPU-second, not faster per wall-second, and less accurate, until a reader applies the
exact calibration (only available from the embedded `vendor/analysis.tdf.gz`).

**B2 CLOSED (14:46): loader batch size does not change content.** specdiff d2_P1 (batch 256) vs d3_B64det
(batch 64): 922,902 spectra, 0 differ; Sage set 12,217 = 12,217, symdiff 0; all 50 `[det]` stage digests
identical. The "SPECTRUM DATA DIFFERS" came from `bench/semantic_digest.py` hashing to EOF, i.e. including
the trailing `<indexList>` byte offsets, which shift when the header changes length (the OpenMS version stamp
carries the build date: 4 bytes longer -> every offset +4). The digest now stops at `</spectrumList>`; with it
the two files hash identically (`ca1ddbbd…`). The 09-01 "differs across thread counts" (B3) was taken with
the old digest and must be re-measured. **dataset A LOAD-only (14:38):** 418 s at 8.0x, MS2 decode 10.1 s,
hand-off 25.0 s -> the ~380 s residual is the MS1 path on this file too.

## A1 LANDED (15:01): the MS1 path parallelised -- LOAD 445 s -> 39 s, total 22:22 -> 15:21
| dataset D, 100 threads | before (P256det) | **after (MS1par)** |
|---|---|---|
| LOAD(stream) wall / CPU / parallel | 444.7 s / 1,785 s / 4.0x | **38.8 s / 3,122 s / 80.5x** |
| of which MS1 decode+sort (parallel region) / hand-off / pick | serial ~429 s | 6.9 s / 15.4 s / 27.2 s (pick incl. MS2 ~11 s) |
| MS2 decode / hand-off | 4.4 / 11.4 s | 4.4 / 11.8 s |
| WINDOW_LOOP | 774 s | 759 s |
| **total wall** | **22:22** | **15:21 (-31%)** |
| peak RSS | 107.0 GB | 105.6 GB |
| spectra / digest | 922,902 / `ca1ddbbd…` | 922,902 / **`ca1ddbbd…` (byte-identical spectrum data)** |
| Sage @1% | 12,217 | **12,217, symdiff 0** |
Pre-registered gate (LOAD <= 45 s, stage digests identical, Sage set identical): **PASSED**. The change is
in the OpenMS loader patch (`[SpeXtractor ms1-par]`: batched parallel MS1 decode with per-thread ZSTD ctx and
FrameCentroider, ordered hand-off) and the consumer (`flushMS1_`: parallel MS1 pick). Remaining wall is
83% window loop; LOAD is now 4%. Baseline for the next measures: 15:21.
**B3 CLOSED (15:14): output is identical across thread counts.** Re-digested with the fixed script: final_t8 ==
final_t100 (`2a1dc5fd…`) and sdkcal_det8 == sdkcal (`0cacb566…`). The 2026-09-01 "spectrum data still differs
between threads 8 and 100" was the `<indexList>` offset artefact. Per-spectrum check on the threads-8 vs threads-100 pair: 923,713 spectra, 0 differ. Claim to ship: byte-identical
spectrum data at any thread count and any loader batch (three independent pairs today).
**mzPeak input, final (15:23, parallel MS1 pick + thread-local reader):** LOAD 104.2 s (80.5x, CPU 8,392 s, system
time still 1 h 04 min -> the reader's decode path, not re-opening, burns it), total **17:04 vs 15:21 from `.d`**,
Sage 10,785 vs 12,217 (unchanged −11.7%). Provenance stamp verified in the output:
`spx:mz_calibration = mzpeak_two_point_transform (...) archive=S30_v092_archive.mzpeak`. Verdict stands: keep `.d`
as the primary input; mzPeak input is 11% slower and 12% less sensitive until the archive carries the exact
calibration.
**Why the mzPeak reader burns 1 h of system time (15:24, from `/usr/bin/time -v`):** voluntary context switches
52.9 M (mzPeak) vs 0.69 M (`.d`) at the same page-fault count (131 M vs 135 M) and zero major faults. That is
futex ping-pong inside the library's decode path (its per-reader mutex and/or the `util/executor` thread pool
behind `get_spectra_batch`), i.e. ~50 M lock hand-offs for 33 k spectra -- not I/O, not allocation. Library-side
fix: a synchronous single-threaded decode mode for callers that parallelise externally, or one shared executor.
**R3 + early short-trace skip (17:41): byte-identical, no measurable wall gain.** Band partition by two
`MZBegin` lookups instead of a per-band scan, and skipping fragments whose support < min_corr_pts before the
dot product: digest identical to P1; WINDOW_LOOP 763 s (vs 759 s), total 15:43, peak RSS 101.8 GB (vs 105.6).
The loop is not bound by either loop; kept as the cheaper code path, not counted as a win.
**A4(i) allocator A/B, first arm (17:57, `LD_PRELOAD=libtcmalloc_minimal.so.4` on the R3 binary):** byte-identical;
total 15:49 (vs 15:43 glibc, noise), **peak RSS 94.3 GB vs 101.8 GB (−7.5 GB, −7%)**, system time 5:12 vs 10:11.
A memory win, not a wall win at 100 threads. Full `libtcmalloc.so.4` arm running.
**A4(i) second arm (18:13, full `libtcmalloc.so.4`):** byte-identical; peak RSS 97.5 GB (glibc 101.8, tcmalloc_minimal
94.3); wall 15:45 (noise). Verdict: tcmalloc_minimal is the better allocator arm (−7% RSS, wall-neutral, half the
system time); adopt via `LD_PRELOAD` in the benchmark harness, link in CMake only when a release build follows.
**A3(iv) `perf:ms1_trace_bands=48` (18:31): WORSE and not output-identical.** MS1_TRACE 144 s at 2.9x (vs 61-80 s at
12 bands, 5.1x) and 922,899 spectra vs 922,902 -- the MS1 band partition is NOT exact (traces near band edges),
so the knob is quality-affecting (BC-3 was right). Keep 12; drop A3(iv).
**D1 first attempt (18:32) aborted:** the cluster copy of the dataset D archive is the --no-vendor one, so
`vendor/analysis.tdf.gz` is absent; the loader now accepts `SPEXTRACTOR_MZPEAK_TDF=<sidecar>`; rerun in flight.

## D1 ISOLATED (18:58): the two-point transform is 91% of the mzPeak loss; exact calibration recovers it
| dataset D input | Sage @1% | vs `.d` (12,217) | paired precursor / fragment ppm vs the `.d` run |
|---|---|---|---|
| `.mzpeak`, archive two-point transform | 10,785 | −11.7% | (vs dt: +8.60 / +10.05) |
| **`.mzpeak`, exact calibration recovered** (tof from the inverted transform + MzCalibration and Frames.T1 from the tdf) | **12,082** | **−1.1%** (symdiff 1,186 / 1,051) | **+0.00 / +0.00** — the m/z axis is identical to the `.d` path |
The remaining −1.1% is the IM axis / band reconstruction (converter-recalibrated 1/K0 vs our rational calibration; ~1.3%
vendor-IM effect measured on 09-01). `SPEXTRACTOR_MZPEAK_EXACT` is now ON by default and fails closed without a tdf
(embedded `vendor/analysis.tdf.gz` or `SPEXTRACTOR_MZPEAK_TDF=` sidecar); `=0` accepts the two-point m/z explicitly.
Provenance: `spx:mz_calibration = tdf_table_modeltype1 (recovered from the archive's two-point transform + ...)`.
