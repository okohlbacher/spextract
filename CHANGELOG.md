# Changelog

## [Unreleased]

### Fixed
- **Runs on the public diaPASEF HeLa benchmark (PXD017703).** Three defects the in-house cohort never
  exposed: the calibration loader demanded exactly one `MzCalibration` row (these files have two;
  it now uses the row the frames reference), rejected a stored `C2 == 0` as missing (a real linear
  calibration on 2020 timsTOF Pro firmware; NULL is now told apart from zero), and -- on a stock
  OpenMS install -- silently fell back to the `openms` detector because `loadTofAxis` sat behind a
  `__has_include` for an in-tree OpenMS header. Both tdf readers now use the sqlite3 C API and a
  missing SQLite is a compile error. First public-data extraction: 136,818 spectra in 57 s.

### Measured
- **`perf:stream_load` is a sensitivity/resource trade, not a free default.** Both readers searched
  with both engines on two datasets plus entrapment: equivalent on the smaller file, but on the
  larger one the one-shot reader finds +223 Sage peptides (+2.0%) at unchanged entrapment FDR
  (1.13% vs 1.15%) -- real signal the streaming path loses, for 1.75x the memory and 1.6x the wall.
  The default stays `true`; `false` is now documented as the more sensitive setting.
- **Entrapment FDR on all six datasets, shipping configuration: 1.13-1.31% at a nominal 1%**, every
  interval overlapping. The peptide counts on record are honest to about a quarter point of FDR.


## v0.3.0 — 2026-09-04

### Benchmark, all six datasets, shipped defaults only

`spextract -in <file>.d -out pseudo.mzML -threads 100` and nothing else:

| file | spectra | wall | peak RSS | Sage @1% | MSFragger @1% |
|---|---|---|---|---|---|
| dataset A | 862,716 | 8:26 | 108.7 GB | 10,909 | 11,935 |
| dataset B | 585,503 | 5:37 | 84.6 GB | 10,272 | 9,691 |
| dataset C | 723,314 | 6:27 | 90.3 GB | 12,149 | 12,516 |
| dataset D | 655,776 | 5:21 | 79.0 GB | 12,482 | 12,337 |
| dataset E | 542,533 | 5:38 | 67.9 GB | 11,217 | 10,585 |
| dataset F | 597,267 | 5:59 | 73.8 GB | 11,362 | 11,049 |

dataset C and dataset F had never been benchmarked before. Peak RSS spans 67.9-108.7 GB across the cohort,
which is where the README's stated 80-125 GB requirement comes from.

### Changed
- **`assembly:require_isotope_support` and `perf:stream_load` now default to TRUE.** Both were
  `registerFlag_`, so they defaulted to OFF -- while every benchmark figure this project has ever
  published passed them explicitly. The shipped default configuration had therefore never been
  measured. It has now been, and it was **6.9x slower and worse on both engines**: with
  `require_isotope_support` off, dataset D emits 1,255,577 spectra in 36:25 at 6.0x window-loop occupancy
  for Sage 12,212 / MSFragger 12,148; on, it emits 655,776 in 5:15 at 65.7x for Sage 12,482 /
  MSFragger 12,337. (An earlier measurement in July had put the gate's cost at −0.78% peptides; the
  sign flipped somewhere behind the charge gate, the apex estimator and the integer detector.)
  Without `stream_load` the loader holds the entire run at a measured 90.3 GB floor, which by itself
  exceeds the tool's current total peak. Both are now `<true/false>` options, since a flag cannot
  default to true; pass `false` to restore the old behaviour. **`perf:stream_load` also CHANGES
  OUTPUT** -- checked rather than assumed, and the assumption was wrong: true gives 655,776 spectra
  in 5:11 at 78.4 GB, false gives 656,371 in 10:30 at 121.8 GB, and the spectrum lists are not
  identical. The two readers pick peaks on different frame groupings. (Which one is closer to the
  truth was settled on 2026-09-05 -- see "Measured" at the top of this file: the one-shot reader is
  slightly better on a large file, at unchanged entrapment FDR, for 1.75x the memory.) Every figure
  this project has published used `true`. This is the same class of defect as
  the `ms2_noise_threshold_int` / `ms2_split_valleys` mismatch fixed on 2026-09-03 -- that sweep
  checked value-bearing options and missed the flags.

- **`charge:min_charge` now defaults to 2**, dropping singly-charged precursor hypotheses before
  extraction. They are ~30% of emission and ~1.7% of peptides, and removing them GAINS peptides on
  both search engines on all three benchmark files (Sage +0.8/+2.7/+1.4%, MSFragger +1.2/+5.4/+0.9%)
  while cutting emission ~30% and runtime 17%. Entrapment FDR stays inside the previous 95% interval
  on every file. Set `charge:min_charge=1` to restore the previous behaviour.

- **Shipped defaults now match the benchmarked configuration.** `trace:ms2_noise_threshold_int`
  100 -> 10 and `trace:ms2_split_valleys` 0 -> 7.0. Every published figure was produced with 10 and
  7.0 passed on the command line, so a default run could not reproduce any of them.
- **`trace:mz_estimator` is a real parameter** (default `apex`), replacing an undocumented
  environment variable that silently set every reported precursor and fragment m/z.

### Added
- Parameter bounds on the numeric options. A negative count used to wrap through `Size`
  (`min_fragments -1` became ~1.8e19 and nothing was emitted) and a negative tolerance gated
  everything out; both failed as an empty but valid-looking result.
- `bench/iso_dup.py`, `bench/iso_collapse.py`, `bench/iso_sim.py`, `bench/iso_loss.py`,
  `bench/iso_guard.py`: the isotope-duplication census, collapse-rule evaluation, content-similarity
  test and loss trace behind the emission findings.

### Changed
- **`trace:detector` now defaults to `integer`.** Mass-trace detection runs on the instrument's own
  integer axes -- flight-time bin for m/z, frame index for retention time -- following OpenMS
  MassTraceDetection step for step, and never materialises a PeakMap of doubles. On dataset D it is
  faster and uses ~40% less memory than the OpenMS path (dataset D 86 vs 150 GB, dataset A 126 vs 222 GB; wall
  6:21 vs 6:54 and 10:35 vs 11:14). The peptide effect is **mixed and not consistent in sign**: Sage
  -1.3% on dataset D but +0.2% on dataset A, MSFragger +2.1% on dataset D but -3.5% on dataset A, entrapment indistinguishable
  on both (dataset D 1.27% [1.03-1.52] vs 1.38% [1.14-1.63]). The detectors **agree on only ~85% of the
  union** of identified peptides -- 1,148 OpenMS-only and 981 integer-only on dataset D -- so the memory
  saving, not a peptide gain, is the case for the default. Without the vendor flight-time calibration
  the tool **falls back to `openms` and says so**; the emitted mzML records which detector ran as
  `spx:detector`.

### Performance (2026-09-04, round 2: three output-neutral changes)
Measured old-vs-new with both arms on the same host, three files concurrently on three machines,
one binary from a shared install. **Every arm is digest-identical to its base**, so peptides are
unchanged by construction and no search was run.

| | dataset D | dataset A | dataset B |
|---|---|---|---|
| RSS at end of the window loop | 65.0 -> **48.3 GB** | 85.6 -> **55.1 GB** | 70.0 -> **49.8 GB** |
| window-loop occupancy | 67.6 -> 69.9x | 56.2 -> 61.3x | 56.9 -> 63.3x |

**These are memory changes, not speed changes.** Replicated as interleaved base/clean pairs on
shared machines with the load recorded before each run -- three pairs on dataset D, two on dataset B, every
digest identical. dataset D: wall 321.3 -> 312.7 s (−2.7%), window-loop CPU 9,229 -> 9,690 (**+5.0%**),
RSS at end of the window loop 63,880 -> 46,922 MB (−26.5%), process peak 81.3 -> 77.9 GB (−4.1%),
all with non-overlapping ranges. dataset B: no wall-clock difference at all, same ~4% CPU cost, same
~30% memory drop. So the trade is **~27% of the window loop's memory for ~5% more CPU at roughly
unchanged wall time** -- worth taking because concurrency is bounded by the free-RAM admission
gate, but it is a memory change. A first pass of single pairs appeared to show 5-6% off wall on
three files; that was load drift on shared nodes and is withdrawn. Do not cite a speed-up.

- **Valley splitting is streamed.** It used to build every trace of a chunk as an OpenMS `MassTrace`
  and then split the whole chunk, so the input payload and the entire split output were alive at
  once -- the two largest structures after the peak slab, 7.5-10.6 GB. One trace is now in flight at
  a time and those lines fall to 0.08 MB. `ElutionPeakDetection` is constructed and parameterised
  once per chunk rather than once per trace.
- **The fragment RT gate rejects on 2 bytes instead of 8.** The gate discards 99.4% of what it
  visits -- 867 billion visits, 6.9 TB of traffic per dataset D run -- and was reading a `double` to do it.
  A parallel `uint16` bucket array does the reject; the exact test is unchanged and runs only on
  survivors, so the candidate set and its order are identical (the score-gate counters match to the
  digit between arms).
- **Trace arenas are reserved before merging.** `absorb()` appended without reserving into a
  just-emptied destination, so 12 band merges and 48 chunk merges regrew geometrically.
- Process peak RSS falls only 2-5 GB because the peak is set outside the window loop; the loop
  ENDING ~30% lighter is what governs how many windows fit in flight.
- Rejected on evidence: compacting the peak slab to above-noise peaks (the seeds line reads 100.0%
  -- there are no sub-noise peaks, so the slab is irreducible), and moving the seed sort into the
  band tasks (<1% of CPU, and it re-partitions seeds across bands, which cost 3.3% of peptides once).

### Performance (2026-09-04, trace representation)
- **A trace's profile is a span, not a point list.** A trace point has one degree of freedom -- its
  intensity -- so the three parallel 4-byte arrays (frame, intensity, flight-time bin) and their
  three allocations per trace are replaced by an entry frame and a range in a per-window arena, with
  zero meaning a missed frame. Verified **byte-identical on both detectors** on dataset D. The per-window
  correlation grid, which was a second copy of every profile, is gone. Window-loop memory falls from
  a 68.6 GB simultaneous peak to 27.9 GB (integer) / 32.4 GB (OpenMS); process peak RSS 105 -> 88 GB
  and 164 -> 147 GB.
- **The `SPEXTRACT_MEM_LEDGER` byte ledger has been REMOVED** (it was temporary scaffolding for the
  memory work above and is gone as of the round-2 changes; ~86 lines). While it existed it was
  corrected twice: the trace line had been cumulative rather than a peak, and a 6.0 GB
  `list<Peak2D>` line was a phantom (that list is per trace and dies each iteration). The structural
  wins it was built to find are permanent: the visited flags are a bitset (4.87 -> 0.61 GB), and the
  slab and per-window preparation arrays are released before valley splitting rather than held
  across it.

### Performance (2026-09-03, window loop)
- **The scoring gate reads its one field from a parallel array instead of striding the 96-byte
  trace record.** The gate rejects 99.4% of the fragments it visits; the scan was moving ~83 TB per
  dataset D run for 8 bytes of payload per visit. dataset D wall 12:53 -> **7:17**, scoring stage 47% -> 4.8% of
  window time, **peptide set identical** on Sage. Replicated three times.
- **The window loop is one OpenMP task pool over all threads** (master/worker) instead of a rigid
  threads/bands x bands grid; window-loop occupancy 65x -> 81x on 100 threads. The number of windows
  resident is bounded by the free-RAM admission gate, whose per-window estimate was re-calibrated
  (10x -> 16x compact bytes; it had been under-booking by 35-55%).
- **The retention-time axis is a global frame index.** A trace point is (frame index, intensity),
  8 B instead of 16, and aligning two profiles is integer equality instead of a binary search per
  point; the per-window correlation grid is built without collecting or sorting ~10^8 RT values.
- `trace:detector=integer`: mass-trace detection on integer arrays with the instrument's flight-time
  bin as the m/z axis, following OpenMS MassTraceDetection step for step, calibration re-applied only
  at export. Introduced here behind the flag at 6:53 wall / 105 GB against 7:17 / 164 GB for the
  OpenMS path; it became the default later the same day, on the numbers under "Changed" above.
  See docs/INTEGER-TRACING-DESIGN.md and docs/reviews/integer-tracing-2026-09-03/.
- Band-count sensitivity measured: 4/6/12/25 bands move the peptide set by −0.02% to +0.10% and
  differ in ~2% of spectra's content. Banding is an approximation, not exact as the help had
  claimed; the halo is now sized per edge.

### Fixed
- MS1 traces are filtered for non-finite coordinates before sorting. A NaN made the comparator
  non-transitive, which is undefined behaviour in `sort()` and silently corrupts every later binary
  search. The fragment path already did this; the MS1 path did not.
- `FragGrid` CSR offsets widened to 64-bit. At 8.5 M fragments and ~510 support points each, the
  32-bit offsets were within a factor of two of wrapping, and an overflow there silently corrupts
  every correlation rather than crashing.
- mzPeak input fails loudly instead of degrading: an unknown frame no longer falls back to the
  reference temperature while still being stamped as exactly calibrated; a short batch result from
  the reader no longer leaves empty tail frames; mismatched m/z, intensity and mobility array
  lengths are rejected instead of truncated to the shortest.

### Removed
- The raw apex-frame fragment backfill and the pre-extraction isotope collapse, with all their
  environment gates. Both were measured on all three files and falsified in every configuration;
  the results are in `docs/dataset D-BASELINE.md` and `docs/ISOTOPE-DUPLICATION-2026-09-03.md`, and the
  implementations are in the git history. 160 lines.
- An unused per-window band-count heuristic that a comment described as kept for reference.

### Runtime
- **MS1 path parallelised** (loader `[SpeXtract ms1-par]`: batched parallel MS1 decode with per-thread ZSTD contexts and
  `FrameCentroider`, ordered hand-off; consumer `flushMS1_`: parallel MS1 pick). dataset D/100 threads: LOAD 445 s -> 39 s,
  total 22:22 -> **15:21**; spectrum data byte-identical, Sage peptide set identical.
- Decode-once frame-major parallel MS2 loader (`[SpeXtract par-load]`), batched parallel MS2 pick; per-phase `[perf]`
  table and `[perf-load]` decode/hand-off/pick timers; `SPEXTRACT_LOAD_ONLY=1` for load-only profiling.
- Falsified and recorded: `MALLOC_ARENA_MAX=4` (10x slower at 100 threads), loader batch depth (64/256/1024 flat).

### Input
- **Streaming `.mzpeak` input** (`src/MzPeakStreamLoad.h`, build-optional `-DMZPEAK_ROOT`) on the mzPeak C++ library
  (OpenMS/mzpeak fork trunk); reads ims-compact archives with per-window mobility bands (works around name-only band
  params and NULL `precursor_index`). Benchmarked: not faster (the floor was the MS1 pick), and -11.7% peptides from the
  archive's two-point TOF transform (fragment error +10 ppm) -- see docs/RUNTIME-PLAN-2026-09-02.md and the converter
  handoff. Output stamps `spx:mz_calibration = mzpeak_two_point_transform ...` for this path.
- Library build recipe for GCC 13 / Arrow 23 (`scripts/build_mzpeak_lib.sh`, `patches/spx_ranges_to.h`,
  `patches/mzpeak-cpp-arrow23.py`).

### Determinism
- `bench/semantic_digest.py` now stops at `</spectrumList>`: the trailing `<indexList>` byte offsets shift with any header
  length change and had falsely flagged two spectrum-identical files. Loader batch size does NOT change content
  (922,902 spectra, 0 differ; Sage set identical), and **neither does the thread count**: threads 8 == threads 100 on two pairs with the fixed digest, retracting the 2026-09-01 "thread-count dependent" finding. Output is byte-identical spectrum data at any thread count and loader batch.

### Scoring / assembly
- **E5:** `apportion` and `rp_max` now go through the same emitted-intensity weights as the default path
  (`weighted_()`); they had silently bypassed `corr_power`/`im_weight`, so every earlier A/B of those flags was
  confounded. First clean apportion run: -3.9% peptides (falsified again; NNLS cut).
- **Reported m/z of every trace is now the APEX member's m/z** (`SPEXTRACT_MZ_ESTIMATOR=apex`, the new default;
  `mean` restores the OpenMS intensity-weighted centroid): dataset D +2.6% Sage (12,537 vs 12,217) and +4.1% MSFragger
  (11,927 vs 11,463) with the paired precursor mass error vs the reference implementation +1.12 -> +0.78 ppm. Both-engines gate passed;
  generalises: dataset A +4.4% (10,789), dataset B +2.7% (10,156); mean ratio vs the reference implementation on Sage 105.1% -> 109.2%;
  entrapment of the apex arm 1.28% [1.03-1.53] at nominal 1%.
- Under test (env switches, defaults unchanged): `SPEXTRACT_PICK_MZ_MODE=seed|top3` (pick-level m/z; the -2.4 ppm
  fragment offset vs the reference implementation is upstream of the trace estimator), `SPEXTRACT_DROP_PREC_ISO=1` (drop the precursor's
  M+1..M+3 from fragment lists; 28%/25% of spectra carried them -- FALSIFIED, -2.3%, keep off), `SPEXTRACT_MZPEAK_EXACT` (now the default; `=0` accepts the two-point m/z),
  `SPEXTRACT_MIN_ISOTOPES=k` (precursor gate by envelope depth; k=3 = -39% emission, -23% wall, MSFragger -0.9%, Sage -6.4%),
  `SPEXTRACT_PRECURSOR_LIST=<tsv>` (emit only precursors matching a reference list), `SPEXTRACT_BACKFILL_RAW=N` /
  `SPEXTRACT_BACKFILL_MAXFRAGS=K` (add the N most intense untraced raw peaks of the apex frame inside the IM band:
  N=50 lifts MSFragger to 90% of the reference implementation at -2.6% Sage; faint-tail-only variant under test).
- Output-identical: band partition by binary search (R3), early skip of fragments below `min_corr_pts` support.
- Help texts of falsified flags (`rp_max`, `consolidate`, `merge`) now say so; `charge:iso_im_tolerance` help
  no longer contradicts its default.

### Benchmark findings
- MSFragger gap is per-precursor CONTENT of the faint tail: same precursor population as the reference implementation -> 85.0%; emission -39% -> 84.8%;
  dt-only peptides are searched but short by ~2 of ~9 matched ions that never become a trace; raw apex-frame peak
  backfill recovers 1,151 of them (MSFragger 85.6% -> 90.0%). FragPipe rescoring (Percolator, MSBooster off for both) gives 13,211 vs the reference implementation
  15,947 peptides (82.8%; raw-hyperscore walk 85.6%) -- rescoring lifts both tools and leaves the ratio; ours has
  8.6 PSMs/peptide vs 2.0. Entrapment of the shipping arm: 1.28% / 0.95% / 1.32% at nominal 1% on dataset D/dataset A/dataset B.
- Falsified today (kept in docs so they are not re-proposed): apportion (clean, -3.9%), precursor-isotope stripping
  (-2.3%), median/seed/top3 m/z estimators, 48 MS1 trace bands (slower, not output-identical), MALLOC_ARENA_MAX=4.

### Docs
- `docs/BACKLOG-ANALYSIS-2026-09-02.md`: every open issue with status, root cause, fix, gate; reviewed by kimi + codex
  (`docs/reviews/`). Confirmed today: apportion/rp_max bypass corr_power (E5); 28% of spectra carry the precursor M+1 as a
  fragment (C9); false mzPeak calibration provenance (fixed).

## v0.2.0 — 2026-09-02

First release with a verified calibration, tests, and CI. Supersedes v0.1.0, which was tagged on the
day-one commit and never represented a working state.

### Headline: exact, license-free TOF -> m/z calibration
The Bruker `.d` reader converted TOF to m/z with a two-point linear-in-sqrt chord that is **-5 to -11
ppm biased** (m/z dependent). We derived the exact ModelType-1 model from the constants in each
file's own `MzCalibration` table and verified it to **2.5e-5 ppm** against Bruker's library:

    t_ns   = tof * DigitizerTimebase + DigitizerDelay
    C1_eff = C1 * (1 + dC1 * (T1_ref - T1_frame) / 1e6)
    t_ns   = C0 + (1e6 / sqrt(C1_eff)) * sqrt(m) + C2 * m

The widely-copied open implementation (timsrust-calibration, adapted then disabled by mzdata as "not
consistently better") **drops the `C2*m` term** -- worth -11..-40 ppm. That omission, not translation
subtlety, is why the port underperformed. Reported upstream.

Measured on three diaPASEF files, closed search: **+6-11% Sage peptide identifications**, no vendor
library. This is an identification-YIELD effect, not better spectra (emission moved -0.2%): mass
error is a feature of Sage's discriminant. Since the reference implementation always carried vendor-calibrated masses,
every earlier head-to-head was biased against the open path -- this removes a self-inflicted
handicap rather than establishing a lead.

### Benchmarks (3 files, both engines, identical settings; docs/BENCHMARK-MATRIX-2026-09-01.md)
| | dataset D | dataset A | dataset B | vs the reference implementation |
|---|---|---|---|---|
| Sage peptides @1% FDR | 11,976 | 10,333 | 9,891 | mean 105.1%, 95% CI [93.0, 117.3] -- consistent with PARITY |
| MSFragger peptidoforms @1% | 11,463 | 11,465 | 9,404 | mean 85.0%, 95% CI [76.7, 93.2] -- a supported DEFICIT |
| median precursor mass error | +1.5 | +3.2 | +2.3 ppm | the reference implementation -1.4 / +0.3 / -0.6 |
SpeXtract also emits ~1.32x more spectra, so per-spectrum efficiency favours the reference implementation. The MSFragger
gap survives perfect masses and is therefore search/detection-side.

### Safety
- **Fails closed.** ModelType != 1, `C2 <= 0` (a NULL/text C2 reads as 0.0 and selects the known-bad
  pure-sqrt law), `dC2`/`C3`/`C4` != 0, implausible `C1`, NaN, multi-row tables, or an unreadable
  `Frames.T1` are all errors, not silent approximations. `SPEXTRACT_ALLOW_CHORD_FALLBACK=1` (exactly
  "1") opts back into the biased chord.
- **Provenance.** Every emitted mzML records `spx:mz_calibration` =
  `tdf_table_modeltype1` | `bruker_sdk` | `legacy_chord_APPROXIMATE`, via an exported accessor that
  works across the shared-library boundary.
- Out-of-model inputs return NaN, never a plausible-looking 0.0.

### Tests and CI (first in this repo)
- `tests/test_calibration_cpp.cpp` -- 60 vendor-derived golden cases against the **shipped** header,
  round-trip inverse, and every rejection path. Ablation guards prove the suite can catch the
  known-bad implementation (39.4 ppm) and the temperature term (0.67 ppm).
- `tests/test_calibration.py` -- an independent re-implementation as a cross-check.
- `tests/test_entrapment.py` -- estimator CI, determinism, NaN-on-empty.
- `.github/workflows/tests.yml` runs all three. None needs a vendor library, cluster, or raw data.

### Determinism (scoped honestly)
Bit-reproducible **at a fixed thread count** -- verified by `bench/semantic_digest.py`, which compares
spectrum data and ignores the wall-clock stamp and the recorded parameters (two runs can never be
byte-identical otherwise; this also explains a July md5 mismatch left unexplained at the time).
Output DOES differ across thread counts, and last-ulp code changes cascade to ~2% in peptide counts.
The 0.06% noise floor in the decision rule is a SAME-BINARY figure.

### Also
- Entrapment estimator corrected (peptide-hypothesis ratio 0.6805, foreign-only counting, bootstrap
  CIs): actual FDR at nominal 1% is ~1.0-1.45%, not the ~2.3% previously reported.
- `assembly:corr_power=2` default (validated +8-10%, three files, both engines).
- Spectrum collapse falsified (-20..-60% peptides at every operating point, intrinsic to both tools).
- `scripts/apply_openms_patches.sh`; `scripts/sync_cluster.sh` refuses to deploy into an unpatched
  tree; evidence chain archived to durable storage.

### Known limitations
- Verified on ONE cohort/instrument; the three golden files share ONE MzCalibration vector.
- Public-data replication not yet done; every claim is single-cohort.
- A +1.5..+3.2 ppm residual remains, 2-4 ppm above the reference implementation on the same files. It survives every
  calibration path, so it is downstream of the m/z axis (centroiding / monoisotope reporting).
- Peptide counts are partly a mass-calibration measurement: compare only at similar residual ppm.
