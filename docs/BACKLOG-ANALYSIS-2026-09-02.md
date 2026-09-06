# Backlog analysis — 2026-09-02 (revised after kimi + codex adversarial review)

Every open item across `docs/` plus today's findings, reduced to: **what the issue is, whether it is confirmed, the root
cause, the fix, the cost, and the gate.** Status vocabulary: **MEASURED** (a number on a named file and build),
**CONFIRMED** (seen in code/data), **HYPOTHESIS** (argued, not measured), **FALSIFIED** / **CLOSED** (kept in the register
at the end so nobody re-proposes or re-does them).

Review corrections applied (kimi 9 findings, codex 13): attribution-by-subtraction is now labelled as such (A1, C1);
mechanisms without a measurement are HYPOTHESIS (B3, C2, D1); mis-stated statuses fixed (C7 was measured, apportion was
falsified — under the E5 confound, D5 is worked around); R2 re-scoped to the byte-identical location; missing rows added
(R4, BC-2, BC-3, OI-2..7, PTM ambiguity, hardening obligations, Q4); duplicates and closed items moved to the register;
codex #13 (false calibration provenance on mzPeak output) **fixed in code today**; every plan step now carries a gate.

## A. Runtime (dataset D, 100 threads; HEAD build with timers; total 22:22)

| # | issue | status | evidence | root cause | fix | gate / cost |
|---|---|---|---|---|---|---|
| A1 | **~429 s of LOAD (444.7 s) is outside the MS2 batched loop** — attributed by control flow to the MS1 path | MEASURED residual (444.7 − 4.4 − 11.4 = 428.9 s); the component split (decode / sort / pick) is **not** timed | dataset D P256det; dataset A LOAD-only: 418 s with MS2 decode 10 s + hand-off 25 s → same shape on a second file | `loadDIAStreaming` decodes + sorts 1,343 MS1 frames (7.7% of frames, 37% of peaks) serially; `consumeMS1Spectrum_` picked each serially; ranked "small" by frame count | loader: batched parallel MS1 decode (`[SpeXtractor ms1-par]`, per-thread ZSTD ctx + `FrameCentroider`, applied 14:45); consumer: parallel MS1 pick `flushMS1_` (in tree) | **DONE 15:01, gate PASSED:** LOAD **38.8 s** (80.5×; MS1 decode 6.9 s, hand-off 15.4 s, pick 27 s), total **15:21 (−31%)**, spectrum data byte-identical to P1 (`ca1ddbbd…`), Sage 12,217 = 12,217 |
| A2 | mzPeak input does not help runtime | **CLOSED 15:23**: with A1 landed, `.mzpeak` LOAD 104 s vs `.d` 39 s; total **17:04 vs 15:21**; CPU 8,392 vs 3,122 s; system time 1 h 04 min = **52.9 M voluntary context switches** (vs 0.69 M from `.d`; same page faults) — futex hand-offs inside the library's decode path (reader mutex / executor pool), a library-side fix (synchronous decode entry point) | format-independent floor removed; the remaining mzPeak cost is the reader itself | keep `.d` primary; mzPeak stays build-optional | closed |
| A3 | WINDOW_LOOP 774 ± 30 s at 59.5×; PRECURSOR_INFER 31 s serial; MS1_TRACE 61 s at 5.1× | MEASURED (±30 s run-to-run on serial phases) | `[perf]` tables, 3 runs | uneven windows; serial precursor inference; MS1 tracing at low width | (i) **R2 re-scoped**: short-trace skip only at the dot-product-only location (byte-identical there; the pre-grid cull changes `g.rt` → every Pearson, WINDOW-LOOP-AUDIT §5); (ii) R3 two-`lower_bound` band construction; (iii) precursor inference is a GREEDY sequential seed walk (`inferPrecursors_`: seeds in score order claim isotope partners; :1089+) — parallelising it changes output; 31 s = 3% of the run, leave it; (iv) `perf:ms1_trace_bands=48` **measured 18:31: slower (MS1_TRACE 144 vs 61-80 s) and NOT output-identical (922,899 vs 922,902 spectra)** — the MS1 band partition is not exact; keep 12, item dropped | byte-identity digest for (i)–(iii); (iv) needs the quality gate |
| A4 | Peak RSS 107–124 GB | MEASURED | RSS 29 GB after load, 81–89 GB loop end | temporal high-water mark: per-window working sets + ~10 GB buffered output; July: ~87% glibc arena retention | (i) allocator A/B **measured 17:57**: tcmalloc_minimal byte-identical, peak RSS 101.8 → 94.3 GB (−7%), wall unchanged (15:49), system time halved — a memory lever only; full tcmalloc 97.5 GB (worse than minimal): `bench/run_joint.sh` already preloads tcmalloc_minimal (line 80) — the cluster prof/ scripts did not; use it for RSS-sensitive runs; (ii) stream the output; (iii) R2 (dot-product location) | (ii) 1–2 days |
| A5 | R4: 2.99% double-scored precursors (window overlap) | MEASURED (07-27) | RUNTIME-QUALITY R4 | overlapping windows both score the same precursor | dedup at admission (hours) | not byte-identical — quality gate (emission Δ≈0, set overlap) |
| A6 | OI-2..OI-7 (MassTrace move ctor, unconditional span diagnostics, RT SoA, `g.rt` capacity, `std::move(tmp_spec)`, reserves) | CONFIRMED still open in code (codex: span calc unconditional :1000–1017, `g.rt` capacity :699–706) | WINDOW-LOOP-AUDIT §3 | — | each hours; OI-2/OI-6 need a libOpenMS rebuild | byte-identical digest |
| A7 | R5 cross-window task graph | HYPOTHESIS | "15–30% of what remains" | window barrier | only after A1/A3, behind the determinism gate | 1–2 weeks; deferred |

## B. Determinism / reproducibility

| # | issue | status | evidence | root cause | fix / claim | gate |
|---|---|---|---|---|---|---|
| B1 | Byte-identical at fixed config | MEASURED 11:36 | S1=S2=P1=P2 | — | README: "deterministic at fixed thread count and fixed loader batch" | closed |
| B2 | Loader batch 64 vs 256 | **CLOSED 14:46 — content identical.** specdiff: 922,902 spectra, 0 differ; Sage set 12,217 = 12,217, symdiff 0; fixed digest identical (`ca1ddbbd…`) | the earlier "DIFFERS" was `bench/semantic_digest.py` hashing the trailing `<indexList>` byte offsets, which shift when the header changes by 4 bytes (build-date version stamp) | digest now stops at `</spectrumList>` (fixed today); both reviewers' "B2 before A1" concern is answered by data, and their protocol (1-thread 64 vs 256, MS1-stage digests first) stays the standing procedure for any future batch/threads claim | closed |
| B3 | Output differs across thread counts (8 vs 100) | **CLOSED 15:14 — it does not.** With the fixed digest, threads 8 == threads 100 on two independent pairs (`final_t8`/`final_t100`: `2a1dc5fd…`; `sdkcal_det8`/`sdkcal`: `0cacb566…`); per-spectrum: 923,713 compared, 0 differ | the 09-01 "thread-count dependent" finding was the same `<indexList>` offset artefact as B2 | README/CHANGELOG: **deterministic across thread counts** (byte-identical spectrum data); BC-5 index tie-break is not needed for this | closed |
| B4 | Cross-binary sensitivity (~2% from a 1e-11 change) | MEASURED once (one perturbation) | +1.9% peptides | greedy thresholds amplify last-ulp differences; **not a noise distribution** — one observation | decision rule stays: set overlap + ppm median beside every count; do not quote ±2% as a bound | policy |

## C. Identification quality (the science)

| # | issue | status | evidence | root cause | next test (with falsifier) | cost |
|---|---|---|---|---|---|---|
| C1 | **MSFragger 85% of the reference implementation (CI [76.7, 93.2])**; Sage at parity (105%, CI [93, 117]) | MEASURED | 3 files × 2 engines | **10% of the gap = wrong charge labels (measured, override_charge)**; H7 (charge discard) falsified; H1 (null broadening) **not** falsified by the z∉{2,3} deletion (codex: lost true targets can exceed the threshold benefit); the remainder is UNRESOLVED — not "tail quality" by subtraction | discriminating arms, in this order: (a) matched-precursor score distributions on shared peptides — **already measured (BENCHMARK-MATRIX H3): on the 9,405 shared peptides, matched ions 12 vs 12 and hyperscore 23.4 vs 23.4, i.e. identical; the deficit lives entirely in the 3,599 the reference implementation-only peptides**; (b) rescoring A/B (MSBooster/Percolator) on BOTH tools — read the spx/dt RATIO after rescoring, not "half closed" — **run 21:09–21:30**: FragPipe 24.0 headless; MSBooster (DIA-NN predictions, with or without spectra) aborts on OUR pepXML ("Prediction missing … REM[ox]DQTM[ox]AANAQK|3"), completes on the reference implementation's: **dt rescored = 18,670 peptides / 37,086 PSMs at 1%** (vs 13,932 raw-hyperscore peptidoforms, +34%); **Percolator-only, like-for-like (21:46): apex 13,211 vs dt 15,947 = 82.8% (raw walk 85.6%) — rescoring lifts both, the ratio does not move; ours 8.6 PSMs/peptide vs 2.0** → the deficit is content (faint-tail coverage / emission competition), not scoring; a MSBooster-tolerant rerun needs the offending PSM class handled (cap var mods for both, or a predictor that emits it); the "07-26 arm" was Sage chimera only (96.0% ratio on the July binary); (c) emission-controlled arm (quality-gate 924k→~700k + precursor-matched sub-arm, step 02); (d) oracle-fragment arm (kimi F7) to separate metadata / content / search | (a) hours; (b) days; (c) 1 week; (d) days |
| C2 | Residual +1.5..+3.2 ppm on identified precursors, 2–4 ppm above the reference implementation on every file | **MEASURED PAIRED 15:20, 3 files** (dataset D-BASELINE): on shared (peptide, charge) pairs our precursor m/z is +2.63 (dataset A) / +1.48 (dataset B) / +1.12 (dataset D) ppm above the reference implementation's — z=2-dominated, absent at low m/z — and our fragment m/z is **−2.0..−2.4 ppm below, uniformly** — selection bias excluded | two reporting conventions differ from the reference implementation's (MS1-trace m/z estimator high, larger at z=2 and mid/high m/z; MS2-trace estimator low), not a calibration offset (that would move both the same way) | **A/B done 17:22**: apex estimator **+2.6% Sage peptides** (12,537 vs 12,217, same binary), precursor bias +1.12 → +0.78 ppm; median −1.0%; the fragment −2.4 ppm does not move (it is upstream: IM-cluster pick centroid). **MSFragger gate PASSED 18:42: apex 11,927 vs 11,463 (+4.1%)** → apex is the new default; **3-file confirmed 20:57**: dataset A +4.4%, dataset B +2.7% (mean ratio vs the reference implementation Sage 109.2%), entrapment 1.28% [1.03–1.53] | **pick-level A/B done 20:00: seed −2.7%, top3 −1.7% (falsified), fragment offset unchanged (−2.4 ppm) → not in our estimator chain; the reference implementation-side reporting; closed** |
| C3 | Entrapment FDR of the shipping table-model arm | **MEASURED 16:41 (dataset D): 1.37% [1.10–1.64] at nominal 1%** (11,489 targets / 107 entrapment), vs the vendor-SDK arm 1.26% [1.00–1.51] — CIs overlap, the calibration gain is not bought with FDR | — | dataset A 0.95% [0.72–1.18] (dt 1.24%), dataset B 1.32% [1.04–1.59] (dt 1.22%) | **done, 3 files** |
| C4 | Emission 1.32× the reference implementation; 4.6 PSMs/peptide; purity 5.8% | MEASURED | dataset D | share-all fan-out by design; all cross-precursor redistributions FALSIFIED ×7 | only a precursor QUALITY gate (C1 c) or within-spectrum weighting; never merging/collapse | inside C1 c |
| C5 | Charge agreement 74.6% vs 82.6%; ~44% of within-cycle duplicates ±1 isotope | MEASURED (dataset A) | — | partner-count charge; mono selection | learned charge/mono predictor (BACKLOG N+1) **only if** the oracle arm (C1 d) isolates a metadata-bound residual; held-out dataset B; must beat 74.6% by > 0.17% | 1–2 weeks, conditional |
| C6 | Native pre-centroiding aggregation | **MEASURED 07-27 (pre-calibration, Sage only)**: `native_ms1_neighbors=1` Sage +7.8%, entrapment 2.31%, at 2.7× emission | RUNTIME-QUALITY Q3 | gain is real but bought with emission | re-test under corrected calibration, both engines, WITH a post-aggregation gate that holds emission ≤ 1.4×; falsifier: peptides@matched-entrapment ≤ baseline | 1 day of runs; only with the gate |
| C7 | `assembly:apportion` / NNLS (BACKLOG 1) | **FALSIFIED CLEANLY 16:46** (E5 fixed, corr_power applied): apportion=1.0 → 11,742 vs 12,217 share-all (−3.9%, symdiff 1,234/759) | cross-precursor apportionment loses peptides, with or without the weights | **NNLS cut**; do not re-propose cross-precursor redistribution | closed |
| C8 | Open-search: step 01 done on dataset D; step 02 emission-controlled arm not run; PTM ambiguity 34.7% vs 14% known-PTM bins and the 1,610-unit Amidation bin uninvestigated | CONFIRMED | dataset D-BASELINE §303–311 | — | step 02 as pre-registered (v3 scorer, full FDR curves); one script on existing output for the Amidation bin | 1 week; hours |
| C9 | Precursor isotopes leak into fragment lists | **MEASURED (28.4% M+1, 24.8% M+2) and the fix FALSIFIED 19:05**: dropping M+1..M+3 → 12,249 vs 12,537 (−2.3%, same binary) while contamination fell to 0.2% | d12_dropiso | the survivors are useful to the engine, not noise | keep them; `SPEXTRACTOR_DROP_PREC_ISO` stays off | closed |
| C10 | Public-PXD replication | not done | patient data unreproducible in principle | — | the reference implementation paper's public diaPASEF data, both tools, both engines | 1–2 weeks; gate for any public number |

## D. mzPeak / converter / reader

| # | issue | status | evidence | root cause | fix | owner |
|---|---|---|---|---|---|---|
| D1 | mzPeak input: −11.7% peptides, fragment error +10.05 vs +3.39 ppm | **ISOLATED 18:58**: with the exact calibration recovered from the tdf (same archive, same reader, same binary) **12,082 peptides (−1.1% vs `.d`)**, paired ppm vs the `.d` run +0.00/+0.00 → the two-point transform was **91%** of the loss; the residual ~1% is the IM axis | d10_MZPexact | archive carries only `mzpeak:transform_params` | exact path now DEFAULT and fail-closed (tdf embedded or sidecar); converter ask (handoff §5) unchanged | closed for SpeXtractor |
| D2 | `--ims-chunked` aborts on TDF | CONFIRMED | dataset D + dataset A | facet-layout mismatch in the writer | converter (handoff §1) | converter |
| D3 | archive 4.5–5.6% larger than `tdf_bin` | MEASURED | — | per-peak f64 1/K0 (0.57) + tof (0.52) | converter question (handoff §2) | converter |
| D4 | band as name-only CV params; `precursor_index` NULL → reader attaches both ions to precursor 0 | CONFIRMED | `spx_band2` | converter metadata + reader attach key | worked around (flatten + match by m/z, 0.05 Da); converter/reader fix (handoff §3–4) | converter + library |
| D5 | library targets GCC 14 / Arrow 24 | **WORKED AROUND** (builds on GCC 13 / Arrow 23 with `patches/spx_ranges_to.h` + `patches/mzpeak-cpp-arrow23.py`) | build log | `std::ranges::to`, `ReadTable()` Result form | upstream as `#if` guards | library |
| D6 | mzPeak output stamped `BrukerTimsFile::lastMzCalibration()` (codex #13: false provenance) | **FIXED today**: mzPeak path stamps `spx::lastMzPeakCalibration()` = "mzpeak_two_point_transform … archive=…" | src :2753 | one global for two loaders | **verified 15:23** in d5_MZPms1/pseudo.mzML (`spx:mz_calibration = mzpeak_two_point_transform (...) archive=S30_v092_archive.mzpeak`) | done |

## E. Engineering / product hygiene

| # | issue | status | fix | cost |
|---|---|---|---|---|
| E1 | determinism smoke test in CI | open (tests + CI exist since v0.2.0) | fixed-config digest on a public-derived fixture with the fixed digest script | 1 day |
| E2 | README 13 commits stale; branch 59 ahead of `main`; binary still `spextractor` | CONFIRMED | README rewrite (dataset D-BASELINE + runtime + mzPeak); fast-forward main; rename executable (compat symlink) | 1 day |
| E3 | falsified flags' help still recommends them (`rp_max` :764, `merge` :796 "summed" vs MAX, `consolidate` :803) | **FIXED 15:15**: all three marked FALSIFIED in help; merge help now states MAX-by-default (`merge:sum_intensity` sums) | — | — |
| E4 | `charge:iso_im_tolerance` default 0.05 contradicts its own help (:840) | CONFIRMED | decide 0.01 behind an dataset D gate, or fix the help | hours + 1 run |
| E5 | `apportion` / `rp_max` bypass `corr_power` / `im_weight` — `assembleOne_` applies them (:1632–1643); apportion emits raw shares (:2345), rp_max raw intensities (:2454) | **CONFIRMED today; FIX APPLIED 15:05** (`weighted_()` is the one place; all three branches use it); invalidates every apportion/rp_max A/B ever run (they compared corr_power 2 vs 0) | **gate PASSED 15:58**: default path byte-identical to P1 (digest) and Sage set identical (12,217); the first apportion=1.0 run with corr_power applied (C7) follows | hours; then C7 |
| E6 | macOS serialises the window loop | NOT VERIFIABLE here (no Mac build of the tool exists; no `__APPLE__` guard) | CMake FAIL without OpenMP unless `-DALLOW_SERIAL=ON`; one Mac build | 1 h |
| E7 | evidence chain / provenance | **tarballs on ceph**: `prof-2026-09-02.tar.gz` (718 MB, c6b2e6b8…, runs to 16:45) and `prof-2026-09-02b.tar.gz` (1.35 GB, f91ed5dd…, runs 16:45–21:05 incl. apex/entrapment/D1/pick/iso/MSFragger TSVs), `prof-2026-09-02c-rescoring.tar.gz` (175 MB, 382a4167…, FragPipe rescoring outputs) | stamp src sha in the mzML (open) | hours |
| E8 | OpenMS patch drift (3 incidents) | process | regenerate from `.orig` + grep-verify (commit template); CI check that the patch applies to the pinned tarball | hours |
| E9 | hardening obligations from dataset D-BASELINE §332–334: benchmark must fail closed on a missing/wrong SDK path; opentims discards `tims_index_to_mz`'s return code | open | fail-closed check in `bench/`; return-code check in the vendored converter | hours |
| E10 | Q4 mono-offset stamp | **CLOSED**: `spx_n_isotopes` (isotope-offset hypothesis) and `spx_guessed` are stamped per spectrum (src :1610–1616) | — | — |

## Recurring root causes

1. **Ranking by the wrong denominator** (A1 frame-vs-peak share, "write is the cost", "MS1 loop is small"): every mis-ranking was an inferred share; every correction was a timer. No runtime measure is ranked without a `[perf]` line, and residuals are labelled residuals.
2. **Instruments that lie**: the semantic digest hashed byte offsets (B2), the `[t=..]` markers hid the load (09-01), the raw-hyperscore harness is mass-error-blind (C1). Validate the instrument before trusting a "FAIL".
3. **Confounded A/Bs**: apportion/rp_max ran with corr_power silently off (E5); native aggregation was compared at 2.7× emission (C6); mzPeak changes m/z and IM at once (D1). Every arm needs one variable and a pre-registered falsifier.
4. **Calibration is the yield lever** (C2, D1): +6–11% from the exact model, −12% from the two-point one; the residual +2–3 ppm is still ours and unexplained.

## Plan (reviewed; ordered by value/effort; each step with its gate)

**Week 1 — land the runtime win, restore the instruments, cheap science first**
1. ~~A1~~ **DONE** (LOAD 39 s, total 15:21, byte-identical, set identical); patch regenerated from `.orig`, reverse dry-run verified, committed.
2. ~~B3~~ **CLOSED**: threads 8 == threads 100 with the fixed digest (two pairs). Ship the claim as "byte-identical spectrum data at any thread count and loader batch".
3. ~~C3~~ **MEASURED on 3 files**: dataset D 1.37%, dataset A 0.95%, dataset B 1.32% at nominal 1%; the reference implementation 1.22–1.24%. Parity claim now FDR-controlled.
4. ~~C2 paired bias~~ **MEASURED**; ~~estimator A/B~~ **apex = +2.6% Sage** (MSFragger arm queued before it becomes the default); fragment bias is at the pick level (next A/B).
5. ~~E5/E3/E4~~ **DONE** (identity gate passed 15:58); ~~C7 apportion re-run~~ **FALSIFIED CLEANLY** (−3.9%): NNLS cut from the backlog.
6. ~~A3(i)/(ii)~~ **applied, byte-identical (17:41), no measurable wall gain** (WINDOW_LOOP 763 vs 759 s; RSS −4 GB); A4(i) allocator A/B running.
7. ~~E7 tarball~~ **DONE** (718 MB on ceph, sha c6b2e6b8…).

**Week 2 — resolve the MSFragger attribution instead of naming it**
8. ~~C1(a)~~ done (identical scores on shared peptides); ~~C1(b)~~ **done: rescoring does not move the ratio (82.8%)** → the deficit is content-side. Next: C1(c) emission-controlled arm (step 02) — **sub-arm 1 DONE 22:34: k=3 cuts emission −39% (below the reference implementation's 700k) and the MSFragger ratio stays 84.8% (was 85.6%) while Sage loses 6.4% → emission competition refuted; the deficit is per-precursor faint-tail content**; **sub-arm 2 DONE 22:49: with the reference implementation's own precursor population (432,579 spectra) the MSFragger ratio is 85.0% → per-precursor CONTENT; **content candidate 2 measured 01:28: on dt-only peptides we have a PSM for 64% with 7 vs 9 matched ions (deficit +2 at every charge; shared peptides 12 vs 12) → under-covered faint tail; **coverage arms FALSIFIED 03:10 (min_corr 0.2: no change; noise 5: byte-identical) → the missing ions are not traced at all; **raw apex-frame peak backfill N=50 CONFIRMS the content hypothesis 05:39: MSFragger 11,927 → 12,539 (85.6% → 90.0% of dt; 1,151 dt-only peptides recovered) at −2.6% Sage; N=150 too noisy; faint-tail-only (≤K fragments) backfill 07:11: K=150 = Sage +1.5% / MSFragger +0.6% (clean, small), but the dt-only deficit is untouched — the MSFragger +5% lives in the rich spectra; precursor-intensity-gated backfill arm queued;** **funnel: of the never-searched third, 46% have a hypothesis (assembly loss), 35% none (MS1 detection)**; **corr_power=0 FALSIFIED 01:52: Sage −8.3%, MSFragger −4.2% — corr_power=2 helps both engines**; C1(d) oracle-fragment arm after; C5 stays conditional on (d).
9. ~~D1 isolation A/B~~ **DONE**: exact calibration recovers 91% of the loss (12,082 vs 10,785; `.d` 12,217); exact path is the mzPeak default.
10. ~~C9 isotope exclusion~~ **FALSIFIED** (−2.3%): the precursor isotopes in the list help the engine.

**Weeks 3–4 — the arms that change the product**
11. C1(c)/C8 step 02 emission-controlled open-search arm (pre-registered, v3 scorer, full FDR curves).
12. C6 native aggregation re-test with the emission gate (1 day of runs).
13. C10 public-PXD replication + v0.3.0 (E2 README rewrite, rename, main fast-forward, E1 determinism CI).

**Cut / not doing:** merging/collapse; NNLS unmixing (its prerequisite failed cleanly); arena capping; ring buffers/reader processes; the 2D Gaussian detector before a demonstrated content gap; NNLS before apportion re-measurement + P conditioning; the learned predictor before the oracle arm; closed-search parity as a goal; ±2% as a generic uncertainty bound.

## Register: closed / falsified today (do not re-open)
- **B2** loader batch size does not change content (specdiff 0/922,902; Sage set identical; fixed digest identical). The digest tool was wrong, not the loader.
- **A5** `MALLOC_ARENA_MAX=4` at 100 threads: ~10× slower window loop; aborted.
- **A6** loader batch 64/256/1024: LOAD 436/432/434 s; keep 256; no ring buffer or reader process.
- **B1** byte-identical at fixed config (4 runs); parallel pick == serial pick.
- **D6** false mzPeak calibration provenance — fixed in code.
- **E1 (tests/CI)** exist since v0.2.0; only the determinism smoke test remains.
- OpenMS in-tree `MzPeakFile` (0 spectra on ims-compact) — closed for SpeXtractor (bypassed); open only as a separate OpenMS product.

## Status at 20:10 (end of the first day of the plan)

Done and committed today: A1 (MS1 path, 22:22 → 15:21, byte-identical), B2/B3 (both retracted — digest artefact; output is
byte-identical at any thread count and batch), C2 (paired bias measured on 3 files; apex estimator +2.6% Sage / +4.1%
MSFragger → new default; pick-level modes falsified; the −2.4 ppm fragment offset is the reference implementation-side reporting), C3
(entrapment 0.95–1.37% on 3 files), C7 (apportion falsified cleanly; NNLS cut), C9 (isotope exclusion falsified, −2.3%),
D1 (two-point transform = 91% of the mzPeak loss; exact path recovers it and is the default), D6 (provenance), E3/E4/E5/E7/E10,
R3 + early skip (identical, no gain), A4(i) (tcmalloc_minimal −7% RSS, wall-neutral), A3(iv) (48 MS1 bands: worse, not identical).

Running: apex on dataset A/dataset B + entrapment of the apex arm (`bench_apex23.sh`), then the FragPipe MSBooster/Percolator
rescoring A/B on both tools (`bench_rescore.sh`, C1 b).

Next (in order): read the rescoring ratio → decide between the emission-controlled arm (step 02) and the oracle-fragment
arm for the MSFragger attribution; C6 native aggregation re-test WITH an emission gate; E1 determinism smoke test in CI
(needs a public-derived fixture); C10 public-PXD replication + v0.3.0 (README rewrite, executable rename, `main`
fast-forward). Not doing: cross-precursor redistribution of any kind, MS1 band count changes, isotope stripping,
median/seed m/z estimators.
