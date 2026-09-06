# Runtime + quality assessment & prioritized plan (adversarial review: codex + vibe + kimi)

## Premise corrections the review forced (my earlier framing was partly wrong)
- **"SCORE 55% / TRACE 40%" is a WALL split, not compute.** The CPU-second instrumentation shows
  **TRACE ~69% / SCORE ~29%** of compute (kimi, confirmed by codex). SCORE dominates the *wall* only
  because it owns the Amdahl tail — so the fat to cut is in TRACE, not SCORE.
- **The wall tail is a BUG, not a law.** The dynamic team-apportionment the PERF doc claims exists is
  **absent**: when the fat window (window 0, 23.7% of precursors, ~2,000 s tail) is the last alive it
  still holds ~12 threads while ~90 cores idle. **My "43/128 cores / SCORE-55%" numbers are polluted
  by this bug** — the 19:30 baseline is measured against a ceiling that shouldn't exist.
- **The Pearson is already sparse** (CSR over each fragment's own support). Cutting SCORE compute /
  "sparse Pearson" is NOT a runtime lever (kimi + codex; vibe got this wrong). It IS a *quality* issue
  (the variance is over the full grid G≈1300 → G-dependent, inflates sparse junk on the faint tail).

## RUNTIME — prioritized (cheap waste first; task graph last)
| # | change | mechanism | impact | effort | risk |
|---|---|---|---|---|---|
| R1 | **Dynamic team apportionment** — retiring windows hand threads to live ones; recompute p_team at score-phase entry; make `omp_set_max_active_levels(2)` unconditional | kills the ~2,000 s near-serial tail on window 0 (the SCORE-55% wall) | **−20–30% wall (19:30→~14 min)** | 1–2 days | **Low** (slot-emit already race-free; require byte-identical mzML) |
| R2 | **Cull dead MS2 traces** — reject traces with support < min_corr_pts (median MS2 span 0.0 s; ≥50% can never pass the corr gate) before MTD→EPD→grid→score | TRACE is 69% of CPU; half of it is dead work | TRACE −30–50% CPU → **wall −15–25%** | hours | Low (prove recall-neutral ±0.3% + emission Δ≈0) |
| R3 | **Band-construction O(bands×peaks) fix** — the "partition once" loop scans every peak for every band; use two `lower_bound`s per band on the m/z-sorted spectra | −10–30% TRACE CPU | −5–12% wall | 0.5–1 day | Low (byte-identical) |
| R4 | dedup 2.99% double-scored precursors (window overlap) + grid merge-join | free compute | −2–3% | hours | Low |
| R5 | **Cross-window task graph** (phases as tasks, no window barrier) | the mapped ~2× lever; overlaps materialize/trace/EPD/score/emit | 19:30→12–15 min (subsumes R1 residual) | 1–2 weeks | High (determinism; PERF-ATTEMPT-1 died here) — only after R1–R4 |
| — | ~~sparse/bounded Pearson, SIMD SCORE~~ | already sparse | <2% | — | **Don't** |

**Realistic 2.5×→~1.3–1.5×: R1+R2+R3+R4 (days), then R5 (weeks). No single patch reaches it.**
3× compute gap ≈ ~1× intrinsic (share-all fan-out 6.45, open-search-safe emission — do NOT try to
remove, falsified 7×) + ~1.5× waste (dead MS2 traces, EPD SVD, grid, double-score) + ~0.5× OpenMS
constant factors.

## QUALITY — prioritized
| # | change | mechanism | impact | effort | risk |
|---|---|---|---|---|---|
| Q1 | **Fix the zero-padded Pearson variance** — compute over the union support, not G≈1300 implicit zeros | a documented-wrong statistic shipping as default; miscalibrated exactly on the sparse faint tail where the gap lives; upstream of the gate, rank key, and IM weight | +1–4% both engines (targeted) | 1 day | Low–med (retune min_correlation). CAVEAT: logoverlap [a different zero-pad fix] already fell −9%; test carefully |
| Q2 | **Correlation-power intensity weight** `inten *= c^k` (within-spectrum, share-all preserved) | engine-agnostic generalization of the IM weight (which was Sage-specific) → targets MSFragger's 79% | +2–5% MSFragger | 2 days | Med (kill if Sage regresses) |
| Q3 | **Native pre-centroiding aggregation** (running) → then decide 2D detector | the reference implementation's actual upstream signal enhancement; MS2-aggregation arm most likely to help MSFragger | 0–8% MSFragger | running / multi-week | the running test IS the experiment |
| Q4 | **Open-search correctness as a release gate** — charge abstention, monoisotope offset, mass-window; only `spx_guessed` is stamped though `open_search_safe` promises isotope-offset metadata | wrong charge/mono → lost IDs or false delta-mass structure; existential for SpeXtractor's actual use | little closed gain, big open-search validity | 3–5 days | Low algo risk |
| Q5 | export extraction evidence (IM residual, XIC corr, isotope quality) as a rescoring sidecar (MSBooster/Percolator) | bare engines discard the multidimensional evidence a 1D pseudo-spectrum drops | +3–8% workflow | 1–2 weeks | Med |

Honest ceiling without the reference implementation's upstream detector: ~90% Sage (near it), MSFragger ~79% with maybe
4–7 pts recoverable (Q1/Q2/Q4); the residual is detection-side. No information-theoretic 90% wall
(the reference implementation recovers it from the same raw data).

## SINGLE HIGHEST-ROI: R1 (dynamic team apportionment) — 2 days, −20–30% wall, output-identical.

## THE STRATEGIC MISTAKE (all three reviewers, unanimous)
**We are optimizing peptide parity with the reference implementation on CLOSED search — a corrupt metric.** Closed-search
peptide count structurally rewards over-emission (our worst failure; 73.5% of spectra pinned at the
500-peak cap IS over-emission), and some of the seven "falsified" fragment levers may have been
falsified on a metric that rewards the behaviour they were fixing. The binding metrics are:
(a) **entrapment-controlled open/blind PTM-site discovery rate** (SpeXtractor's reason to exist),
(b) **total pipeline time = extract + search** (924k fat spectra cost the search engine dearly),
(c) **quantification accuracy / FDR calibration**.
And the real differentiator is NOT open-search capability (the reference implementation does unrestricted/PTM/HLA too) —
it is that SpeXtractor is an **open, auditable BSD-3 implementation** independent of a closed JAR. The
product claim should be *"an auditable open diaPASEF extractor whose open-search FDR, charge/mono
correctness, PTM localization and quant are demonstrably trustworthy."* **Build the entrapment-
controlled open-search benchmark and re-baseline; until it exists, both 7:47 and 11,517 are the wrong
target.**

---
## R1 (dynamic p_team) FALSIFIED as implemented (2026-07-27)
Naive `p_team = max_threads/windows_live` computed independently per window at score-entry
OVERSUBSCRIBES: wall 19:30 -> 32:53, CPU ~43 -> ~31 cores. Windows enter the score phase at
different times, each divides the whole node by the current live-count, so several concurrent score
phases request more threads than exist -> thrash. The tail fix needs a GLOBAL thread budget (a real
task-graph / thread pool), not per-window division — kimi under-rated it (not 1-2 days). Reverted.
Runtime win remains the gate (-28%); the task-graph (R5) is the real remaining lever.

## Q3 (native pre-centroiding aggregation) — a QUALITY/RUNTIME TRADEOFF, not a clean loss
native_ms1_neighbors=1 (+ gate + IM 0.005): Sage 10,770 (best recall yet, +7.8% over gate 9,989),
entrapment target 9,929 @ 2.31% (vs gate 9,304 @ 2.33% — real gain, entrapment-confirmed) BUT
emission 2.48M spectra (2.7x — big runtime + search cost). So aggregation DOES surface real faint
precursors (contra the pure MS1-floor-failure prediction), it just over-emits. Worth only if the
recall matters more than the runtime; a tighter post-aggregation gate would be needed to keep the
gain without the 2.7x emission.

---
## Q1 + Q2 adversarial review (codex + kimi + vibe, all three, 2026-07-27) — UNANIMOUS
### Q1 (gate:variance_support = union-support Pearson): DROP — do not ship.
All three independently derived the SAME counter-example (precursor 8 pts, fragment 5, overlap 3):
G-Pearson = +0.47 (passes gate 0.3), union-Pearson = -0.50 (hard-rejected). Math is correct and
scale-invariant (NOT the logoverlap-v1 magnitude leak), but the statistic is SEMANTICALLY WRONG:
union-support Pearson penalises PARTIAL SUPPORT OVERLAP as anti-correlation -- it measures
co-LOCALISATION, not co-elution SHAPE. Faint fragments have only their apex above the trace
threshold (support censoring); union-Pearson reads the missing shoulders as real zeros against a
non-zero precursor = shape disagreement, and gates them out. It would WORSEN the low-abundance
recall gap it was meant to fix. Also: G-Pearson with G>>|support| has ~0 mean-correction, so it
already degrades to a robust cosine-over-intersection -- the [B0] "zero-padding defect" is largely
cosmetic, not load-bearing. Decision: DROP Q1 on unanimous analysis (the flag stays default-off in
code for reproducibility, marked falsified). codex also found + we fixed an option-selection bug
(var_support silently overrode logoverlap: now `var_support_ && !log_overlap_`).

### Q2 (assembly:corr_power = emit intensity *= c^k): clean code, mechanism hostile.
All three validated the emit-only implementation (cap key stays base*c; NOT rank_by=intensity) and
all three predict it HURTS MSFragger (-2 to -5%): c is ABUNDANCE-CORRELATED (faint precursor ->
noisy XIC -> low c near the 0.3 floor), so scaling emitted intensity by c^k attenuates exactly the
faint fragments the recall gap is made of. Contrast im_weight_sigma (+4.1% Sage): it down-weights
by IM distance, an axis ORTHOGONAL to abundance (interferent purity), so it BOOSTS true faint
fragments; Q2's axis is ALIGNED with abundance, so it suppresses them. Expected pattern: small/flat
Sage, negative MSFragger. Running q2k1 (k=1) + q2k2 (k=2) as the empirical confirmation the
reviewers themselves proposed (base + Q2-k1 was the one arm all three wanted); non-zero surprise
chance is why it's worth the single decisive arm rather than dropping unbenched.
