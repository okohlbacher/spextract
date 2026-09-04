# dataset D is the decision file

**All progress/regression decisions are made on dataset D alone.** Other samples are for confirmation
once something has already won here — cross-sample comparison has produced three wrong
conclusions in this project, and single-sample iteration removes that class of error entirely.

## Reference numbers (frozen 2026-07-22)

| measure | value | how |
|---|---|---|
| **peptides @ peptide_q ≤ 0.01** | **10,072** | `sage_deiso.json`, current defaults |
| peptides, common-FDR procedure | 9,977 | DEPRECATED 2026-09-01 (procedure's primary outputs no longer exist) |
| the reference implementation, same procedure | 11,552 | DEPRECATED (doc-only; canonical ref = 11,517) — legacy ratio 86.4% not citable |
| MSFragger, same file | 10,470 | DEPRECATED (the reference implementation 13,014 doc-only; canonical = 13,932 form-level / 12,752 seq-level) |
| PSMs/peptide | 8.98 | the reference implementation 1.66 |
| within-cycle multiplicity | 2.85 | the defect |
| across-cycle spread | 2.66 | **correct** — MS1 FWHM is 2.61 cycles |
| MS1 FWHM | 3.61 s | flat with RT (1.04× across gradient) |

**Current defaults** (registered defaults, verified against the code 2026-09-03 --
`ms2_noise_threshold_int` and `ms2_split_valleys` had been benchmarked at 10 and 7.0 while the
registered defaults were still 100 and 0, so a default run could not reproduce any figure in this
file; the defaults were moved to the benchmarked values): `trace:ms1_split_valleys 7.0`,
`trace:ms2_split_valleys 7.0`, `trace:ms2_noise_threshold_int 10`, `trace:ms2_chrom_peak_snr 1.0`,
`trace:ms2_min_length_sec 0`, `trace:mz_estimator apex`, `charge:scoring count`,
`charge:min_charge 2`, `assembly:corr_power 2`, `perf:trace_bands 12`,
`assembly:require_isotope_support true`, `perf:stream_load true`, `trace:detector integer`.

**The 2026-09-03 defaults sweep only checked VALUE-BEARING options, not FLAGS, and two flags that
every benchmark passes were left defaulting OFF** -- found 2026-09-04 and fixed the same day. Both
are now `<true/false>` options defaulting to `true`, because `registerFlag_` cannot default to true.

`assembly:require_isotope_support` measured on dataset D for the first time in either direction:

| | spectra | wall | window-loop occupancy | Sage | MSFragger |
|---|---|---|---|---|---|
| **true** (every benchmark on record) | 655,776 | **5:15** | 65.7x | **12,482** | **12,337** |
| false (what the shipped default did) | 1,255,577 | **36:25** | **6.0x** | 12,212 | 12,148 |

**The shipped default was 6.9x slower AND worse on both engines.** The only prior measurement
(2026-07-26) put the flag's cost at −0.78% peptides; today it GAINS 2.2% Sage and 1.6% MSFragger.
The sign flipped somewhere between the charge gate, the apex estimator and the integer detector,
which is why it was re-measured rather than inherited. Note the occupancy collapse to 6.0x with the
gate off: 2.15 M precursor hypotheses instead of 1.10 M puts the window loop at 91.8% of the run
and it barely parallelises -- so the guessed singletons cost far more than their share of the work.

`perf:stream_load` is the frame-by-frame reader; with it off, `loadExperiment()` holds the whole run
at a MEASURED 90.3 GB floor at t=0, which alone exceeds the tool's entire current peak (80 GB).
Every figure in this file was taken with it on.

**And it is NOT output-neutral, contrary to the assumption every benchmark has rested on.** Measured
2026-09-04 on dataset D, both arms on one node:

| `perf:stream_load` | spectra | wall | peak RSS | digest |
|---|---|---|---|---|
| true (shipped, and what every figure here used) | 655,776 | 5:11 | 78.4 GB | -- |
| false | 656,371 | 10:30 | 121.8 GB | **DIFFERS** |

0.09% more spectra from the one-shot reader. The two paths pick peaks on different frame groupings,
so this is a second detector-like choice hiding in an option named for performance. **Which reader
is closer to the truth is UNMEASURED** -- neither arm has been through a search engine. What can be
said is that the shipped default is the arm every number in this file was taken with, and that the
alternative costs 2x the wall and 1.6x the memory. The help text now says CHANGES OUTPUT.

## Measurement methodology (2026-07-28) — variance model + estimator fix

* **BOTH engines are DETERMINISTIC: run-to-run sigma = 0.** Re-searching one mzML (dataset D corr_power=2)
  3x: Sage 10,802/10,802/10,802 AND MSFragger 11,560/11,560/11,560. So ALL benchmark peptide counts
  are EXACT; n=1 is justified for both. CONSEQUENCE: the MSFragger k-curve non-monotonicity
  (11,402/11,560/11,482/11,751 across k=1..4) is NOT run-to-run noise -- it is a REAL (non-smooth)
  response to the different inputs (codex called this in review; now measured). The earlier "MSFragger
  is a noisy guardrail" framing was WRONG on the noise: MSFragger numbers are exact, just non-monotonic
  in k (its optimum k differs from Sage's). The ONLY stochastic uncertainty in the whole benchmark is
  the entrapment FDR estimate (small n), which the bootstrap CI below now quantifies.
* **Entrapment estimator was wrong TWO ways (fixed in `bench/entrapment.py`, codex #12):**
  (1) it counted peptides SHARED between human and Arabidopsis as entrapment -> inflated the fraction
  ~2.5x (dataset D base entrap 222 -> correctly 86 foreign-ONLY). (2) it used the PROTEIN ratio (0.800); the
  correct PEPTIDE-hypothesis ratio (in-silico tryptic, Sage rules) is 0.6805, correction 1/r = 1.469.
* **Corrected absolute FDR at nominal 1% peptide_q is ~1.0-1.45%, NOT the ~2.3% previously reported** --
  the tool is WELL-CALIBRATED. With 95% bootstrap CIs (n_entrap ~ 60-90 -> noisy):
  dataset D base 1.34% [1.06-1.65] / corr2 1.20% [0.95-1.46]; dataset A base 1.01% / corr2 1.08%; dataset B base 1.45%
  [1.15-1.77] / corr2 1.18% [0.91-1.47]. **The arm FDR CIs OVERLAP** -> the earlier "FDR improved
  2.33->2.13%" was an estimator+noise artifact. Correct statement: **corr_power adds +9-10% target
  peptides at an entrapment FDR statistically INDISTINGUISHABLE from base** (more peptides, same FDR).
* **Consequence for future levers:** the +1-3% regime (im_weight increments, cap sweeps, rescoring
  features) needs this CI machinery + MSFragger replicate sigma to be falsifiable; do not decide small
  levers on point estimates.

## Decision rule

* **Progress** = Sage peptides ≥ 10,802 (standard bench config = tool defaults + `-assembly:require_isotope_support`; re-based 2026-09-01, was 10,072 pre-gate/pre-corr_power) and PSMs/peptide not worse (guardrail retained).
* **Regression** = peptides fall ≥ 1%.
* **Noise floor** = 0.06% same-node (from four runs of an effectively identical config:
  8,411 / 8,408 / 8,406 on one node). Differences below ~0.2% mean nothing.
* **Screen on Sage** (4.2× faster, signs reproduce); confirm winners on MSFragger before any
  default change. Sage is unreliable on effect SIZE — measured compression up to 3× on one arm.

## Confirmed WINS post-freeze (gate baseline: Sage 9,989 / MSFragger 11,033 / entrap 9,304 @ 2.33%)

* **`assembly:corr_power` (Q2) + `im_weight_sigma` STACK — WIN on BOTH engines, entrapment-confirmed (2026-07-27).**
  Emitted fragment intensity *= co-elution^k. corr_power is EMIT-ONLY (excluded from the 500-cap key),
  so the corr_power arms preserve peak MEMBERSHIP (identical 924,255 spectra + same peak m/z set; pure
  intensity re-ranking). im_weight IS in the cap key, so with 73% of spectra at the 500-cap the STACK
  also REORDERS which fragments survive (its intended interferent-removal channel) -- same spectrum +
  peak COUNT, but not the same membership. So the stack's +310 over corr_power=2 has TWO channels
  (intensity reweight + cap membership), not one [kimi 2026-07-27]. the reference implementation on dataset D: Sage 11,517,
  MSFragger 13,932. MSFragger here is a NOISY guardrail (no replicate variance measured; screen on
  Sage+entrapment). The 2.33->2.13% entrap frac is WITHIN Poisson noise (n~85 entrapment) -- a guardrail
  (it passed: extra targets are not disproportionately entrapment), NOT an FDR win.
  | arm | Sage | MSFragger | entrap target @ actual-FDR |
  |---|---|---|---|
  | base (gate) | 9,989 | 11,033 | 9,304 @ 2.33% |
  | corr_power=1 | 10,622 | 11,402 | 9,894 @ 2.18% |
  | corr_power=2 | 10,802 | 11,560 | 10,235 @ 2.21% |
  | corr_power=3 | 10,901 | 11,482 | 10,288 @ 2.16% |
  | corr_power=4 | 10,875 | 11,751 | 10,262 @ 2.15% |
  | **corr_power=2 + im_weight 0.005** | **11,112 (96.5% of dt)** | **11,813 (84.8% of dt)** | **10,574 @ 2.13%** |
  corr_power alone PEAKS at k≈3 on Sage; the im_weight STACK is best on EVERY metric incl. the LOWEST
  entrap FDR (2.13%) -- the two are COMPLEMENTARY purity axes (RT co-elution vs IM), they ADD not erase.
  Sage 86.7% -> 96.5% of the reference implementation; MSFragger 79.2% -> 84.8%. Mechanism: a true fragment co-elutes with
  its precursor by construction (high c), so c^k down-weights mediocre-c INTERFERENTS. This OVERTURNED a
  UNANIMOUS 3-CLI review that predicted -2..-5% MSFragger (premise "faint precursor -> low-c fragment"
  was wrong). MSFragger is NOISY here (non-monotonic) -- screen on Sage+entrapment.

  **GENERALISATION TEST PASSED for corr_power (frozen k=2/sigma=0.005 on untouched dataset A, dataset B; 2026-07-27):**
  Sage nominal 1%: dataset A base 8,822 -> k2 9,735 (+10.4%) -> stack 10,428; dataset B base 8,311 -> k2 9,032 (+8.7%)
  -> stack 9,159. So **corr_power=2 robustly generalises: +8-10% over base on ALL THREE files, both
  engines (MSFragger dataset A +5.6%, dataset B +7.6%), at nominal 1% AND <=2% empirical FDR.** The im_weight STACK
  increment is FILE-VARIABLE (+310 dataset D, +693 dataset A, +127 dataset B over k2; flat on MSFragger) -- a Sage-leaning
  bonus, NOT a robust default increment (dataset B +127 < the pre-registered +200 threshold). The <=1% empirical
  FDR numbers are erratic (k2 dips on dataset A, stack dips on dataset B) = estimator noise at ~50 entrapment n (the
  <=2% and nominal metrics are clean). **DECISION (pending 3-CLI review of the result): ship corr_power=2
  as the new default; im_weight_sigma=0.005 as a documented optional stack. Caveat: dataset A/dataset B are the same
  study series as dataset D, so this de-risks the dataset D family, not diaPASEF broadly -- a different-cohort file is
  the stronger (unavailable-today) test.**

## What is already falsified on this file or dataset B

Do not re-propose without new evidence: `competitive`, rank-pruning (`rp_max` 2/4/8),
`assembly:apportion`, `trace:frame_aggregation_n` (3, 5), `charge:ambiguity_margin` (0.5/1.0/2.0
— and it was *unreachable* in count mode, so it is untested rather than falsified),
`consolidate:delta_rt` (0.7 → −12.7%, 3.0 → −24.0%), **`merge:rt_window` 1.4 → −4.4%**,
mass recalibration, a monotonic ion-mobility charge prior, AlphaDIA-style log-sum co-elution
scoring (three implementations), **`gate:variance_support` (Q1, union-support Pearson — 2026-07-27:
DROP on unanimous 3-CLI analysis; it measures support co-LOCALISATION not co-elution SHAPE and reads
faint fragments' censored shoulders as anti-correlation, e.g. precursor-8/frag-5/overlap-3 scores
G-Pearson +0.47 but union −0.50. G-Pearson already ≈ robust cosine, so [B0] was cosmetic).**

## Spectrum COLLAPSE (same precursor across cycles) — FALSIFIED for quality (2026-07-28)
Tested on the OPEN-search entrapment results (simulation on existing PSMs, no re-search): group rank-1
PSMs by precursor coordinate ACROSS RT (charge + expmass 10 ppm + IM 0.01, single-link chained over a
10 s RT gap), keep one (or top-N) per group, RECOMPUTE the corrected entrapment FDR from scores.

| rule | spx kept | reduction | targets @1% FDR | dt targets |
|---|---|---|---|---|
| full            | 558,255 | 1.00x | **8,273** | 9,743 |
| apex (ID-agnostic) | 86,842 | 6.43x | 3,342 (-60%) | 4,491 (-54%) |
| best (ORACLE)   | 86,842 | 6.43x | 4,390 (-47%) | 5,089 |
| top3 / top5 / top10 (collapse + chimeric report_psms) | 132k/155k/190k | 4.2x/3.6x/2.9x | 5,448/5,767/6,603 | 6,068/6,672/7,422 |

**The multiple-testing-burden hypothesis is REFUTED.** Every arm's FDR CI sits at ~0.6-1.3%: collapse
removes true AND entrapment PSMs PROPORTIONALLY, so the FDR does not improve -- peptides are simply
lost. Emitting fewer spectra does NOT buy open-search sensitivity.
**It is INTRINSIC, not a SpeXtract defect:** the reference implementation collapses the same way (-24% at top10, -54% at
apex) and BOTH tools reduce to the SAME ~87k precursor features (spx 86,842 / dt 87,112) -- SpeXtract
just emits 1.33x more spectra per feature (6.43 vs 4.84). The across-cycle multiplicity is INDEPENDENT
EVIDENCE (each cycle = a different noise realisation + co-isolation mixture, another shot at a different
co-eluting peptide), not redundancy. Consistent with the long-standing "across-cycle spread 2.66 is
CORRECT (MS1 FWHM = 2.61 cycles)" finding.
**Residual option, untested:** TRUE merging (union of peaks -> better S/N) could in principle ADD IDs
rather than select among existing ones; this simulation only bounds the SELECT-among-existing case. But
it must beat a -20%..-60% deficit and it blends interferents across RT, so it is not promising.
**Practical use:** collapse is a SPEED/recall knob, not a quality lever -- top10 gives 2.94x fewer
spectra (much cheaper search) for -20% peptides, if total pipeline time ever matters more than recall.

## Monoisotope -1.003 Da open-search artefact: STOPPING guard FALSIFIED, scoring fix under test (2026-07-28)
The open-search re-baseline found our largest delta-mass artefact is a **-1.003 Da population** (top bin
7,616 PSMs; 8.20% isotope-shifted vs the reference implementation 3.46%): the `charge:scoring=count` walk defines the
monoisotope as the FURTHEST-LEFT peak it can reach, and `findPartner` matches m/z + RT + IM but **never
intensity**, so it latches onto noise or a co-eluting species one isotope below the true mono. Closed
search HIDES this (Sage enumerates `isotope_errors [-1,2]` and corrects it for free); open/blind search
turns it into a phantom -1 Da modification.

**Attempt 1 — `charge:mono_averagine_guard` (a STOPPING rule): FALSIFIED.** Require the leftward
candidate to reach the averagine minimum 1/lambda times the previous peak. At slack 0.5 on dataset D:
| metric | control | guard 0.5 |
|---|---|---|
| -1 Da bin | 7,616 | **2,734 (-64%, the target WORKED)** |
| +1 Da bin | 4,697 | **23,994 (+411%)** |
| total isotope-shift | 8.20% | **16.43% (WORSE)** |
| open targets @ FDR | 8,887 @1.29% | 8,939 @1.28% (flat, CIs overlap) |
| closed Sage | 10,802 | 11,697 (+8.3%) |
| emission | 924,255 | 1,005,353 (+8.8%) |
**A stopping rule can only trade over-reach for under-reach**, and worse: a walk that halts early leaves
the M+1 peak unconsumed by `used[]`, so it later seeds its OWN precursor one isotope too heavy -- that is
BOTH the +1 explosion and the +8.8% emission (one cause, two symptoms). Open search did not improve.
NOTE the closed-search 11,697 (> the reference implementation 11,517) is NOT a quality win to bank: peptides/1k spectra is
11.7 -> 11.6, i.e. the gain tracks the extra emission -- exactly the over-emission-rewards-closed-search
corruption the reviewers warned about. Do not ship the guard for it.

**Attempt 2 — `charge:mono_averagine_select` (a SCORING decision): under test.** Walk the full run, then
choose WHICH peak is the mono by averagine cosine of the envelope starting there. Corrects BOTH
directions, and every found peak is still consumed, so de-isotoping and emission stay unchanged (only the
reported precursor mass moves) -- which also removes the emission confound from the experiment.

### Attempt 2 result + the finding that STOPS this line of work (2026-07-29)
`charge:mono_averagine_select` raw cosine repeated the envelope scorer's documented short-vector trap
(`len` shrinks as the candidate mono moves right -> short vectors score high -> too-HEAVY mono). Adding
the same evidence weight the envelope scorer uses, (len-1)/3 capped at 1, largely fixed that:
| metric | control | guard 0.5 | msel raw | **msel + evidence weight** |
|---|---|---|---|---|
| emission | 924,255 | 1,005,353 | 924,716 | **924,552** |
| -1 bin | 7,616 | 2,734 | 2,314 | **4,015 (-47%)** |
| +1 bin | 4,697 | 23,994 | 21,591 | **12,188** |
| unmod | 31.2% | 24.3% | 26.4% | **32.0% (best)** |
| isotope-shift | **8.20%** | 16.43% | 15.24% | **9.48%** |
| open targets @FDR | 8,887 @1.29% | 8,939 @1.28% | 8,916 @1.33% | 8,932 @1.33% |
| closed Sage | 10,802 | 11,697 | 11,458 | 11,279 |
Pre-registered criterion was total isotope-shift BELOW 8.20% with both +-1 bins down: **NOT MET (9.48%)**.
Stopped rather than running a 4th variant. Both flags stay default-OFF, documented, for reproducibility.

**THE FINDING THAT MATTERS — the -1 artefact is NOT the open-search recall lever (my hypothesis was WRONG).**
Open-search targets are FLAT across all four arms (8,887 / 8,939 / 8,916 / 8,932; every CI overlaps)
even though the -1 population was cut 47-70%. Reason, obvious in hindsight: **a -1-shifted precursor is
still IDENTIFIED** -- it appears in the delta-mass histogram AS an identified PSM at -1.003 Da. Open
search does not lose the peptide, it mislabels it with a phantom modification. So the monoisotope
artefact is a **PTM-INTERPRETATION problem (real, and it matters for blind-PTM correctness, our actual
differentiator), NOT a recall problem.** It cannot explain the 8,887 vs 9,989 open-search gap; do not
spend more effort here expecting recall. The open-search recall gap needs a different explanation --
candidates not yet tested: our 924k vs the reference implementation's 700k emission competing in the huge open search
space, fragment quality on the faint tail, or precursor mass ACCURACY (ppm) rather than isotope offset.

## STEP 01 (2026-09-01): the metric re-denominated + denominators reconciled

### A. Delta-mass-BIN-level open-search scoring (bench/open_ptm_score.py v2, corrected estimator)
Unit = (peptide, delta-bin); classes scored separately: unmod (|d|<=8 mDa), nearzero (8-900 mDa =
mass-error tail, NOT a mod), artifact (iso lattice +-k*1.00335, k=1..9), known (curated Unimod list),
other (off-lattice, per-BIN FDR walks for top bins -- pooled class FDR let junk bins ride; the spx
-485..-499 Da window-edge wall at ~36% raw entrap collapses to ~30-80 @1% under per-bin walks).
| class @1% corrected FDR | SpeXtract | the reference implementation | ratio |
|---|---|---|---|
| unmodified | 4,885 | 7,015 | 70% |
| known-PTM (the raison-d'etre number) | 2,692 | 3,092 | **87%** |
| isotope-artifact (k<=9) | 4,073 | 3,611 | 1.13x (v1's "2.4x artifact" was a k<=3-grid artifact of its own: dt's big envelope ladder is ON-lattice at k>=4; ours skews to iso-1) |
| **nearzero mass-error (8-900 mDa)** | **7,099** | **2,814** | **2.5x -- NEW finding** |
| other (off-lattice) | 44,143 | 78,244 | 56% (interpret with care: correlated bins; contains e.g. -0.984 amidation-type masses off the iso lattice) |
**Findings:** (1) In the honest currency the reference implementation leads every legitimate class; the old
"89% targets@FDR" UNDERSTATED the unmod gap (70%) and roughly matched known-PTM (87%).
(2) NEW: spx's near-zero mass-error tail is 2.5x dt's, and spx unmod+nearzero = 11,984 vs dt 9,829 --
much of our "missing unmodified" population is IDENTIFIED but MASS-DEGRADED (precursor 8-900 mDa off),
i.e. an interpretation-quality defect (like the -1 Da class), pointing at PRECURSOR MASS ACCURACY.
NOTE: "mass recalibration" sits in the falsified list, but it was falsified ON THE CORRUPT CLOSED
METRIC -- exactly the class of falsification the 2026-09-01 review said may need re-audit. Flagged for
re-evaluation in THIS currency; not re-litigated yet. (3) The phantom-mod gap (artifact class) is
+13%, not 2.4x -- the earlier headline overstated it by using a too-short isotope grid.

### B. Denominators reconciled (the "no ratio is citable" blocker)
- **CANONICAL the reference implementation references from 2026-09-01 on: Sage 11,517** (dt_s30/results.sage.tsv,
  nominal peptide_q<=0.01 -- re-verified today; TSV archived to ceph) and **MSFragger 13,932**
  (msf_dt/fixed.tsv via score_msf_td target-decoy; archived).
- **11,552 / 13,014 are DEPRECATED**: products of the 2026-07-22 joint common-FDR procedure whose
  primary outputs (joint_results.json et al.) exist in NO file on the cluster -- doc-only numbers.
  Not to be used in ratios. The frozen table's 86.4%/80.5% ratios are legacy.
- **Decision-rule re-base:** the frozen threshold 10,072 predates the isotope-support gate AND
  corr_power (superseded twice). New reference for progress/regression: **Sage 10,802** = current
  tool defaults (corr_power=2) + the standard bench config (which passes -assembly:require_isotope_support
  explicitly; that flag remains OFF in the tool). Noise floor and +-1% rules unchanged.

### STEP 01 CLOSURE (v3 scorer, post kimi+codex review, 2026-09-01) — SUPERSEDES the v2 table above
The v2 class table is NOT citable (review: Da-constant tolerance made every ratio a mass-accuracy
measurement biased against spx; per-bin @1% on zero-entrapment slivers vacuous; ontology gaps).
v3 (bench/open_ptm_score.py @ e9f5b64b): nearest-candidate assignment over {0, iso-lattice k<=15,
28 KNOWN masses}, ppm-scaled tol max(8 mDa, 12 ppm), conservative (e+1)/r tie-grouped walks with
bootstrap CIs, min-evidence rule (per-name @1% only when accepted e>=10), provenance headers.

**Citable v3 table (dataset D open search, @1% corrected class-FDR, CIs in ptm_score_v3.out):**
| class | spx | dt | reading |
|---|---|---|---|
| unmod | 6,907 | 7,890 | 87.5% at hypothesis level; PEPTIDE level ~parity (96-106% by definition: accepted-union 7,744 vs 8,079; best-hypothesis 6,547 vs 6,164) |
| knownPTM | 4,878 | 4,651 | "105%" is CONTAMINATED: spx's top-3 knownPTM bins are all iso-AMBIGUOUS (Didehydro/Amidation/Deamid = 34% of spx knownPTM units vs dt 14%) = our lattice spill relabeled. Ambiguity-excluded units: spx ~78%. Honest: BROADLY COMPARABLE SCALE, per-mod 57-131% (Ox 93%, Acetyl 131%, Phospho 57%), NO single mod FDR-certifiable (all per-name e<10) |
| isotope lattice (k<=15) | 7,737 | 8,826 | dt carries MORE lattice hypotheses (its +k up-ladder); spx signature = iso-1 (1,790 @1%, e=11, certified). NEITHER tool is clean |
| nearzero (mass-error) | 5,415 | 1,432 | 3.8x — THE spx-specific defect, mechanism nailed (below) |
| other | 40,160 | 72,747 | DESCRIPTIVE ONLY: dominated in BOTH tools by off-lattice integer-Da bins (dt's ladder to +21); not PTM discovery |

**The nearzero mechanism is PROVEN (three independent legs, all on existing data):**
(1) 70.1% of spx nearzero peptides ARE dt's unmod peptides; (2) 89.4% same-scan concordance with our
OWN closed search (same scan -> same peptide; 0.7% different -> chimeric explanation dead);
(3) deltas 93.8% negative, median -16.2 mDa, slope -18.3 ppm vs calcmass (dt control ~symmetric).
=> SpeXtract's reported precursor masses carry a SYSTEMATIC ~-18 ppm CALIBRATION BIAS introduced by
our pipeline (same raw data as dt). Fixable (per-run linear recalibration of reported mono m/z);
distinct from the old falsified "mass recalibration" lever (different target, corrupt metric).
Scheduled as its own arm AFTER step 02 (both reviewers: do not confound the pre-registered arm).

**MSFragger denominator reconciled (codex #14):** 13,932/11,560 are PEPTIDOFORM-level (pair them);
sequence-level recount from archived TSVs: dt 12,752 vs spx 11,055 = 86.7%; legacy 13,014 was
sequence-level and is within 2% of 12,752 — unit mismatch closed, ratios now like-for-like.

**Step-01 bottom line:** the old "89% on open search" dissolves. At peptide level we are at ~parity
on unmodified; known-PTM scale is comparable (ambiguity-limited both ways); BOTH tools' open-search
delta dimension is dominated by non-PTM structure (lattice + integer-Da junk = 5-9x the legitimate
classes) — a tool-CLASS finding worth publishing; and SpeXtract's one distinctive open-search defect
is the -18 ppm mass bias (fixable), vs the reference implementation's up-ladder envelope artifacts. Step 02 (emission
arm) runs AS REGISTERED with a nearzero-stratified secondary readout, scored under this v3 metric.

## THE -8 ppm MASS BIAS: mechanism found + fixed (2026-09-01, follows the step-01 nearzero finding)
**Every SpeXtract-reported m/z to date carries a systematic, m/z-dependent -5..-11 ppm error**
(median -8.3 ppm on identified precursors; dt control -1.1 ppm on the same raw data).
**Mechanism (proven by direct comparison against Bruker's timsdata library via ctypes):** the loader
(BrukerTimsFile -> opentims++) converts TOF->m/z with `OpenSourceTof2MzConverter (linear-in-sqrt)` --
a TWO-POINT chord fit to just (MzAcqRangeLower=100, MzAcqRangeUpper=1700, DigitizerNumSamples): error
~0 at the anchors, bowing to -10.8 ppm mid-range, frame-independent -- reproducing every observed
feature (flat in RT and mass-in-ppm; the "IM trend" is the m/z-dependence in disguise; temperature was
a red herring, dC2=0 and dT1 tiny). The loader's SDK branch (OPENMS_BRUKER_SDK_PATH) upgraded ONLY the
ion-mobility converter -- the m/z converter was never swapped (and the env var was never set anyway).
**Fix: patches/openms-brukertims-sdk-mz.patch** -- 3 lines: the SDK branch now also installs
`BrukerTof2MzConverterFactory` (vendor per-frame calibration via tims_index_to_mz). Requires
OPENMS_BRUKER_SDK_PATH pointing at a local Bruker timsdata .so (present on the cluster from fp-tools;
used locally, never redistributed; without it the tool falls back to the open model as before).
The open-model inaccuracy is itself reportable upstream (opentims; also explains why mzdata disabled
its own m/z model port).
**Why mzPeak was NOT the answer (user directive "pull latest mzPeak", researched 2026-09-01):** spec
is draft 0.9 / prototype v0.1.0 / breaking layout change 2026-07-14 / no releases; the current
TDF->mzpeak converter stores RAW TOF indices and its own m/z path uses the same two-point model
(mzdata implemented the true vendor model 0.66.3 but disabled it 0.66.5 as "not consistently better");
the OpenMS in-tree MzPeakFile is unmerged/dormant and keyed to a superseded spelling. mzPeak re-poses
the calibration problem; it does not solve it.
**Impact on prior results:** [CORRECTED by review: the original "closed counts unaffected" claim
here was falsified by sdkcal itself, +12.3% -- candidate-window inclusion was never the only channel];
step-01 v3 open-search numbers stand AS MEASURED but the nearzero class (5,415 vs 1,432) was this bug.
Decisive arm `sdkcal` running: vendor-calibrated re-extraction -> expect ppm median -> ~-1, nearzero
to collapse toward dt's level, and possibly small closed gains (tighter effective tolerance).

### sdkcal CONFIRMED (2026-09-01): vendor m/z calibration is the new benchmark configuration
| metric | pre-fix (corr2) | **sdkcal (corr2 + vendor m/z)** | the reference implementation |
|---|---|---|---|
| ppm median (open IDs) | -8.29 [-12.5,-4.0] | **+1.74 [-2.4,+6.1]** | -1.10 [-3.8,+1.2] |
| closed Sage | 10,802 | **12,128 (+12.3%)** | 11,517 (scoped claim only -- see review verdict below) |
| closed Sage entrapment | 10,381 @1.20% | **11,435 @1.26% [1.00-1.51]** | -- (FDR unchanged -> gain is real) |
| closed MSFragger | 11,560 | 11,460 (-0.9%, FLAT) | 13,932 (82.3%) -> Sage-specific gain |
| emission | 924,255 | 922,132 (~unchanged) | 700,434 |
| open unmod @1% | 6,907 (87.5%) | **7,566 (95.9%)** | 7,890 |
| open nearzero @1% | 5,415 (3.8x dt) | **2,769 (1.9x)** | 1,432 |
| open knownPTM @1% | 4,878 | 6,058 (total; ambiguity-excluded units ~85% of dt) | 4,651 |
**Decision: OPENMS_BRUKER_SDK_PATH + patches/openms-brukertims-sdk-mz.patch = REQUIRED benchmark
config from now on; closed-Sage reference re-based 10,802 -> 12,128.** Caveats carried into review:
(1) the +12% is Sage-specific (MSFragger flat -- its gap is search/detection-side, consistent with
prior reviews); (2) EVERY pre-2026-09-01 number was measured at -8 ppm, INCLUDING the corr_power
validation and all falsifications -- corr_power=0 ablation under vendor-cal launched (sdkcal_cp0);
dataset A/dataset B holdout re-validation owed; falsified-lever re-audit only where a lever plausibly interacted
with mass accuracy; (3) knownPTM ambiguous fraction unchanged (34.7% vs dt 14%) but the bins are now
physically separable at +-4 ppm -- the 1,610-unit Amidation bin (1.4% raw entrap) deserves real
chemical scrutiny rather than automatic dismissal as spill.

### Calibration-step adversarial review verdict (kimi + codex, 2026-09-01) — 12,128 is PROVISIONAL
Convergent verdict: the chord mechanism is the proven DOMINANT cause of the old raw-axis bias, but the
step over-claimed in four ways, now corrected:
1. **Attribution not isolated:** setting the env var flipped BOTH converters (m/z AND ion mobility) --
   the SDK IM swap engaged for the first time ever alongside the m/z patch. Part of the +12.3% could be
   IM-side; ALL prior IM-dependent results (incl. im_weight's validation) silently ran on the open IM
   model. IM-isolation arm required (m/z vendor + IM open).
2. **12,128 = provisional dataset D/vendor-SDK regression oracle, NOT a validated reference** until: cp0
   ablation (running), the IM-isolation arm, and dataset A/dataset B re-validation under vendor-cal.
3. **"Leads the reference implementation" is not citable.** Maximum defensible wording (codex): "On dataset D, using Sage and
   the vendor-SDK calibration configuration, SpeXtract yielded 5.3% more closed-search peptides; under
   MSFragger it yielded 17.7% fewer." Per-1k-spectra efficiency: spx 13.2 vs dt 16.4 (dt +25%) -- the
   emission column stays in any public row. MSFragger flatness is UNINFORMATIVE (its default
   calibrate_mass auto-corrects input bias), not evidence against the fix.
4. **Residual defect remains:** sdkcal centers +1.74 ppm vs dt -1.10 (2.8 ppm apart, quartiles
   symmetric -- earlier "right-skewed" claim wrong) and nearzero is halved, not closed (1.9x dt).
   Candidates: our centroiding/mono reporting point. Paired same-feature spx-vs-dt reported-m/z
   comparison is the discriminator. The -18 ppm (step-01 regression slope on the truncated nearzero
   subset) vs -8.3 ppm (median, all IDs) are different statistics of the same defect -- reconciled here.
Hardening obligations recorded: benchmark must FAIL CLOSED on a missing/wrong SDK path (silent
fallback = known-wrong masses); opentims' BrukerTof2MzConverter DISCARDS tims_index_to_mz's return
code (silent incomplete-buffer risk -- fix when the patch is vendored properly); SDK identity pinned
in evidence/SHA256SUMS (sha a9708613..., PRERELEASE TDF-SDK 2.21.0 -- Bruker notes a tims_index_to_mz
implementation change in this line); patch base = OpenMS 3.6.0-pre-exported-20260717 tarball; runtime
headlines must be RE-MEASURED under the SDK path before publication. Determinism check (threads-8 vs
threads-100 semantic/byte hash) running. Two-tier calibration product conflict stands: license-blocked
users (the stated market) get -8 ppm masses from the open fallback -- an open true-polynomial port or
loud two-tier disclosure is a v0.2.0 gate.

### cp0 ablation result (2026-09-01): corr_power SURVIVES the calibration fix — provisionality condition 1/3 cleared
corr_power=0 + vendor-cal = 11,017 vs corr_power=2 + vendor-cal = 12,128. The 2x2 (closed Sage):
old-cal cp0 9,989 / cp2 10,802 (+8.1%); vendor-cal cp0 11,017 / cp2 12,128 (+10.1%). The levers are
independent and additive (mild positive synergy); corr_power's +8-10% validation was NOT a
calibration artifact. Note the cp0 arm also ran vendor IM, so the m/z-vs-IM attribution of the
calibration gain itself still awaits the IM-isolation arm (condition 2/3); dataset A/dataset B = condition 3/3.

### Table-model calibration: 3-FILE RESULT (open path, no vendor .so; 2026-09-01)
| file | old chord | **TDF-table model** | gain | vendor-SDK | ppm (table) |
|---|---|---|---|---|---|
| dataset D | 10,802 | **11,976** | +10.9% | 12,128 | +1.49 [-2.4,+5.4] |
| dataset A |  9,735 | **10,333** | +6.1%  | -- | +3.19 [-0.8,+7.1] |
| dataset B |  9,032 | **9,891**  | +9.5%  | -- | +2.29 [-1.7,+6.3] |
**The calibration gain GENERALISES: +6-11% closed Sage on all three files with NO vendor library.**
Provisionality conditions on the re-based reference: (1) cp0 ablation CLEARED; (2) IM-isolation
ANSWERED FOR FREE -- the table arm is vendor-exact m/z + OPEN (rational) IM, the sdkcal arm was
vendor both: dataset D 11,976 vs 12,128 => ~93% of the +12.3% is m/z-side, ~1.3% (152 peptides) is the
vendor-IM contribution; (3) dataset A/dataset B re-validated above. **New open-path reference: dataset D = 11,976**
(the 12,128 SDK number stays as the vendor-oracle upper bound, not the shipping config).
Residual +1.5..+3.2 ppm persists across all files and BOTH calibration paths -> consistent with the
centroiding/mono-reporting hypothesis (Option I), not the m/z axis. the reference implementation per-file ppm + dataset A/dataset B
head-to-head running to test whether part of it is instrument/per-file rather than ours.

### FIRST 3-FILE HEAD-TO-HEAD vs the reference implementation (Sage closed, corrected calibration, 2026-09-01)
the reference implementation had ONLY ever been measured on dataset D; dataset A/dataset B references created here under identical config.
| file | SpeXtract (table model) | the reference implementation | ratio | our ppm | dt ppm |
|---|---|---|---|---|---|
| dataset D | **11,976** | 11,517 | **104.0%** | +1.49 | -1.37 |
| dataset A | **10,333** | 10,242 | **100.9%** | +3.19 | +0.27 |
| dataset B | **9,891**  |  8,948 | **110.5%** | +2.29 | -0.60 |
**SpeXtract reaches 100.9-110.5% of the reference implementation on closed Sage across the three files** -- no vendor
library, open BSD path. STATISTICALLY (n=3 paired): mean 105.1%, 95% CI [93.0, 117.3] = consistent
with PARITY, not an advantage; the earlier phrasing "matches or exceeds on ALL THREE files" was
over-claiming on dataset A's +91. CORRECTION (review F4): the "ours 1.26% [1.00-1.51]" entrapment figure
quoted here came from the **sdkcal** arm, a different converter -- **the shipping table-model arm's
entrapment FDR has NOT been measured on any file**; that is an owed measurement, not a passed check.
CAVEATS THAT STAND: (a) emission is still ~1.32x dt, so per-spectrum efficiency still favours dt --
publish the emission column with any count row; (b) MSFragger arm pending (dt led 82% there; that
engine's gap is search/detection-side); (c) single cohort/instrument -- public-PXD replication is
still the external-claim gate; (d) our ppm residual (+1.5..+3.2) is 2-4 ppm ABOVE the reference implementation's
(-1.4..+0.3) on every file: dt's per-file spread (-1.37/+0.27/-0.60) is small and centred, so the
residual is OURS (centroiding/mono-reporting, Option I), NOT a per-file instrument property.

### C2 PAIRED m/z bias vs the reference implementation (2026-09-02, 15:14) -- the residual is ours, in BOTH directions
Same (peptide, charge) identified by both tools at 1% (best PSM per tool), so selection bias is excluded:
| file | shared pairs | precursor ppm ours / dt / **paired ours−dt** | fragment ppm **paired ours−dt** | by charge (paired, precursor) | by m/z tercile |
|---|---|---|---|---|---|
| dataset A (tbl_s08 vs dt3_s08) | 10,091 | +5.17 / +2.13 / **+2.63** [p10 −2.6, p90 +10.3] | **−1.99** [−4.5, +0.5] | z2 +3.12, z3 +1.77, z4 +1.96 | low +1.87, mid +3.35, high +2.62 |
| dataset B (tbl_s23 vs dt3_s23) | 9,094 | +3.78 / +1.97 / **+1.48** [−3.4, +8.4] | **−2.08** [−4.4, +0.3] | z2 +1.96, z3 +0.53, z4 +0.29 | low +0.42, mid +1.88, high +1.87 |
| **dataset D** (d2_P1 = HEAD table arm vs dt_s30_closed, fresh the reference implementation Sage run = 11,517 peptides, reproducing the matrix reference) | 11,191 | +3.58 / +2.05 / **+1.12** [−3.7, +7.7] | **−2.41** [−4.7, +0.1] | z2 +1.59, z3 +0.31, z4 +0.54 | low +0.05, mid +1.62, high +1.49 |
Reading: on the same raw file and the same exact m/z axis, our REPORTED precursor m/z sits +1.5..+2.6 ppm
above the reference implementation's and our REPORTED fragment m/z sits ~2 ppm below -- two different centroid/reporting
conventions (precursor: MS1 trace m/z estimator, larger at z=2 and in the mid/high m/z terciles; fragment:
MS2 trace m/z estimator), not a calibration offset (which would move both the same way). Both feed Sage's
discriminant (precursor_ppm, average_ppm) and the 20 ppm fragment tolerance. Next: locate the two estimators
(`Trace.mz` for MS1 and MS2 traces: intensity-weighted mean over the trace vs apex) and A/B an apex /
top-k-weighted estimator behind the dataset D set gate, reporting ppm beside the count. dataset D pair running
Three-file pattern: the precursor bias is z=2-dominated (z3/z4 near zero on dataset D/dataset B) and absent in the lowest m/z tercile -- an isotope/mono-envelope effect (our reported mono for z=2 mid/high-m/z precursors sits a fraction of a ppm-scaled isotope offset high?) rather than a uniform centroid shift; the FRAGMENT bias is uniform −2.0..−2.4 ppm on all three files and is the one that touches the 20 ppm search tolerance and Sage's average_ppm feature (−3 ppm uniform shift cost −165 peptides on 09-01). Estimator A/B (apex / median vs weighted mean) queued behind the E5 gate.

### C3 MEASURED (2026-09-02, 16:41): entrapment FDR of the SHIPPING table-model arm
`entrap_apply.py` (peptide-hypothesis ratio 0.6805, corrected estimator) on a Sage entrapment search of the HEAD
output (d2_P1, table-model calibration, parallel loader -- byte-identical to every run today):
| arm | targets @1% | entrapment | raw % | **FDR %** | 95% CI |
|---|---|---|---|---|---|
| **table model (shipping, HEAD)** | 11,489 | 107 | 0.92 | **1.37** | 1.10-1.64 |
| vendor-SDK calibration (09-01 arm) | 11,435 | 98 | 0.85 | 1.26 | 1.00-1.51 |
| **dataset A table model (shipping)** | 9,736 | 63 | 0.64 | **0.95** | 0.72-1.18 |
| dataset A the reference implementation | 9,953 | 84 | 0.84 | 1.24 | 0.99-1.54 |
| **dataset B table model (shipping)** | 9,462 | 85 | 0.89 | **1.32** | 1.04-1.59 |
| dataset B the reference implementation | 8,581 | 71 | 0.82 | 1.22 | 0.96-1.53 |
The owed measurement (review F4) is now made: the shipping arm's true FDR at nominal 1% is ~1.4%, inside the
same band as the vendor-calibrated arm (CIs overlap), so the open-path calibration gain is not bought with
FDR. On all three files the shipping arm's true FDR at nominal 1% is 0.95-1.37%, inside the reference implementation's band on the two files where it exists (1.22-1.24%): the 3-file closed-search parity claim now carries its FDR control.

### C7: `assembly:apportion` FALSIFIED CLEANLY (2026-09-02, 16:46) -- the first A/B with corr_power actually applied
Every earlier apportion/rp_max A/B compared corr_power=2 (share-all) against corr_power=0 (the variant) because those
branches bypassed the emit weights (E5, fixed 15:05). Re-run with the weights applied: apportion=1.0 -> Sage 11,742
vs 12,217 share-all (**-3.9%**, only-share-all 1,234 / only-apportion 759). Same verdict as July, now clean:
cross-precursor intensity apportionment loses peptides. Consequence: BACKLOG item 1 (NNLS unmixing onto the MS1
basis) is cut -- its prerequisite ("apportion shows a signal") failed under the correct measurement. Side note: the
apportion path ran the window loop at 6.1x (3,658 s) -- an hour-long run -- irrelevant now that it is dead.

### C2 estimator A/B (2026-09-02, 17:22): APEX m/z for traces is +2.6% peptides (same binary, env switch)
`SPEXTRACT_MZ_ESTIMATOR` selects the reported m/z of every trace in `toTrace()`; default = OpenMS centroid
(intensity-weighted mean over the trace's members).
| estimator | spectra | Sage @1% | vs share-all default (12,217) | paired precursor ppm vs dt | paired fragment ppm vs dt |
|---|---|---|---|---|---|
| mean (default) | 922,902 | 12,217 | — | +1.12 | −2.41 |
| **apex** (max-intensity member's m/z) | 927,813 | **12,537 (+2.6%)** | only-default 757 / only-apex 1,077 | **+0.78** | −2.41 |
| median | 920,021 | 12,099 (−1.0%) | 756 / 638 | +1.18 | −2.36 |
Same binary -> the 0.06% same-binary floor applies; +2.6% is real on Sage. The precursor bias moves a third of
the way to the reference implementation; the FRAGMENT bias does not move at all -> the −2.4 ppm fragment offset is upstream of the
trace estimator (the IM-cluster pick's centroid, or the reference implementation reporting the apex RAW peak), a separate A/B.
Default change is gated on the standing rule: **MSFragger on d7_apex before apex becomes the default** (queued).

### APEX estimator passes the both-engines gate -> NEW DEFAULT (2026-09-02, 18:42)
| dataset D | Sage @1% | MSFragger @1% (raw-hyperscore walk) | precursor ppm paired vs dt |
|---|---|---|---|
| table model, mean estimator (until today) | 12,217 | 11,463 | +1.12 |
| **table model, APEX estimator** | **12,537 (+2.6%)** | **11,927 (+4.1%)** | **+0.78** |
| the reference implementation | 11,517 | 13,932 | — |
MSFragger ratio vs the reference implementation 82.3% -> 85.6%. Same binary (env switch), so the same-binary floor applies; both
engines move the same way; ppm beside the count improved. Default changed in code (`SPEXTRACT_MZ_ESTIMATOR=mean`
restores the OpenMS centroid). Owed next: dataset A/dataset B apex arms (Sage) and entrapment on the apex arm (queued).

### C9 FALSIFIED (2026-09-02, 19:05): dropping the precursor's M+1..M+3 from fragment lists LOSES peptides
`SPEXTRACT_DROP_PREC_ISO=1` on the apex-default binary: isotope contamination of the emitted lists falls from
28.4% / 24.8% (M+1 / M+2 within 10 ppm) to 0.2% / 0.2%, and Sage drops **12,537 -> 12,249 (-2.3%**, only-apex
1,061 / only-dropiso 773; same binary, env switch). The unfragmented precursor isotopes in a pseudo-MS/MS spectrum are
useful to the engine (or their removal lets weaker peaks refill the 500-cap). Keep them; switch stays off.

### Pick-level m/z modes FALSIFIED (2026-09-02, 20:00) -- and the fragment offset is not in our estimator chain
`pickIMCluster:mz_mode` (OpenMS PeakPickerIM patch; default `weighted` = unchanged, verified byte-identical to the
apex-default run):
| pick m/z | spectra | Sage @1% | vs apex default (12,537) | paired precursor / fragment ppm vs dt |
|---|---|---|---|---|
| weighted (default) | 927,813 | 12,537 | identical | +0.78 / −2.41 |
| seed (most intense point) | 896,143 | 12,196 | −2.7% | +0.73 / −2.38 |
| top3 (weighted mean of 3) | 926,359 | 12,325 | −1.7% | +0.70 / −2.42 |
The −2.4 ppm fragment offset vs the reference implementation survives BOTH the trace estimator (apex/median/mean) and the pick centroid
(weighted/seed/top3) -- it is not produced by anything in our reporting chain after the TOF->m/z conversion, which is
the same exact model on both sides (2.5e-5 ppm). What remains: the reference implementation's own fragment m/z reporting/recalibration.
That is not a defect to chase on our side; the precursor axis is where our reporting differs (apex fixed a third).

### APEX estimator: 3-FILE RESULT + entrapment (2026-09-02, 20:57) -- the new open-path reference
| file | mean estimator (09-01 reference) | **apex (new default)** | gain | the reference implementation (Sage) | ratio | paired precursor ppm vs dt (mean → apex) | fragment |
|---|---|---|---|---|---|---|---|
| dataset D | 12,217 | **12,537** | +2.6% | 11,517 | **108.9%** | +1.12 → +0.78 | −2.41 |
| dataset A | 10,333 | **10,789** | +4.4% | 10,242 | **105.3%** | +2.63 → +1.47 | −1.95 |
| dataset B | 9,891 | **10,156** | +2.7% | 8,948 | **113.5%** | +1.48 → +1.17 | −2.13 |
Mean ratio 109.2% (n=3 paired; 09-01 was 105.1%). Entrapment of the dataset D apex arm: **1.28% [1.03-1.53]** at nominal 1%
(table/mean arm 1.37% [1.10-1.64]; the reference implementation 1.22-1.24% on dataset A/dataset B) -- the gain is not bought with FDR. MSFragger on
dataset D: 11,927 vs 11,463 (+4.1%; 85.6% of the reference implementation's 13,932). Runtime note: dataset A (1.23 M spectra, the largest file)
26:51 wall / 156 GB peak; dataset B 15:28 / 109 GB. The apex arm outputs are the new benchmark reference
(`d7_apex`, `d13_apex_s08`, `d13_apex_s23`); the mean-estimator numbers stay as the 09-01 row.

### C1(b) rescoring A/B, first two attempts (2026-09-02, 21:30)
FragPipe 24.0 headless (Basic-Search workflow, 25 ppm, our human_decoy FASTA) on the dataset D apex arm and on the reference implementation.
MSBooster (DIA-NN predictions) aborts on OUR pepXML both with and without spectrum prediction: "Prediction missing
in file for REM[15.9949]DQTM[15.9949]AANAQK|3" (a PSM whose peptide DIA-NN's predictor does not emit). the reference implementation's
arm completed: **18,670 peptides / 37,086 PSMs at 1% (Percolator + Philosopher sequential/picked)** vs 13,932
peptidoforms in the raw-hyperscore walk -- rescoring is worth +34% on the reference implementation's spectra, which is exactly why the
comparison must be made on both. Third attempt running: MSBooster off, Percolator on MSFragger features only, both
tools (like-for-like).

### C1(b) RESCORING A/B, like-for-like (2026-09-02, 21:46): rescoring does NOT close the MSFragger gap
FragPipe 24.0 headless, Basic-Search (MSFragger 25 ppm, our human_decoy FASTA) -> Percolator -> Philosopher
(sequential, picked, 1% peptide), MSBooster OFF for both (it aborts on our pepXML, see above):
| dataset D arm | peptides @1% | PSMs @1% | PSMs / peptide | ratio spx/dt |
|---|---|---|---|---|
| SpeXtract apex | 13,211 | 113,685 | 8.6 | — |
| the reference implementation | 15,947 | 31,700 | 2.0 | **82.8%** |
| (raw-hyperscore walk, for reference) | 11,927 vs 13,932 | | | 85.6% |
| (the reference implementation + MSBooster RT+spectra) | 18,670 | 37,086 | 2.0 | ours unavailable |
Rescoring lifts both tools (+11% / +14%) and the ratio does not move (85.6 -> 82.8%). So the MSFragger deficit is
not a scoring-function artefact that a rescorer repairs: the 3,600 the reference implementation-only peptides have no competitive PSM
in our output at all. Together with H3 (identical scores on shared peptides) this points at CONTENT/coverage of the
faint tail -- emission competition or missing fragments -- not at search-side handling. Our 8.6 PSMs per peptide vs
the reference implementation's 2.0 is the over-emission signature: Percolator's peptide-level FDR gives redundancy nothing.
**Decision (plan step 8):** the next attribution arms are content-side: (c) the emission-controlled arm (pre-registered
step 02: quality-gate 924k -> ~700k + the reference implementation-precursor-matched sub-arm, both engines, full FDR curves) and (d) the
oracle-fragment arm. The learned charge/mono predictor stays conditional on (d). A MSBooster-tolerant rerun (var mods
capped identically for both tools) is owed for the +17% predictions bring the reference implementation.

### STEP 02, sub-arm 1 -- PRE-REGISTERED (2026-09-02, 22:00): precursor-quality gate by isotope-envelope depth
`SPEXTRACT_MIN_ISOTOPES=k` keeps precursors with >= k isotope peaks (k=2 is today's default via
require_isotope_support; k=3 and k=4 are the arms). Emission falls with k; fragment sharing is untouched.
Prediction if the MSFragger deficit is EMISSION COMPETITION: the spx/dt MSFragger ratio rises as emission falls
toward the reference implementation's 700k (peptides@1% hold or rise while spectra drop). Prediction if it is FAINT-TAIL CONTENT: the
ratio is flat or falls (the gate removes exactly the faint precursors). Gates: Sage set overlap vs the apex arm,
MSFragger raw walk vs msf_dt, emission count; both engines before any default change.

### STEP 02, sub-arm 2 -- PRE-REGISTERED (2026-09-02, 22:03): the reference implementation-precursor-MATCHED emission
`SPEXTRACT_PRECURSOR_LIST=<rt_sec mz z>` keeps only our precursors that match a listed the reference implementation dataset D precursor
(|dRT| <= 10 s, |dm/z| <= 10 ppm, same charge). The list is the reference implementation's own 700,434 pseudo-spectra (extracted from
its mzML; RT converted from minutes). This holds the PRECURSOR POPULATION fixed and varies only the spectra.
Prediction if the MSFragger deficit is EMISSION COMPETITION: with the same precursor set the spx/dt MSFragger ratio
moves toward parity (>= 95%). Prediction if it is PER-PRECURSOR CONTENT (fragment coverage of the faint tail): the
ratio stays near 85% even with matched precursors. Gates: emission count, Sage set vs the apex arm, MSFragger raw
walk vs msf_dt. Runs after sub-arm 1 (`bench_match.sh`, d15_matched).

### STEP 02 sub-arm 1 RESULT (2026-09-02, 22:34): EMISSION COMPETITION IS NOT THE MSFRAGGER MECHANISM
| dataset D arm | precursors kept | spectra | Sage @1% (vs apex 12,537) | MSFragger @1% (vs dt 13,932) | ratio | wall |
|---|---|---|---|---|---|---|
| apex (k=2, today's default) | 1,105,532 | 927,813 | 12,537 | 11,927 | 85.6% | 15:21 |
| **k=3** (>= 2 isotope partners) | 655,297 | **566,537 (−39%)** | 11,739 (−6.4%; only-apex 1,010 / only-k3 212) | **11,820 (−0.9%)** | 84.8% | 11:53 |
| k=4 | 399,837 | 353,689 (−62%) | 10,201 (−18.6%) | 10,538 (−11.6%) | 75.6% | 9:51 |
Pre-registered reading: cutting emission by 39% -- to BELOW the reference implementation's 700k -- moves the MSFragger ratio from 85.6%
to 84.8%. It does not rise. So the extra spectra are not what costs us on MSFragger (the emission-competition
hypothesis is refuted for the second time, now with a precursor-quality gate rather than merging), and the gate removes
real faint peptides on Sage (−6.4%). The deficit is per-precursor CONTENT of the faint tail, consistent with the
rescoring A/B (ratio unchanged after Percolator) and H3 (identical scores on shared peptides). Sub-arm 2 (matched
precursor population) is the direct test of that and is running. Side result: `SPEXTRACT_MIN_ISOTOPES=3` is a
legitimate speed/emission knob (−39% spectra, −23% wall, −0.9% MSFragger, −6.4% Sage) -- not a default.

### STEP 02 sub-arm 2 RESULT (2026-09-02, 22:49): SAME PRECURSORS, SAME 85% -- the deficit is per-precursor content
Reference: the reference implementation's 700,434 dataset D pseudo-spectra (RT/m/z/z). Of our 1,105,532 precursor hypotheses 445,008 match one
(10 s / 10 ppm / z); emitting only those gives 432,579 spectra.
| dataset D arm | spectra | Sage @1% | MSFragger @1% | MSFragger ratio vs dt |
|---|---|---|---|---|
| apex (all precursors) | 927,813 | 12,537 | 11,927 | 85.6% |
| **the reference implementation-precursor-matched** | 432,579 | 11,903 (103% of dt's 11,517) | **11,842** | **85.0%** |
| the reference implementation | 700,434 | 11,517 | 13,932 | — |
With the precursor POPULATION held to the reference implementation's, the MSFragger ratio does not move (85.0%). Together with sub-arm 1
(emission −39%: 84.8%), the rescoring A/B (Percolator: 82.8%) and H3 (identical scores on shared peptides), every
search-side and population-side explanation is now excluded. What is left is the CONTENT of the spectra we emit for the
precursors the reference implementation also emits: for ~3,600 peptides our spectrum of the same precursor does not reach MSFragger's
1% threshold while the reference implementation's does, and on Sage it does. Candidates, in test order: (1) the intensity reshaping
`corr_power=2` (a Sage win that MSFragger's rank-based hyperscore may dislike -- never measured on MSFragger; cp0 arm
launched); (2) fragment coverage of the faint tail (our ~55 fragments vs the reference implementation's ~500 per spectrum, top-150 cap in
MSFragger); (3) the oracle-fragment arm (true fragments at true charge/mono) to bound what content can give.
Note also: only 62% of the reference implementation's emissions have a matching precursor of ours -- the other 38% are either duplicates
on the reference implementation's side or precursors we never hypothesise (the MS1 funnel), a separate coverage question.

### CONTENT CANDIDATE 2 MEASURED (2026-09-03, 01:28): the faint tail is missing ~2-3 matched fragment ions per spectrum
MSFragger raw-walk TSVs (dataset D, best PSM per peptide+charge, target only; thresholds 15.2 ours / 14.0 dt):
| | apex arm (all precursors) | the reference implementation-precursor-matched arm |
|---|---|---|
| the reference implementation-identified / shared / **dt-only** | 14,974 / 11,037 / **3,937** | 14,974 / 11,367 / **3,607** |
| dt-only for which we HAVE a PSM (searched, below threshold) | **2,501 (64%)** | 1,567 (43%) |
| dt-only with NO PSM at all (never emitted or < min_frags) | 1,436 (36%) | 2,040 (57%) |
| dt-only, ours vs dt: hyperscore median | 12.5 vs 16.2 | 13.1 vs 16.2 |
| dt-only, ours vs dt: matched ions median (of 26 theoretical) | **7 vs 9** (deficit +2; p25 +1, p75 +4; same at every charge) | 6 vs 9 (deficit +3) |
| shared, ours vs dt: hyperscore / matched ions | 23.6 vs 22.9 / 12 vs 12 | 24.4 vs 22.4 / 12 vs 12 |
Reading: on the peptides we both find, our spectra are as good as the reference implementation's. On the ~3,900 we miss, two thirds
are searched and fall short by ~2 of ~9 fragment ions -- the faint tail is under-covered: the weaker fragments of
faint precursors are gated out (`gate:min_correlation` 0.3, `min_correlation_points`) or never traced
(`trace:ms2_noise_threshold_int` 10, `ms2_chrom_peak_snr` 1.0). The remaining third we never emit at all (the MS1
funnel: no hypothesis or < min_frags). This is the first mechanistic, per-peptide account of the MSFragger deficit.
**Pre-registered arms (running next, both engines):** (a) `gate:min_correlation 0.2`, (b) `trace:ms2_noise_threshold_int 5`,
(c) both. Prediction: the dt-only matched-ion deficit shrinks toward 0 and the MSFragger ratio rises above 90%;
falsifier: the deficit stays at +2 (then the missing ions are not in our traces at all -> MS2 aggregation / detector).

### CONTENT CANDIDATE 1 FALSIFIED (2026-09-03, 01:52): corr_power=2 helps MSFragger too
`-assembly:corr_power 0` on the apex arm: Sage 11,502 (−8.3% vs 12,537; only-cp2 1,640 / only-cp0 605) and
**MSFragger 11,429 (−4.2% vs 11,927)**. The correlation-power intensity reshaping is a win on both engines (first
MSFragger measurement of it; the July +8-10% was Sage-only). Intensity structure is not what MSFragger dislikes in our
faint-tail spectra. Remaining content candidate: fragment coverage (the +2 matched-ion deficit) -- arms running.

### COVERAGE ARMS FALSIFIED (2026-09-03, 03:10): the missing fragment ions are not in our traces at all
| dataset D arm | spectra | Sage vs apex 12,537 | MSFragger vs apex 11,927 | dt-only: ours vs dt matched ions |
|---|---|---|---|---|
| `gate:min_correlation 0.2` | 928,572 | 12,509 (145/117) | 11,935 | 7 vs 9, deficit +2 |
| `trace:ms2_noise_threshold_int 5` | 927,813 | **12,537, symdiff 0 (byte-identical)** | 11,927 | 7 vs 9, deficit +2 |
| both | 928,572 | 12,509 | 11,935 | 7 vs 9, deficit +2 |
Loosening the correlation gate adds nothing to the faint tail, and the intensity threshold 5 vs 10 changes NOTHING
(the SNR / min-length gates already decide). Per the pre-registered falsifier, the ~2 missing ions per faint spectrum are
not in our fragment TRACES: they are peaks that never become a trace (single-frame or sub-SNR) -- exactly the population
the reference implementation's ~500-peak spectra carry and our ~55-trace spectra do not. Next arm (content candidate 3): backfill each
pseudo-spectrum with the raw picked MS2 peaks of the precursor's apex frame inside its IM band (top-N not already
present), the DIA-Umpire/the reference implementation-style "peak-level" spectrum, both engines; falsifier: dt-only deficit stays +2.

### CONTENT CANDIDATE 3 -- PRE-REGISTERED (2026-09-03, 03:35): raw apex-frame peak backfill
`SPEXTRACT_BACKFILL_RAW=N` adds to each pseudo-spectrum the N most intense raw picked MS2 peaks of the frame nearest
the precursor's RT, inside its 1/K0 band (delta_im), not already present within 10 ppm -- untraced, single-frame
content (the population the coverage arms showed we never trace). Arms: N=50, N=150. Predictions: if the faint tail's
missing ions are these peaks, the dt-only matched-ion deficit (+2) shrinks and the MSFragger ratio rises (>= 90%);
Sage may fall if the added peaks are noise (the 150-peak search cap protects MSFragger more than Sage). Falsifier:
deficit stays +2 with N=150 -> the missing ions are not in the picked frame either (pick-level / MS2 aggregation).

### Backfill arm, first run INVALID (2026-09-03, 04:42) -- and the never-searched third attributed
N=50 and N=150 came back byte-identical to apex (Sage symdiff 0, MSFragger 11,927): the backfill never fired because
`detectTraces_` partitions the window's PeakMap into bands and clears it, so the map is empty by assembly time. Fixed
(a retained copy when the arm is on); rerun in flight. Not a result.
**MS1 funnel (funnel_dtonly.py on the apex arm's hypothesis dump, 1,105,532 hypotheses):** of the 1,436 the reference implementation-only
peptides for which we have NO PSM, 654 (46%) have a precursor hypothesis of ours within 10 s / 10 ppm / z (918 = 64%
within 30 s) -> hypothesised but no spectrum reached the search (assembly: < min_frags or no correlated fragments);
~500 (35%) have no hypothesis at all (MS1 trace / envelope). So the 3,937 dt-only peptides split ~64% searched-but-
under-covered / ~17% hypothesised-not-emitted / ~19% never hypothesised.

### ASSEMBLY-LOSS ARM -- PRE-REGISTERED (2026-09-03, 04:45): `assembly:min_fragments` 2
17% of the the reference implementation-only peptides have a precursor hypothesis of ours but no spectrum reached the search. The gate that
drops a hypothesised precursor at assembly is `assembly:min_fragments` (a spectrum with fewer correlated fragments is
not emitted). Arm: min_fragments 2 (the default is already 3), both engines. Prediction: the dt-only "we have SOME PSM"
fraction rises from 64% and MSFragger gains if those spectra are identifiable; falsifier: the added spectra are
unidentifiable noise (Sage/MSFragger flat or down, emission up). Runs after the backfill rerun (`bench_minfrag.sh`).

### CONTENT CANDIDATE 3 CONFIRMED (2026-09-03, 05:39): raw apex-frame peak backfill lifts MSFragger to 90% -- at a Sage cost
| dataset D arm | Sage @1% (vs apex 12,537) | MSFragger @1% (vs apex 11,927) | ratio vs dt | shared with dt / dt-only | dt-only matched-ion deficit |
|---|---|---|---|---|---|
| apex | 12,537 | 11,927 | 85.6% | 11,037 / 3,937 | +2 |
| **backfill N=50** | 12,208 (−2.6%; 1,215 / 886) | **12,539 (+5.1%)** | **90.0%** | **12,188 / 2,786** | +3 (on the smaller remaining set) |
| backfill N=150 | 11,190 (−10.7%) | 12,042 (+1.0%) | 86.4% | 11,437 / 3,537 | +2 |
The untraced peaks of the precursor's apex frame ARE the missing content: adding the 50 most intense of them inside the
IM band converts 1,151 the reference implementation-only peptides into shared ones and moves the MSFragger ratio from 85.6% to 90.0% --
the first lever that moves it since the calibration fix. It costs Sage 2.6% (the added peaks are noise for the
spectra that were already rich: Sage's scoring is hurt by them, MSFragger's top-150 cap is not); N=150 is too much
for both. Not a default as measured -- the engines disagree. **Pre-registered next arm:** backfill only the faint
tail: `SPEXTRACT_BACKFILL_RAW=50` + `SPEXTRACT_BACKFILL_MAXFRAGS=K` (apply only to spectra with <= K assembled
fragments; K = 50 and 150 -- the first launch of this arm collided with a second chain on the same output paths at
06:18 and was stopped and relaunched cleanly at 06:19). Prediction: MSFragger keeps most of +5% while Sage returns to ~12,500; falsifier: Sage still
drops -> the noise cost is in the faint spectra themselves, and the selection must be per-peak (IM/RT co-elution of the
raw peak, i.e. MS2 aggregation across neighbouring frames, C6).

### ASSEMBLY-LOSS ARM FALSIFIED (2026-09-03, 06:16): `min_fragments 2` changes nothing
928,378 spectra (+565), Sage and MSFragger identical to apex, dt-only PSM coverage unchanged (64%). The hypothesised-
but-not-emitted precursors (17% of the reference implementation-only) have NO correlated fragment trace at all, not too few -- the same
untraced-content problem as the searched tail. **Pre-registered follow-up:** `SPEXTRACT_BACKFILL_EMPTY=1` with
`SPEXTRACT_BACKFILL_RAW=50`: a hypothesised precursor with no correlated trace still gets a spectrum made of its apex
frame's IM-band peaks (min_fragments honoured). Prediction: the dt-only "no PSM at all" count (1,436) falls and
MSFragger rises further; falsifier: those spectra are unidentifiable (counts flat, emission up). Runs after the
faint-tail-only backfill arm.

### FAINT-TAIL-ONLY BACKFILL (2026-09-03, 07:11): removes the Sage cost, but also most of the MSFragger gain
`SPEXTRACT_BACKFILL_RAW=50` applied only to spectra with <= K assembled fragments:
| dataset D arm | Sage @1% (vs apex 12,537) | MSFragger @1% (vs 11,927) | dt-only deficit |
|---|---|---|---|
| unconditional N=50 (05:39) | 12,208 (−2.6%) | **12,539 (+5.1%)** | +3 on the remaining set |
| **K=150** | **12,724 (+1.5%; only-apex 11 / only-K150 198)** | 12,000 (+0.6%) | +2 |
| K=50 | 12,591 (+0.4%; 6 / 60) | 11,961 (+0.3%) | +2 |
Reading: the MSFragger gain of the unconditional arm comes from backfilling the RICH spectra (73% of spectra sit at
the 500-cap; share-all gives even faint precursors hundreds of fragments), which is exactly where Sage pays. Gating by
fragment count therefore does not reach the the reference implementation-only spectra (deficit still +2). K=150 is a clean small win on
both engines (+1.5% / +0.6%, superset-like) -- a candidate default pending dataset A/dataset B + entrapment, not a mechanism.
**Pre-registered next arm:** gate the backfill by PRECURSOR intensity instead (`SPEXTRACT_BACKFILL_MAXQ=q`: only
precursors below the window's q-quantile of MS1 intensity), N=50, q = 0.5 and 0.25. Prediction: MSFragger keeps most
of +5% (the faint precursors are the dt-only ones) while Sage stays >= 12,500; falsifier: the gain tracks the rich
spectra regardless of precursor intensity (then MSFragger simply likes more peaks, and the two engines want
different spectra -- a per-engine output option, not a default).

### BACKFILL-EMPTY FALSIFIED (2026-09-03, 08:35): 53,410 extra spectra buy nothing on either engine
`SPEXTRACT_BACKFILL_EMPTY=1` (a hypothesised precursor with no correlated fragment trace still gets a spectrum built
from its apex frame's IM-band peaks) + `SPEXTRACT_BACKFILL_RAW=50`. Emission rises 927,813 -> **981,223 (+5.8%)**, so
the feature does what it says; the identifications do not follow.
| dataset D arm | spectra | Sage @1% (apex 12,537) | MSFragger @1% (apex 11,927; dt 13,932) | dt-only |
|---|---|---|---|---|
| empty + K=0 (backfill all) | 981,223 | 12,218 (−2.5%) | 12,600 (+5.6%, 90.4% of dt) | 3,309; some PSM 58% |
| empty + K=60 (faint tail) | 981,223 | 12,593 (+0.4%) | 11,982 (+0.5%) | 3,892; some PSM 64% |
| no-empty N=50 (reference) | 927,813 | 12,208 (−2.6%) | 12,539 (+5.1%) | some PSM 64% |
| no-empty K=150 (reference) | 927,813 | **12,724 (+1.5%)** | 12,000 (+0.6%) | +2 |
Reading: at equal backfill setting the empty arm reproduces the no-empty arm to within replicate noise on both engines
(12,600 vs 12,539 MSFragger; −2.5% vs −2.6% Sage), and the K-gated empty arm is *worse* than K=150 alone (12,593 vs
12,724 Sage). The 53k trace-free precursors are not identifiable: the pre-registered falsifier ("counts flat, emission
up") is met exactly. **Verdict: `SPEXTRACT_BACKFILL_EMPTY` stays off and is not a default candidate.** It also
confirms the assembly-loss finding from the other side -- the the reference implementation-only precursors we never emit are not lost to
a fragment-count threshold, they have no correlated MS2 signal for us to find at all.

### PRECURSOR-INTENSITY-GATED BACKFILL (2026-09-03, 09:05): prediction held; q=0.5 is the best point on the frontier
`SPEXTRACT_BACKFILL_RAW=50` applied only to precursors below the window's q-quantile of MS1 intensity. Emission is
unchanged (927,813 spectra in every arm), so this is purely a content lever.
| dataset D arm | Sage @1% (apex 12,537) | MSFragger @1% (apex 11,927; dt 13,932) | dt-only | shared w/ dt |
|---|---|---|---|---|
| **q=0.5** | 12,442 (−0.8%; only-apex 460 / only-q 365) | **12,492 (+4.7%, 89.7% of dt)** | 2,940 | **12,034** |
| q=0.25 | 12,504 (−0.3%; 141 / 108) | 12,008 (+0.7%) | 3,913 | 11,061 |
| N=50 unconditional | 12,208 (−2.6%) | 12,539 (+5.1%) | — | — |
| K=150 fragment-gated | **12,724 (+1.5%)** | 12,000 (+0.6%) | — | — |
Reading: the pre-registered falsifier is NOT met. Backfilling only the faint HALF of precursors retains 92% of the
unconditional MSFragger gain (+4.7 of +5.1 points) while cutting the Sage cost by two thirds (−0.8% vs −2.6%), so the
gain does track precursor intensity rather than peak count alone. It is concentrated in the second quartile band:
q=0.25 (the faintest quarter only) collapses to +0.7%, i.e. the very faintest precursors' apex frames carry nothing
identifiable -- the same wall the backfill-empty arm hit. Sage degrades monotonically with backfill volume in all
arms. q=0.5 also gives the highest overlap with the reference implementation measured to date (12,034 shared, dt-only down to 2,940),
but the "no PSM at all" residue is 1,406, unmoved from baseline: this is rich-spectrum content, not the deficit
mechanism.
**Two candidate defaults now stand, and they are not the same trade:** K=150 is strictly better than baseline on both
engines but small (+1.5% / +0.6%); q=0.5 is much larger on MSFragger (+4.7%) at a −0.8% Sage cost. Per the both-engines
rule neither is adopted on one file. **Pre-registered confirmation arm:** both settings on dataset A and dataset B, both engines,
plus entrapment FDR. Adopt as default only if the dataset D sign holds on both files AND entrapment stays <= 1.4%
(shipping arm: dataset D 1.37%, dataset A 0.95%, dataset B 1.32%); a Sage regression on dataset A/dataset B larger than dataset D's −0.8% kills q=0.5 and
leaves K=150 as the only candidate.

### ISOTOPE DUPLICATION: HYPOTHESIS (a) CONFIRMED, (b) REFUTED (2026-09-03, 09:35)
Full analysis in [ISOTOPE-DUPLICATION-2026-09-03.md](ISOTOPE-DUPLICATION-2026-09-03.md). 31.4% of
emitted spectra have a co-eluting same-charge partner 1-3 isotope steps below (decoy offset 3.8%,
excess 27.7%); removing them leaves 671,270 spectra against the reference implementation's 700,434, i.e. isotope-offset
duplication accounts for the entire emission excess. But they are NOT mergeable: content cosine
0.505 vs 0.446 for an arbitrary co-eluting neighbour (only 6.5% above 0.9), and every id-agnostic
collapse rule loses 5,011-5,173 of 14,944 peptides (oracle still loses 1,387). "Keep the lightest" is
backwards -- where only the heavy member is identified, Sage says the mass we reported was already
correct 70% of the time. Mechanism: `findPartner` skips `used[]` peaks, so a leftover one step above a
consumed run becomes its own monoisotope (heavy members carry fewer isotope partners, 3.53 vs 4.59).
Upstream ownership gate pre-registered; downstream merging stays falsified.

### BACKFILL CONFIRMATION ARM: BOTH CANDIDATES FALSIFIED (2026-09-03, 11:47) -- backfill stays OFF
Pre-registered arms `q=0.5` and `K=150` re-run on dataset A and dataset B, both engines, plus entrapment.
| file | arm | Sage @1% | vs apex | MSFragger @1% | vs apex |
|---|---|---|---|---|---|
| dataset D | apex | 12,537 | -- | 11,927 | -- |
| dataset D | q=0.5 | 12,442 | −0.8% | 12,492 | +4.7% |
| dataset D | K=150 | 12,724 | **+1.5%** | 12,000 | +0.6% |
| dataset A | apex | 10,789 | -- | 11,444 | -- |
| dataset A | q=0.5 | 10,516 | **−2.5%** | 11,620 | +1.5% |
| dataset A | K=150 | 10,618 | **−1.6%** | 11,483 | +0.3% |
| dataset B | apex | 10,156 | -- | 9,738 | -- |
| dataset B | q=0.5 | 10,134 | −0.2% | 10,178 | +4.5% |
| dataset B | K=150 | 10,166 | +0.1% | 9,751 | +0.1% |
Entrapment is flat everywhere (dataset A apex 0.91% / q05 0.91% / K150 0.97%; dataset B 1.13% / 1.13% / 1.14%),
so neither arm is bought with false positives -- but neither survives its own falsifier:
* **K=150 dies on the sign rule.** dataset D said +1.5% on Sage; dataset A says −1.6%. The "strictly better on both
  engines" claim was a one-file artefact.
* **q=0.5 dies on the magnitude rule.** The pre-registration killed it if any dataset A/dataset B Sage regression
  exceeded dataset D's −0.8%; dataset A is −2.5%. The MSFragger gain does replicate (+4.7 / +1.5 / +4.5), so the
  two engines genuinely want different spectra -- but that makes it a per-engine option at best, not a
  default, and nothing in the evidence says which engine to optimise.
**Verdict: raw apex-frame backfill stays OFF by default in every form tested** (unconditional,
fragment-gated, precursor-intensity-gated, empty-precursor). The MSFragger deficit is not closed by
adding peaks to spectra we already emit.

### NEW DEFAULT: `charge:min_charge` = 2 (2026-09-03, 15:05) -- confirmed on all three files, both engines
Singly-charged precursor hypotheses are ~30% of emission and ~1.7% of peptides. A z=1 pseudo-spectrum
is identified on 0.42% of its own spectra against 15.3% for z=2; tryptic peptides are essentially
never 1+ in ESI and `charge:scoring count` breaks ties toward the LOW charge, so the mis-assignments
land there. Dropping them removes their share of the multiple-testing burden too, so peptides go UP.
| file | spectra | Sage @1% | MSFragger @1% | entrapment |
|---|---|---|---|---|
| dataset D | 927,813 -> **656,254** (−29.3%) | 12,537 -> **12,642** (+0.8%) | 11,927 -> **12,073** (+1.2%) | 1.28% -> 1.38% |
| dataset A | 1,228,875 -> **863,319** (−29.7%) | 10,789 -> **11,084** (+2.7%) | 11,444 -> **12,061** (+5.4%) | 0.91% -> 0.99% |
| dataset B | 844,755 -> **586,069** (−30.6%) | 10,156 -> **10,294** (+1.4%) | 9,738 -> **9,822** (+0.9%) | 1.13% -> 1.28% |
dataset D wall time 15:21 -> 12:48 (−17%). Every pre-registered adoption condition is met: peptides rise on
BOTH engines on ALL THREE files, every entrapment estimate stays inside the corresponding apex 95% CI
and below the 1.4% bound. Against the reference implementation, MSFragger goes 85.6 -> 86.7% (dataset D), 83.7 -> 88.3% (dataset A),
91.8 -> 92.6% (dataset B); Sage on dataset D is 109.8% of the reference implementation. This is the first change in the emission line
that improves both engines while cutting emission AND runtime.
**Adopted as the default.** `charge:min_charge=1` restores the old behaviour; 3+ is catastrophic
(charge 2 carries 56.9% of peptides: dataset D Sage 5,088). The isotope collapse is NOT adopted in any form
-- on top of this gate it turns +0.8% into −3.8% (Sage) and +1.2% into −2.5% (MSFragger).


### WINDOW-LOOP PERFORMANCE LINE (2026-09-03, evening) -- one adopted, one held
| dataset D, 100 threads | wall | peak RSS | window-loop occupancy | Sage @1% |
|---|---|---|---|---|
| morning default (charge gate) | 12:53 | 103 GB | 64.8x | 12,642 |
| task pool alone | 15:50 | 152 GB | 86.2x | 12,642 (set identical) |
| **+ scoring gate reads a parallel array** | **7:17** | 164 GB | 80.7x | **12,642 (set identical)** |
| + integer detector (`trace:detector=integer`) | **6:53** | **105 GB** | 67.6x | 12,466 (−1.4%) |
**ADOPTED: the parallel-array gate and the task pool.** The gate rejected 99.4% of the fragments it
visited (867 billion visits for 5.5 billion survivors, measured) and read the field it needed out of
a 96 B record: ~83 TB of memory traffic per run. Set-identical output, replicated three times. The
task pool alone was a regression only because it amplified that traffic; with the gate fixed it is
neutral-to-positive, and the in-flight cap I briefly added on the wrong diagnosis is removed again.
**NOT ADOPTED (yet): the integer detector.** Faster and 36% lighter, but −1.4% Sage peptides against
the reference, with a symmetric ~8% set churn that is the residual of a faithful-not-identical
reimplementation. The only semantic difference between its 12,650 and 12,466 versions is whether a
frame with only sub-noise peaks counts as empty (the OpenMS rule) or as a miss; that A/B, MSFragger,
and entrapment are running. Adoption needs: peptides within the replicate spread on BOTH engines,
entrapment inside the charge-gate arm's interval, and dataset A replication of the runtime.


### WINDOW-LOOP PERFORMANCE, ROUND 2 (2026-09-04) -- three output-neutral changes, THREE files

Measured old-vs-new with both arms on the SAME host, three files concurrently on three machines,
all running one binary from the shared ceph install (`/path/to/shared/spextract`) so the
arms cannot differ by build. **Every arm is digest-identical to its base**, which is the gate:
these changes are output-neutral by construction, so peptides are identical by definition and no
search was run.

| | dataset D (node-1, 128c) | dataset A (node-2, 224c) | dataset B (data, 224c) |
|---|---|---|---|
| wall base -> perf (SINGLE PAIR -- see the retraction below) | 6:31 -> 6:07 | 8:36 -> 8:09 | 5:40 -> 5:22 |
| process peak RSS | 85.7 -> 83.3 GB (−2.8%) | 119.6 -> 114.2 GB (−4.5%) | 91.3 -> 88.6 GB (−2.9%) |
| **RSS at end of window loop** | 65.0 -> **48.3 GB** (−25.8%) | 85.6 -> **55.1 GB** (−35.6%) | 70.0 -> **49.8 GB** (−28.9%) |
| window loop wall | 209.5 -> 192.7 s (−8.0%) | 317.6 -> 303.5 s (−4.4%) | 175.6 -> 162.7 s (−7.3%) |
| window-loop occupancy | 67.6 -> 69.9x | 56.2 -> 61.3x | 56.9 -> 63.3x |
| system time | 7:12 -> 4:55 (−32%) | 11:34 -> 8:08 (−30%) | 8:17 -> 6:05 (−27%) |
| EPD MassTrace (ledger) | 7,502 -> 0.08 MB | 10,578 -> 0.08 MB | 9,192 -> 0.08 MB |
| sum of per-structure peaks | 36.2 -> 28.3 GB | 58.5 -> 46.9 GB | 41.3 -> 31.2 GB |
| digest vs base | **IDENTICAL** | **IDENTICAL** | **IDENTICAL** |

### RETRACTION (same day): the wall-clock claim does NOT survive replication

The single-pair table above was replicated as an interleaved **base/clean/base/clean** on dataset B,
`data`, with the load recorded before each run. Four runs, all four digest-identical:

| run | load before | wall | window loop | window-loop CPU | RSS at end of loop |
|---|---|---|---|---|---|
| base rep1 | 3.4 | 5:29 | 174.3 s | 9,850.8 | 69,999 MB |
| clean rep1 | 29.5 | 5:19 | 159.8 s | 10,269.2 | 48,215 MB |
| base rep2 | 41.8 | 5:25 | 168.3 s | 9,988.9 | 69,536 MB |
| clean rep2 | 44.7 | 5:36 | 173.8 s | 10,400.5 | 50,135 MB |

**Wall clock: base mean 5:27, clean mean 5:27.5. There is no wall-clock difference.** The within-arm
spread (base 5:25-5:29, clean 5:19-5:36) is larger than the −5% the single pairs appeared to show.
The window loop is 171.3 s vs 166.8 s (−2.6%) with overlapping ranges. **Window-loop CPU-seconds are
consistently ~4% HIGHER for the new code** (9,850/9,988 vs 10,269/10,400), which is the one timing
signal that reproduces, and it points the wrong way: per-trace `detectPeaks` calls and the gate's
extra build pass cost CPU that batching did not.

**What DOES replicate, tightly, is the memory.** RSS at the end of the window loop: 69,999 / 69,536
(base) vs 48,215 / 50,135 (clean) -- **−20.6 GB, −29.5%**, with the two reps of each arm within
2 GB of each other. That matches the −25.8% / −35.6% / −28.9% seen on the first pass and is
corroborated structurally by the EPD ledger lines falling 7.5-10.6 GB -> 0.08 MB, which is a change
in what exists, not a measurement.

**So the defensible claim is: identical output, ~30% less memory in the window loop, no measured
wall-clock change, and possibly ~4% more CPU.** The earlier "−5 to −6% wall on three files" was
three single pairs taken while another user's load drifted underneath them, and it is withdrawn.

#### THREE-PAIR REPLICATION, dataset D on node-2 (the settled numbers)

Six runs, interleaved base/clean x3, load recorded before each, **all five digest comparisons
identical**. Every metric below has NON-OVERLAPPING ranges between the arms, which is the bar the
single pairs failed:

| metric | base mean [range] | clean mean [range] | delta |
|---|---|---|---|
| wall | 321.3 s [320-323] | 312.7 s [310-317] | **−2.7%** |
| window-loop CPU-s | 9,228.9 [9,182-9,309] | 9,689.7 [9,624-9,802] | **+5.0%** |
| RSS at end of loop | 63,880 MB [63,715-64,137] | 46,922 MB [46,459-47,313] | **−26.5%** |
| process peak RSS | 81.3 GB [80.9-81.9] | 77.9 GB [76.4-78.9] | **−4.1%** |

**Final reading, both files together.** The memory result is large, tight and certain: ~27-30% off
the window loop on both files, ~4% off process peak, arms never overlapping. The wall-clock result
is a REAL BUT SMALL win on dataset D (−2.7%, half what the single pairs claimed) and ZERO on dataset B --
so "up to ~3% on one file, none on another", not "−5 to −6% on three files". The CPU cost is real
and reproduces on both files (+5.0% dataset D, +4.2% dataset B): per-trace `detectPeaks` calls and the gate's
extra build pass cost CPU that batching did not. **The trade is ~27% of the window loop's memory
for ~5% more CPU at roughly unchanged wall time.** For a tool whose concurrency is bounded by the
free-RAM admission gate, that is worth taking -- but it is a memory change and must be cited as one.

**THE MACHINES ARE SHARED, AND THIS WAS NEARLY MISREPORTED.** The peer nodes were called "idle"
on the strength of `nproc` and `free` alone -- load average was never checked. They were not idle:
node-1 was carrying another user's three python jobs at load 42, node-2 three `jackhmmer` at
load 21. A later unpaired run of the shipping binary came back at dataset D 6:28 / dataset A 8:40, i.e. back at
baseline, purely because it ran in a more contended window than the arms it was being compared to.
**Rule from here: quote a wall-clock delta only from arms run BACK TO BACK on one node, check
`uptime` before and after, and prefer the interleaved A/B/A/B below to a single pair.** The table
above satisfies the back-to-back condition; the unpaired 11:02 runs did not and are not quoted.

**Read the memory row that matters.** Process peak RSS barely moves (2-5 GB) because the peak is
set outside the window loop; the window loop itself now ENDS 26-36% lighter. That is the number
that governs how many windows can be in flight, not the process peak.

**CPU-seconds RISE slightly on two of three files** (dataset A 5:40:40 -> 5:51:55, dataset B 3:21:01 ->
3:23:16) while wall falls everywhere. That is not a contradiction and it is not free: occupancy
rose 2-6x, so the tool is using more of the machine for less elapsed time. Do not quote the CPU
column as an improvement.

The three changes:
1. **Arena reserves.** `TraceStore::absorb` appended without reserving into a destination that had
   just been swapped empty, so 12 band merges and 48 chunk merges regrew geometrically.
2. **Streamed valley splitting.** Was: build every trace of a chunk as an OpenMS `MassTrace`, then
   split the whole chunk, so the input payload and the entire split output were alive together.
   Now one trace is in flight at a time. `ElutionPeakDetection` is hoisted to chunk scope -- doing
   its four `setValue` calls per trace would have been ~4e4 `Param` round-trips per chunk and would
   have turned this into a slowdown.
3. **Quantised fragment RT gate.** The gate rejects 99.4% of what it visits and was reading an 8 B
   double to do it: 867e9 visits = 6.9 TB per dataset D run. A `uint16` bucket array makes the reject
   2 B; the exact test is unchanged and only runs on survivors. The score-gate counters are
   byte-identical between arms (862,781,140,969 visits, 5,467,207,106 survivors), which is the
   gate-level proof the candidate set did not move.

**Adversarial review (codex/kimi/vibe) — all three verdicts OUTPUT-NEUTRAL, with two real defects
found in change 3 that the safety proof had missed.** vibe's proof assumed finite intermediates;
codex and kimi independently found that a denormal `delta_rt` gives `1/delta_rt = inf`, fragment
RTs spanning +/-DBL_MAX give an infinite span, and a NaN `pc.rt` reaches `(int)NaN` -- undefined
behaviour that on x86 rejects every fragment where the old test passed all of them. Fixed by making
`bucketOf` NaN-safe and degrading the whole field to disabled (`inv = 0` => every bucket 0 => every
`dq` 0 => the exact test alone decides), which costs no branch in the scan. Also fixed: the
`MEM_EPD_MT` guard was constructed after the loop and measured an already-destroyed payload, and
`absorb()` left `bins` short if a binless source followed a binful one.

**REJECTED on evidence, not opinion:**
* **Slab compaction to above-noise peaks** (was the largest proposed memory item). The existing
  `[mem] window ... seeds N of P peaks` line reads **100.0%** -- every peak is above noise, so
  there is nothing to drop and the slab's ~12 GB is irreducible. Killed for free, from a number the
  binary already printed.
* **Moving the seed sort into the band tasks.** Priced by its own supporting evidence at <1% of CPU,
  and it re-partitions seeds across bands -- the exact failure mode that cost 3.3% of peptides once.


## v0.3.0 FINAL BENCHMARK — six datasets, shipped defaults (2026-09-04)

Every dataset in the cohort, run with `spextract -in <file>.d -out pseudo.mzML -threads 100` and
**nothing else**. If a figure below needed a flag, the defaults would be wrong; that is now a test.
Three nodes in parallel, both engines, `spx:detector=integer` and
`spx:require_isotope_support=1` recorded in every output.

| file | spectra | wall | peak RSS | window-loop occupancy | Sage @1% | MSFragger @1% |
|---|---|---|---|---|---|---|
| dataset A | 862,716 | 8:26 | 108.7 GB | 62.8x | 10,909 | 11,935 |
| dataset B | 585,503 | 5:37 | 84.6 GB | 61.3x | 10,272 | 9,691 |
| dataset C | 723,314 | 6:27 | 90.3 GB | 65.5x | 12,149 | 12,516 |
| dataset D | 655,776 | 5:21 | 79.0 GB | 60.7x | 12,482 | 12,337 |
| dataset E | 542,533 | 5:38 | 67.9 GB | 61.4x | 11,217 | 10,585 |
| dataset F | 597,267 | 5:59 | 73.8 GB | 66.7x | 11,362 | 11,049 |

**dataset C and dataset F had never been benchmarked before**; the cohort had only ever been exercised on
dataset A/dataset B/dataset D. Both behave like the rest, which is the first evidence that the defaults generalise
beyond the three files every decision in this file was made on.

Ranges across the cohort: wall 5:21-8:26, peak RSS 67.9-108.7 GB, occupancy 60.7-66.7x, emission
0.54-0.86 M spectra. Memory tracks acquisition size, and 108.7 GB on dataset A is what sets the
"80-125 GB" requirement now stated in the README.

**Read the two engines separately, as always.** Sage and MSFragger disagree on which files are
easy: dataset C is MSFragger's best (12,516) and Sage's second (12,149), while dataset B is the weakest on
both. Do not average them, and do not quote one as "the" peptide count.

Timings are NOT comparable between rows: the three nodes differ (128 vs 224 cores) and are shared
with other users, whose load is recorded in each run's own output file. The numbers that ARE
comparable across rows are spectra, peptides and peak RSS.
