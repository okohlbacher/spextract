# The recall gap is low-abundance peptides we EMIT but fail to IDENTIFY (dataset D, matched Sage)

Direct matched comparison (same dataset D .d, both searched with sage_deiso.json / human.fasta / ±25 ppm):

- the reference implementation **11,517** peptides, ours **10,036**. **SHARED 8,485, MISS 3,032, ours-only 1,551.**
- Not a subset: we find 1,551 peptides the reference implementation misses (our cleaner spectra win the high-abundance
  end); we miss 3,032 it finds. Net the reference implementation +1,481 (~87% recall vs the reference implementation on dataset D).

## What distinguishes the peptides we MISS — abundance, and nothing else

`bench/coverage.py`, MISS vs SHARED medians (ratio):

| property | ratio | | property | ratio |
|---|---|---|---|---|
| **precursor intensity** | **0.31** | | peptide length | 0.93 |
| ms2_intensity | 0.35 | | precursor m/z | 0.95 |
| matched_intensity_pct | 0.49 | | charge | 1.00 |
| matched_peaks | 0.75 | | ion mobility | 0.98 |
| hyperscore | 0.77 | | RT | 0.93 |

**The peptides we miss are ~3× lower abundance. Length, charge, m/z, mobility, RT do NOT
discriminate at all** (ratios 0.93–1.00). Missed peptides concentrate at hyperscore 15–25 (6.1×
enriched over SHARED) — borderline IDs. It is purely a low-abundance / borderline-ID gap.

## A/B split (chance-floor controlled): it's ASSEMBLY/ID, not extraction

For each missed peptide, did we emit a spectrum at its (m/z, RT, IM) coordinate? Decoy floor =
same RT/IM, physically impossible m/z offset (+7.3), to subtract coincidental matches (the control
the retracted MS1 funnel lacked).

| RT window | A not emit | B emitted | A% | decoy-floor B% |
|---|---|---|---|---|
| **0.2 min** | 318 | 2,714 | **10.5%** | **16.0%** |
| 0.5 min | 152 | 2,880 | 5.0% | 27.4% |
| 1.0 min | 75 | 2,957 | 2.5% | 37.3% |
| positive control (SHARED) @0.5 | 78 | 8,407 | 0.9% | — |

At the tightest, most-reliable window (chance floor 16%): B_observed 89.5%, chance-corrected
**~87% of missed peptides are GENUINELY EMITTED but not identified**; only ~13% are truly not
emitted. The SHARED positive control (0.9% A) validates the method.

## Verdict + direction

**The recall gap is an IDENTIFICATION-quality problem on low-abundance precursors we already emit,
NOT an extraction/sensitivity problem.** We emit a spectrum for ~87% of the low-abundance peptides
the reference implementation identifies; ours just don't cross the 1% FDR line. Since abundance is the only
discriminator, the fix must make our emitted spectra for FAINT precursors more identifiable —
candidates to test:
1. **Precursor annotation on faint precursors** — wrong charge/monoisotope → Sage searches the wrong
   mass → no ID even though a spectrum exists there (charge is 74.6% vs the reference implementation 82.6%). The A/B
   test matched on m/z only, so a right-m/z / wrong-mass spectrum still counts as "emitted".
2. **Fragment richness/cleanliness at low abundance** — do our faint-precursor spectra carry fewer
   of the true fragments (MS2 trace sensitivity) or more co-isolation dilution than the reference implementation's?
3. **Per-spectrum FDR burden** — we emit 1.53M spectra (924k gated); more spectra = heavier
   multiple-testing at a per-spectrum threshold, which pushes borderline low-abundance IDs under q.
   The gate (924k) should already help this; worth re-measuring recall on the gated arm.

Next: take OUR emitted spectra at the ~2,714 emitted-but-unidentified coordinates and check whether
they carry the right precursor mass and enough correct fragments — separating cause 1 from cause 2.

---

## CAUSE SEPARATION (decode our vs the reference implementation spectra at the missed coordinates)

`recall_cause2.py` — for each missed peptide, locate OUR emitted spectrum (right precursor mass) and
the reference implementation's, compute b/y ladder coverage the SAME way, with a SHARED positive control and a
decoy-ladder (+13 Da) chance floor. the reference implementation's patched.mzML has stale byte offsets (patched after
indexing) so its spectra are decoded by a sequential pass; our mzML index is valid.

**Positive control (SHARED, both identify):** our ladder cov 0.401, the reference implementation 0.462 — both high,
method validated.

**Missed peptides (n=2,837 with our right-mass spectrum):**

| signal | ours | the reference implementation | chance | verdict |
|---|---|---|---|---|
| precursor mass+charge right | **95.7%** | — | — | **NOT charge/mono** (4.3% wrong-charge) |
| peak count (median) | **500** | **500** | — | **NOT peak count** (both at the cap) |
| b/y ladder coverage (median) | **0.083** | **0.278** | 0.039 | **fragment PURITY** — 3.3× fewer true ions |
| ladder cov on SHARED | 0.401 | 0.462 | — | our quality COLLAPSES at low abundance |
| missing-ion recovery from our ±2-cycle neighbours | 0.095 | — | — | modest per-spectrum |
| misses liftable to ≥the reference implementation cov by completion | **990 (~35%)** | — | — | adjacent completion is a real lever |

**Verdict.** The recall gap is NOT precursor annotation (charge is right 95.7%) and NOT peak count
(both capped at 500). It is **fragment purity/selectivity**: at low abundance our 500 emitted peaks
are diluted with noise / co-isolation fragments, matching only 8.3% of the true b/y ladder vs
the reference implementation's 27.8% — barely above the 3.9% chance floor, so below Sage's min-matched-fragments. Our
per-spectrum quality collapses from 0.401 (shared) to 0.083 (missed) while the reference implementation degrades
gracefully (0.462→0.278). Separately, ~35% of missed peptides could reach the reference implementation-level ladder
coverage by completing the ladder from our OWN adjacent-cycle spectra — a completion lever, distinct
from the falsified coordinate-merge (that fused unrelated spectra; this gathers the SAME precursor's
fragments across cycles to fill ladder gaps).

---

## ADVERSARIAL REVIEW of the cause-separation (codex + vibe + kimi, 2026-07-26)

**The diagnosis DIRECTION survives, the MAGNITUDE is soft.** The SHARED positive control (both high,
same method) rules out a pure artifact, but the 3.3× "purity collapse" is confounded three ways
(all three reviewers):
- **Asymmetric spectrum selection** (codex, kimi): the reference implementation's side uses its exact identified scan;
  our side uses `near()`+first-match — not the best member of our multi-spectrum family — deflating
  our number. Fix: peptide-blind best-member selector + an oracle upper bound; also run the reciprocal
  ours-only set.
- **dt decoy floor never computed** (kimi): the 0.039 chance floor is ours only; if the reference implementation's 500
  peaks are denser (it pre-deisotopes), its floor is ~0.10 and the gap shrinks to ~1.8×.
- **Coverage over 500, not Sage's top-150-by-intensity** (kimi, codex): our correlation-ranked cap
  keeps faint true fragments Sage never scores, so 0.083 over-states our identifiability — the
  effective gap is likely *worse*. Recompute on the intensity top-150, deisotoped.
- The +13 Da null is fair-ish but measures specificity not purity, and understates co-isolation
  (pure-chance term is ~0.013; the 0.039 excess is real other-peptide ions).

**Convergent mechanism (all three): the correlation gate becomes a LOTTERY at low S/N.** Two noise
processes compound: (i) Pearson r between noisy profiles attenuates, dropping true faint fragments
below `min_correlation=0.3`; (ii) the zero-padded full-grid Pearson (`gate:coelution=pearson`
default) + tiny overlap counts (3–6 points) let noise/co-isolated fragments pass by chance (|r|≥0.3
at ~1/3 rate for n=4). The 500-slot cap then fills with lottery winners / co-isolated interferents,
not the true ladder. This is the ONLY mechanism inherently abundance-dependent in the right
direction, and it predicts the observed 0.401→0.083 collapse (the reference implementation smooths its precursor and
deisotopes → degrades gracefully 0.462→0.278). NOT charge, NOT peak count, NOT threshold strictness.

**Wavelet smoothing is NOT the fix** (codex + kimi refute vibe's "enable it"): at 2.6 samples/FWHM,
0.5×FWHM ≈ 1.3 samples is below the à-trous first meaningful scale, so `atrousLevels()` returns 0
(no-op) or smears the apex and *reduces* true r; and it smooths only the PRECURSOR (raises rel_p) —
the faint FRAGMENT side still attenuates, so it is at most half a fix. Measure the pass-rate-vs-
intensity curve before spending here.

### Ranked fixes (converged)
1. **MS2 deisotoping before the 500 cap** — cluster 1.0034/z-spaced fragment peaks, keep the
   monoisotope (sum intensity), THEN cap. Composition-only: zero new peaks/spectra, **FDR-neutral by
   construction**, frees ~150–200 slots (isotopes are ~30–40% of the cap, preferentially from
   abundant interferents at low S/N). NEW feature (none exists). Validate OFFLINE first: envelope
   count + envelope-coverage ours vs the reference implementation on the MISS set (no re-search); if ours rises to ≥70%
   of dt's, implement. **Highest confidence-to-effort; do first.**
2. **Fix the gate STATISTIC, not the threshold** — enable/retune the already-implemented
   `gate:coelution=logoverlap` (AlphaDIA-style, summed only over precursor support), and/or add a
   min-overlap-points floor that RISES as precursor S/N falls (kills the 3–5-point lottery directly).
   Highest ceiling, but NOT FDR-neutral → must pass entrapment (bench/entrapment.py) on both engines.
3. **Free A/B: `assembly:rank_by=intensity`** — one flag, already implemented; Sage keeps top-150 by
   intensity so a correlation-ranked cap can discard exactly the peaks the engine scores. Run same day.
4. **The 13% NOT-emitted (~390 peptides)** — give it its own extraction cause-split (trace too short?
   ms1_split_valleys merge? guessed-mono failure?); an extraction fix has no FDR surface.

### Verdict on abundance-adaptive extraction (the author's direction) — half right, half backwards
- **m/z relaxation for faint precursors: NO.** m/z tolerance barely enters fragment *selection* (the
  gates are IM/RT/correlation); relaxing it in tracing/monoisotope buys merged traces + wrong
  monoisotopes, and admits co-isolated fragments the (already-lotteried) gate can't reject → Sage
  mis-assigns them with CONFIDENT scores (real ladder, wrong peptide) → peptide-FDR inflation.
- **"Weight co-elution more for faint": backwards while the statistic is the degraded one** — at low
  S/N the Pearson is what failed; leaning on it imports noise. Fix the statistic FIRST (fix 2), then
  weighting it is meaningful.
- **IM-consistency: directionally sound, but uncertainty-normalized.** The FRAGMENT IM apex is also
  noisy at low S/N; tighten `gate:delta_im` only to ~0.005–0.007 for faint precursors (below that
  rejects true faint fragments on measurement error). Best done as a CONJUNCTIVE trade (lower min_corr
  only in exchange for a tighter IM gate + raised overlap-points floor), never a unilateral relaxation.

### The trap (all three): adjacent-cycle ladder completion
It optimizes the exact proxy this analysis measures; the 9.5% recovery sits within striking distance
of its own unmeasured chance floor (~5%/ion); at a binding 500/500 cap completed ions must DISPLACE
real peaks; and the imported noise is what a hyperscore engine converts to FDR, not recall — the
same family as the falsified coordinate-merge (−8% PSMs/spectrum). Gate behind a decoy-completion
control before building.

---

## Explored: SUBTRACTING adjacent spectra (remove shared interferents) — REFUTED

Hypothesis: peaks shared with same-cycle co-eluting neighbours are interferents; subtracting them
would let the faint true ladder rise into Sage's intensity top-150. `subtract_test.py` (1,200 MISS
peptides) refutes it — neighbour-sharing marks REAL fragments, not interference:

| peak class | median sharing | mean | fraction shared (≥2 neighbours) |
|---|---|---|---|
| TRUE b/y ions | 1.00 | **1.63** | **42.7%** |
| OTHER (interferent) | 1.00 | 0.98 | **24.4%** |

Top-150 coverage: raw 0.104 → after subtract(≥2) **0.067** (below the random-removal control 0.100),
stripping 8.3% of true ions. A real fragment co-elutes with its precursor AND recurs across
overlapping isolation windows (the co-elution signal DIA correlation methods exploit); noise is
sporadic. So "appears in neighbours" ≈ real — subtracting it removes signal. This is the mirror of
why adjacent-cycle COMPLETION is a trap: both fight the grain of co-elution. The interferents that
bury our true fragments are intense-but-sporadic, not shared; the lever is promoting faint-but-real
fragments (deisotoping / a cleaner gate), not subtracting shared ones.

---

## Fix benchmark: MS2 deisotoping FALSIFIED as a recall lever (dataset D, Sage)

| arm | spectra | peptides | vs gated baseline |
|---|---|---|---|
| gated baseline | 924,255 | 9,989 | — |
| + `assembly:ms2_deisotope` | 924,233 | **9,856** | **−1.3%** |
| + deisotope + `rank_by=intensity` | 924,233 | 9,852 | −1.4% |
| the reference implementation | 700,434 | 11,517 | — |

Deisotoping (the review's #1, "FDR-neutral") did **not** recover recall — it slightly hurt. This
CONFIRMS the mechanism and falsifies the isotope-dilution hypothesis: freeing the ~22% isotope cap
slots does not promote our faint true fragments into Sage's intensity top-150, because the peaks
burying them are intense **non-isotope co-isolation** fragments. Deisotoping only pays the
false-positive cost (a few true fragments at coincidental 1.0034/z spacing) with no offsetting gain.
`rank_by=intensity` is likewise flat (Sage re-ranks to top-150 by intensity regardless of our cap).
Reverted. **The remaining lever is the one the reviewers ranked #2 (higher ceiling, not FDR-neutral):
fix the co-elution GATE STATISTIC so the intense interferents are rejected upstream** — enable/retune
`gate:coelution=logoverlap` (already implemented) and/or add a min-overlap-points floor that rises as
precursor S/N falls, entrapment-validated on both engines. That directly attacks "true fragments
faint relative to interferents" at its source rather than reshuffling a cap already full of them.

---

## CORRECTION: subtraction refutation was too strong; apex-assignment is untested, not refuted

The `subtract_test.py` refutation tested peak *presence* in neighbours (sharing), which conflates a
real fragment's OWN elution (spans ~2.6 cycles) with interferent bleed — exactly the confound the
tool author flagged. The sharper claim is APEX ASSIGNMENT: a fragment belongs to the spectrum at its
elution apex; where it appears at a different RT (a co-isolated neighbour off its apex) it is an
interferent and should be excluded THERE.

`apex_test.py` tried to test this offline and FAILED (degenerate: 0.2% "aligned"). The failure is the
TEST, not the idea: you cannot reconstruct a fragment's XIC apex from the EMITTED pseudo-spectra —
with share-all assignment (default) a fragment's intensity is copied at its constant trace value into
every precursor that claims it, so across neighbours the peak is flat and has no apex. Retracted.

**Apex assignment is testable only IN THE TOOL** (it holds each fragment's real trace apex `f.rt`).
The tool already gates `fabs(f.rt - pc.rt) > delta_rt` at delta_rt=3.0 s (±2.2 cycles) and then
SHARES the fragment across every co-eluting precursor. The author's refinement = assign each fragment
to the precursor NEAREST its own apex and exclude it from off-apex neighbours (tighten the apex gate,
or apex-competitive assignment = winner-by-apex-RT, distinct from the falsified winner-by-correlation
`assembly:competitive`). This attacks the off-apex-bleed mechanism directly; it needs its own
benchmark (competitive-family knobs have lost peptides before, keyed on correlation not apex).

---

## Apex-competitive assignment FALSIFIED (dataset D, Sage) — and the robust lesson

| arm | spectra | peptides | vs gated baseline |
|---|---|---|---|
| gated baseline | 924,255 | 9,989 | — |
| apex-competitive rp=1 (pure apex-winner) | 634,963 | 3,416 | **−66%** |
| apex-competitive rp=3 | 828,958 | 6,597 | **−34%** |
| the reference implementation | 700,434 | 11,517 | — |

The author's apex-assignment idea (assign each fragment to the precursor at its own elution apex,
drop off-apex bleed) is sound in principle but the removal cure is worse than the disease. The
fan-out histogram shows fragments are MASSIVELY shared — 27.9% span 11-25 co-eluting precursors,
14.7% span 26-50, 7.5% span 51+ (wide diaPASEF isolation window; many peptides genuinely share
fragment m/z). Restricting each fragment to its 1 (or 3) nearest-apex precursor SCATTERS every
peptide's b/y ladder across dozens of spectra, so no single spectrum retains a complete ladder and
recall collapses. Fragment apexes are also noisy at low S/N, worsening the scatter.

**Robust conclusion (3 independent falsifications): the recall lever cannot be REMOVING fragments
from spectra.** correlation-competitive (falsified earlier, assembly:competitive), sharing-based
subtraction (refuted), and apex-competitive (falsified, both rp) all fail the same way — depleting
or scattering the true ladder. The only remaining direction is to change HOW co-elution is SCORED
(reject noise at scoring time via a better statistic) WITHOUT stripping shared true fragments:
gate:coelution=logoverlap (no zero-padding) + a min-overlap-points floor that rises as precursor S/N
falls. That is the next (and last untried) lever.

---

## FINAL: greedy deconvolution + logoverlap falsified — the recall gap is not a fragment-level lever

| lever | peptides | vs gated baseline (9,989) |
|---|---|---|
| **gated baseline (share-all)** | **9,989** | — |
| MS2 deisotope | 9,856 | −1.3% |
| logoverlap mc=0.10 / 0.20 / 0.30 | 9,071 / 8,937 / 8,962 | −9 to −11% |
| apex-competitive rp=3 / rp=1 | 6,597 / 3,416 | −34 / −66% |
| greedy peel-bright-first deconvolution | 6,505 | −35% |
| the reference implementation | 11,517 | +15% |

**Seven independent approaches now falsified** (deisotope, adjacent completion, sharing-subtraction,
correlation-competitive, apex-competitive, logoverlap gate-statistic, greedy deconvolution). They all
fail for ONE reason, visible in the fan-out histogram: fragments are shared across 11-50 co-eluting
precursors (wide diaPASEF window; many peptides genuinely share fragment m/z). ANY scheme that
redistributes, reduces, removes, or re-scores fragment intensity by precursor assignment DEPLETES the
faint precursors — they are always "last in line" behind dozens of brighter legitimate co-eluters, so
their residual intensity/fragment count collapses and they stop identifying. The greedy deconvolution
made this explicit: by the time the residual pool reaches a faint precursor, brighter co-eluters that
each correlated (even weakly) have already claimed most of it. **Share-all at full intensity is
robustly optimal**; every attempt to be cleverer about assignment loses the low-abundance tail it was
meant to save.

**Conclusion.** The recall gap (~1,500 peptides, the low-abundance tail, ~87% of the reference implementation on dataset D) is
NOT closable at the fragment-assignment / co-elution-scoring / intensity-deconvolution level. It is
structural to a wide-isolation-window pseudo-spectrum built by per-precursor fragment correlation.
Closing it would require what the reference implementation does UPSTREAM — cleaner MS1 feature/precursor definition and
precursor-profile smoothing (2D Gaussian + Savitzky-Golay) so faint precursors are defined and scored
more sensitively BEFORE fragment assignment — a feature-detection re-engineering effort, not a
scoring tweak. The committed isotope-support gate (9,989 peptides, near-neutral, both engines) stands
as the recall/quality operating point.

---

## FIRST RECALL GAIN: IM continuous weighting (+4.1%, dataset D/Sage)

The adversarial review's untried lever WORKS. `assembly:im_weight_sigma` down-weights each emitted
fragment's intensity by exp(-(fragIM-precIM)^2 / 2 sigma^2). TIMS precedes fragmentation so a true
fragment shares the precursor 1/K0 exactly; off-IM co-isolation interferents are suppressed, so faint
true fragments rise into Sage's intensity top-150. WITHIN-spectrum (no fragment removed/redistributed
across precursors) -> escapes the fan-out depletion that falsified all 7 prior levers.

| arm | spectra | peptides | vs gated baseline (9,989) |
|---|---|---|---|
| IM weight sigma=0.005 | 924,255 | **10,394** | **+4.1% (+405)** |
| IM weight sigma=0.003 | 924,255 | 10,013 | +0.2% (too aggressive) |
| the reference implementation | 700,434 | 11,517 | — |

Closes ~27% of the the reference implementation gap at ZERO emission cost (same 924k spectra; only intensities
reweighted). sigma=0.005 optimal so far. Validation pending: entrapment FDR (is it real?), MSFragger
(both-engines), finer sigma. Also refutes the search-side hypothesis for OUR benchmark (chimera
-2.2%/-11.7%; and we searched the reference implementation's spectra through the SAME bare Sage, so the 9,989 vs 11,517
gap is extraction, not a FragPipe/rescoring confound).

---

## FULL the reference implementation COMPARISON, BOTH ENGINES (dataset D) — IM-weighting is Sage-specific

| | Sage (peptide_q<=0.01) | MSFragger (TD 1% FDR) |
|---|---|---|
| baseline-gate | 9,989 | 11,033 |
| ours-best (IM sigma=0.005) | 10,394 (+4.1%) | 11,060 (+0.2%) |
| the reference implementation | 11,517 | 13,932 |
| ours / the reference implementation | 90% | 79% |

Entrapment (Sage): baseline target=9,304 @2.33% -> imw005 target=9,798 @2.26% (+5.3% target at LOWER
true FDR -> the Sage gain is genuine signal). BUT the +4.1% is SAGE-SPECIFIC: flat (+0.2%) on MSFragger
(Sage keeps top-150 by intensity so the IM reweight rescues faint fragments; MSFragger's scoring is
already robust to them). And the the reference implementation gap is WIDER on MSFragger (79%) than Sage (90%) -- the reference implementation's
spectra carry structure MSFragger exploits that ours lack. So IM-weighting is a real but engine-specific
partial win; the core spectrum-content gap persists and is worse on the reference implementation's native engine.
