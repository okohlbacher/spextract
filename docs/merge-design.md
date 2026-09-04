# Merging split pseudo-spectra: final design and verdict

**Date:** 2026-07-22. **Status:** across-cycle KILLED; within-cycle CONDITIONALLY APPROVED, default
OFF, ship blocked on a pre-registered paired measurement.

This supersedes the proposal in `merge_design.md`. It records only what survived four adversarial
reviews and one direct measurement on dataset B. Several numbers in `docs/benchmarks/split-spectra.md`
are corrected here; that file must be amended (§6).

---

## 0. What the merge is actually for

The proposal justified merging by **spectral richness** — union recovers peaks that splitting
scattered. The measurement does not support that as the primary benefit (§1). What survives is a
different and better-supported claim:

> `consolidate:delta_rt` already collapses within-cycle multiplicity **by deletion**. A union merge
> achieves the same count collapse **without discarding the non-shared peaks**. Merge is
> consolidate done without the deletion.

The value proposition is therefore *count collapse at no peak cost*, not *richer spectra*. Whether
count collapse itself helps anything is **unknown** — six redundancy-reducing levers have been
falsified in this project, and the `consolidate:delta_rt` sweep is the floor for what merge can
achieve on peptide count. State it that way; do not re-inflate the richness claim.

---

## 1. VERDICT — within-cycle merge: **MODIFY AND GATE** (implement in restricted form, default OFF)

### What holds

Direct measurement on dataset B (76,119 PSMs at 1% PSM-FDR, 10,053 peptides; pairs filtered by the merge
predicate; four null arms carried) shows the within-cycle population has the strongest co-elution
signature of anything tested:

| | Jaccard | cosine | bitwise-identical shared peaks |
|---|---|---|---|
| within-cycle (ΔRT 0.000 s) | 0.456 | 0.672 | **0.747** |
| same peptide, 45.7 s apart | 0.098 | 0.194 | 0.000 |
| different peptide, matched | 0.046 | 0.040 | 0.000 |
| random spectrum | 0.020 | 0.009 | 0.000 |

Merging (union, cap re-applied at 500) gains **+1.23 annotated b/y peaks** per pair against a
matched-different-peptide null of **−2.20**. The sign is right and the null separation is large.

The step is also architecturally safe in a way `trace:frame_aggregation_n` was not. That lever ran
`aggregateFrames_` (src:218–285) at src:1312 and src:1498, **before** `detectTraces_` at src:1321
and src:1500 — it corrupted the XICs feeding the Pearson gate, and its own doc comment (src:213–214)
says adjacent output points become autocorrelated, "which inflates Pearson r for signal and noise
alike." A merge over `all_out` runs at src:1727, after `storeExperiment` is the only consumer
(verified: `all_out`/`out_exp` appear nowhere outside src:1683–1781). It cannot alter any inferred
quantity. Every merged peak already passed `gate:min_correlation` inside `scoreCandidates_`
(src:1015–1103, gate at src:1100). This distinction is real, not rebranding.

### What does not hold, and why the verdict is MODIFY rather than IMPLEMENT

**The headline evidence for "split spectra" is partly an estimator artifact.** `Δ(RT) = 0.000 s` was
called "decisive." It is quantisation. With `trace:ms1_split_valleys` on (default 7.0, src:517),
`ElutionPeakDetection` calls `updateSmoothedMaxRT()` at both `ElutionPeakDetection.cpp:452` and
`:535`, and `MassTrace.cpp:507` sets `centroid_rt_ = trace_peaks_[max_idx].getRT()` — the RT of an
actual MS1 frame. Confirmed by measurement: 1,032 distinct emitted RT values, consecutive gap median
**1.3850 s** (p25 1.3850, p75 1.3851) — exactly the MS1 frame grid. Two *unrelated* peptides sharing
a frame also read ΔRT = 0.000. The statistic carries no provenance information. The base arm's
0.366 s is a different estimator (`updateWeightedMeanRT()`, `MassTraceDetection.cpp:664`), not
different biology.

**`union/larger` was read with its poles inverted.** `split-spectra.md` states "≫1 → split; ~1 →
distinct." Measured floors: *distinct* peptides give **1.88–1.95**; identical spectra give 1.0. The
published **1.355 sits at the similar end of the scale**, so it is evidence of *similarity*, not of
complementarity. Risk 4 in the proposal ("base has HIGHER union/larger 1.450 than treat 1.355, so
splitting exists without ms1_split_valleys too") reads backwards: higher means *less* alike.

**"41.1% of peaks unique to one" is mislabelled and near its floor.** `bench/split_test.py:115`
computes `1 − inter/min(|A|,|B|)` and prints it at `:134` as "unique to one." Measured floors for
that statistic: within-cycle 0.359, different-peptide null 0.908, random null 0.958. The number is
near its *minimum*. Drop it.

**The benefit is small and mostly destroyed by the fragment cap.** Union gains +2.50 annotated peaks
pre-cap; re-applying `assembly:max_fragments = 500` (src:494) leaves **+1.23** — the cap destroys
83% of the gain. Only **52.4%** of pairs improve. The gain is ~6% of the annotated peaks already in
the best member. And raw purity *falls* on union (b/y fraction 0.0373 → 0.0351); the apparent
recovery to 0.0438 in the capped column is the intensity cap re-selecting, not the merge.

**Free finding that changes the costing.** An unbiased 4,000-spectrum sample of the emitted file:
**median 500 peaks, mean 435.8, 75.7% exactly at the 500 cap, 0% empty** (86.3% at cap among
PSM-bearing spectra). This resolves the repository's 10× contradiction: `harness/collate.py:124`
("we ~228") is 2× low and `scripts/truthset.py:9` ("we emit ~23") is **20× low**. Both must be
corrected. Consequence: a union merge produces 700–900 peaks that the cap immediately truncates back
to 500, so the merge **cannot deliver its union without also raising the cap** — and raising the cap
raises the open-search chance-match rate, which grows as n⁴ under `min_matched_peaks: 4`
(`bench/sage_open.json`). At n = 500 that cost is already the dominant term. **Do not raise the cap
as part of this change.**

### Verdict

Implement the within-cycle merge in the restricted form in §3, **default OFF**, as a pure
post-processing pass. Ship only if the §4 metric fires, and only if it beats the consolidate null.
Justify it by count collapse, not by richness.

---

## 2. VERDICT — across-cycle merge: **KILL**

Not "kill because it adds interference." **Kill because the claimed mechanism is absent from the
data.** Measured on dataset B, all pairs filtered by the proposed merge predicate (same charge, Δm/z ≤ 20
ppm, ΔIM ≤ 0.01), medians:

| | n | Jaccard | union/larger | cosine | bitwise | cap gain (annotated b/y) | % pairs improved |
|---|---|---|---|---|---|---|---|
| lag 0 (0.00 s) | 2000 | 0.456 | 1.344 | 0.672 | 0.747 | +1.23 | 52.4 |
| lag 1 (1.39 s) | 2000 | 0.346 | 1.454 | 0.590 | 0.629 | +1.03 | 51.7 |
| lag 2 (2.77 s) | 2000 | 0.236 | 1.579 | 0.465 | 0.463 | +1.26 | 57.0 |
| lag 3 (4.16 s) | 782 | 0.156 | 1.684 | 0.350 | 0.239 | **+1.35** | **62.1** |
| same peptide, 45.7 s apart | 210 | 0.098 | 1.802 | 0.194 | 0.000 | +0.49 | 45.2 |
| different peptide, matched | 1895 | 0.046 | 1.878 | 0.040 | −2.20 | −2.20 | 10.9 |

**Four independent reasons, any one sufficient:**

1. **No window can be read off the data.** Every co-elution statistic decays smoothly to its null
   with no shoulder: bitwise 0.747 → 0.629 → 0.463 → 0.239 → 0.000; Jaccard 0.456 → 0.346 → 0.236 →
   0.156; cosine 0.672 → 0.590 → 0.465 → 0.350. There is no breakpoint separating "adjacent" from
   "not adjacent." Any `rt_window` would be arbitrary. This answers proposal question 1: **nothing
   bounds the window.**

2. **The benefit is flat or *increasing* in lag** — the opposite of the design's prediction. Cap
   gain +1.23 / +1.03 / +1.26 / **+1.35**, and 62.1% of lag-3 pairs improve versus 52.4% at lag 0.
   Mechanism: more separation → more distinct peaks → more candidates for the intensity cap to pick
   from. The "benefit" is a selection effect of having more peaks, not of reassembling a split
   feature.

3. **The benefit tracks peptide identity, not chromatography.** At 45.7 s apart — where the two
   spectra provably cannot be one chromatographic feature — merging still gains +0.49 annotated
   peaks, 45.2% of pairs improve, and added peaks annotate at 0.016 against a 0.005 decoy floor. So
   "the neighbour contributes real fragment ions" is true of *any* second spectrum of the same
   peptide anywhere in the run. Accepting that as justification commits you to merging the entire
   run. Nobody proposes that.

4. **The neighbour is largely a copy, not an observation.** 62.9% of shared peaks at lag 1 are
   bitwise-identical in m/z *and* intensity. `Trace::intensity = mt.getMaxIntensity(false)`
   (src:175) is a whole-trace scalar, `frag_traces` is detected **once per isolation window over the
   whole window map** (src:1500), and share-all is the default path (`assembly:apportion` 0.0,
   `assembly:rp_max` 0, src:496–497), so `frags.emplace_back(frag_traces[fi].mz,
   frag_traces[fi].intensity)` (src:1154) emits the *same double* into every gated precursor. Merged
   intensities therefore carry **no new information at any lag**, and summing double-counts one
   trace's maximum.

Supporting arguments, not needed for the verdict but consistent with it: post-assembly merging
cannot recover a fragment that failed *detection*, only one that failed *correlation* — the gate
whose job is rejecting interference (src:1100). It is frame aggregation minus its only real benefit,
the √N pre-detection gain documented at src:205–209. And an `rt_window` growing toward the ~25 s
scale that `trace:ms1_split_valleys` exists to separate (src:513) re-fuses what that lever paid
+36.7%/+37.4% peptides to split.

**Do not re-propose across-cycle merging without a new mechanism and a new measurement.** The
existing measurement is a negative result, not an absence of evidence.

---

## 3. FINAL DESIGN — only what survived

### 3.1 Form: a standalone post-processing pass

Implement as a pass over the emitted mzML, invokable on an existing file, **not** as an inline stage
of the tool run. Three reasons:

- The merge is then a **pure function of a fixed input**, so merge-OFF vs merge-ON is exactly paired
  with **zero run-to-run variance**. This project owns no replicate data and therefore no noise
  floor (§4.4); this construction removes the need for one entirely.
- It cannot touch any inferred quantity, which is the whole basis of the frame-aggregation
  distinction (§1).
- One emitted mzML yields three searchable arms (OFF / merge / consolidate) at the cost of one
  pipeline run.

The `merge` block runs at src:1727, in place of or after the `consolidate` block, on the
already-canonically-sorted `all_out`. **Run once. It is not idempotent** — merging changes spectrum
content, which changes the canonical sort key (src:1693–1709 tiebreaks on fragment m/z and
intensity), which changes the partition on a second pass. Do not iterate; do not re-sort and re-run.
The merge must be a pure function of the sorted `all_out` content and must never depend on trace
indices (trace order is OMP-nondeterministic when `ms1_split_valleys > 0` —
`ElutionPeakDetection.cpp:320–341` writes `single_mtraces` in thread-completion order under a
critical section).

### 3.2 Parameters

```
merge:within_cycle      flag           default OFF     enable the within-cycle merge
merge:mass_ppm          <ppm>   20.0   precursor m/z tolerance
merge:delta_im          <1/K0>  0.01   precursor ion-mobility tolerance
```

That is the complete parameter set. Note what is **absent**:

- **No `rt_window`.** Killed by §2. Grouping is exact-frame only (§3.3).
- **No `mz_tol`** for peak matching. Peaks are deduplicated on **exact m/z equality**, never a
  tolerance (§3.4).
- **No `same_charge_only` flag.** Same charge is **mandatory**, not optional (§3.3).
- **No apex/TIC selector.** Killed by §3.5.

`merge:mass_ppm` default 20.0 matches the established `consolidate:mass_ppm` (src:526) and is
measured-safe: the within-cycle isotope-misassignment population sits at ~1000 ppm ≈ one isotope
spacing (1.00335/2 = 0.5017 Da = **836 ppm at m/z 600**), 40× above the tolerance. **Hard upper
bound 100 ppm** — above that, merge fuses a monoisotopic assignment with a +1-isotope
misassignment, i.e. two different precursor mass hypotheses in one spectrum. The tool must refuse
values > 100.

`merge:delta_im` default 0.01 matches `gate:delta_im` (src:487). Note the control ΔIM for the target
population is 0.0052 — 52% of the budget — and the false-capture population's ΔIM has never been
measured. Do not loosen this.

### 3.3 Grouping rule

Within `all_out` sorted by (RT, precursor m/z, charge, content):

```
group i and j iff
    all_out[i].getRT() == all_out[j].getRT()                       // EXACT equality
and |mz_j - mz_i| <= mz_i * merge:mass_ppm * 1e-6
and  charge_j == charge_i                                          // mandatory
and  im_i and im_j both present  and  |im_j - im_i| <= merge:delta_im
```

Four points, each fixing a confirmed defect:

**RT equality is exact, not windowed.** Emitted RT is quantised to the MS1 frame grid (measured
1.3850 s spacing, §1), so "same cycle" *is* bitwise RT equality. A tolerance would only admit
adjacent cycles, which §2 kills. Note the corollary: ΔRT = 0 does *not* by itself imply shared
provenance — unrelated peptides share frames too. The m/z, IM and charge predicates do all the
discriminating.

**The IM gate must read the precursor, and must fail closed.** `consolidate`'s IM gate **has never
fired.** src:1745 and src:1755 call `all_out[i].getDriftTime()` — `MSSpectrum::drift_time_`, which
is never written anywhere in the tool and defaults to −1 (`MSSpectrum.h:693`). The only
`setDriftTime` call is src:1137, `prec.setDriftTime(pc.im)` on the `OpenMS::Precursor`. The guard
`im_i > 0 && im_j > 0` then short-circuits on every pair. Every consolidate sweep run to date merged
on (RT, m/z, charge) only, and the parameter's own help text (src:527) describes exactly the
IM-sub-range merge it cannot perform. `merge` must read `getPrecursors()[0].getDriftTime()`, and
must **decline the merge** when IM is unavailable rather than skipping the check. Fix the
`consolidate` block the same way in the same commit.

**Same charge is mandatory.** src:955–957 records the reference implementation retaining multiple charges when
ambiguous; merging deletes that hypothesis. `inferPrecursors_` marks partners `used[]` and emits one
precursor per seed trace, so two charge-disagreeing emissions come from genuinely different seed
traces — different inferences, not one feature. The ~6.8% charge-discordant within-cycle pairs stay
unmerged **by design**, and the pass must **report their count**: they are enriched post-merge by
construction, which is a confound for any downstream fraction-of-PSMs statistic (§5.10).

**Grouping is seed-anchored and will under-merge.** `mz_i`/`im_i` are fixed on the seed and never
updated (as at src:1743–1745). With 20 ppm windows, 7.5 ppm median pairwise separation and mean
multiplicity 2.85, groups of 4+ members will fragment. This is the **conservative** direction and is
accepted deliberately. The pass must **report the number of groups whose last accepted member sat
within 10% of the window edge**, so under-merging is measured rather than assumed. Do not "fix" it
with transitive closure — a ppm tolerance chains without bound.

Also fix, while in this code: src:1749 reads `all_out[j].getRT()` before the `dropped[j]` guard at
src:1750, and `best` may exceed `i` and be moved at src:1766. Benign today (`retention_time_` is a
plain double surviving the move) but it is a moved-from read.

### 3.4 Merge rule

```
merged peaks = union of member peaks, deduplicated on EXACT m/z equality,
               each surviving peak keeping its own intensity (no combination)
then truncate to assembly:max_fragments by descending intensity
```

**Intensities are never summed and never averaged.** Two members that share a fragment emit the
*same double* — `Trace::intensity` is `mt.getMaxIntensity(false)` (src:175), `frag_traces` is
detected once per isolation window (src:1500), and share-all is the default path so no per-precursor
scaling occurs (src:496–497, src:1154). Measured: **74.7%** of shared within-cycle peaks are
bitwise-identical in m/z *and* intensity. Summing multiplies one trace's maximum by its membership
count, and membership correlates with fragment intensity and fan-out, so it systematically re-ranks
the spectrum toward high-fan-out (chimera-prone) fragments before truncation. Union without
combination is exact, not a compromise.

**Exact m/z equality, no tolerance.** Members of one within-cycle group draw from the identical
`frag_traces` array, so genuinely shared peaks *are* bitwise identical. Two peaks that differ at all
are two traces that `MassTraceDetection` deliberately separated; do not second-guess it. The single
case where a tolerance would be defensible — a precursor on a shared isolation-window boundary,
assembled from two independently detected `frag_traces` arrays (src:1536–1537 uses `lower_bound` /
`upper_bound`, inclusive at both ends) — is rare and is left unmerged.

**The cap stays at 500, applied once after merge, by descending intensity.** This knowingly discards
83% of the pre-cap peak gain and is the right trade: spectra are already at the cap 75.7% of the
time, and n⁴ chance matching in open search makes any increase in n expensive. Two honest
consequences to record: (a) the correlation-based rank score from `assembleFromList_`
(src:1113–1126) is **not comparable across members** — a shared peak has a different correlation in
each — so intensity is used instead, and post-hoc merging therefore has no principled ranking; (b)
the principled fix is to group precursor *hypotheses* **before** `assembleFromList_` and assemble
once from the union of gated fragments with a well-defined per-fragment score. That is the correct
implementation. The post-hoc pass is the cheap, low-risk, exactly-paired version that should be
measured first. Raising the cap is a **separate experiment with a separate justification**; it is
not part of this change.

### 3.5 Apex / representative rule

Within-cycle members share a frame, so **RT is identical by construction and there is nothing to
choose.** Merged RT = the common RT.

For the merged precursor m/z, charge and IM: take them from the member whose **MS1 precursor trace
apex intensity** is highest. This requires stamping that scalar at assembly time — `dedupPrecursors_`
already uses exactly this quantity (`ms1[pcs[j].trace_idx].intensity`, src:1003) — as a metavalue on
the emitted spectrum, since `trace_idx` is not otherwise recoverable post-hoc. One line at src:1129.

Both criteria the proposal offered are rejected:

- **Highest TIC is a selection effect.** TIC is a function of how many fragments passed the gate.
  Selecting the representative by the fragment content that is then scored contaminates every
  per-spectrum quality number with a selection no decoy floor controls. `consolidate` already does
  this wrong (src:1758–1762 tiebreaks on `calculateTIC()`); this is a correction to inherited
  behaviour, not a new choice.
- **Highest MS1 trace intensity cannot discriminate *cycles*.** It is a whole-trace apex scalar
  (src:175), constant across every cycle the trace is emitted in. It works here only because
  within-cycle members are *different traces*. It would have been blind for the across-cycle merge —
  a further reason that half was unimplementable as specified.

Corollary worth stating plainly: because merged peak intensities are cycle-invariant by
construction, a merge can only ever add **fragment identities**, never intensity information. That
is a much smaller claim than the proposal made.

---

## 4. METRIC — how this is judged, and what kills it

### 4.1 Why the obvious metrics are all disqualified

| candidate | why it cannot judge this change |
|---|---|
| closed peptide count alone | `peptide_q` makes redundancy free score. It rewards the thing merge removes. It is reported, never gating. |
| "peptides flat while PSMs/peptide → 1.69" (the proposal's candidate) | **Circular and non-discriminating.** `sage_closed.json` and `sage_open.json` both set `report_psms: 1` and `chimera: false`, so one PSM per spectrum makes PSMs/peptide → 1 an *arithmetic identity* of collapsing n spectra into 1. Worse, it falls **identically** whether the extra PSMs were merged or deleted — so it cannot certify the design's one distinguishing claim against the consolidate null. |
| `joint_bench.py` PSMs/peptide | Mixed acceptance sets: `joint_bench.py:176–183` runs `target_decoy_fdr` separately at PSM and peptide level, so the numerator's PSMs are not drawn from the denominator's peptides. Single-set values (`cycle_redundancy.py:57–77`) are 7.57 / 1.62, not 9.42 / 1.71 — the published 5.6× excess overstates by ~18%. The correct scale is **4.67×**. |
| `truthset.py` RECALL | **A correct merge reduces it mechanically.** `truthset.py:124–127` credits one precursor per spectrum (`best = min(cand, …)`), so two members crediting two different DIA-NN precursors become one; and RT reassignment can push a spectrum out of the stratified `in_sel(rt)` bins. Separately its decoy is partly void: when `pmz + 11.003 > hi` no candidate can match (`:114–116` vs `:121–123`), voiding 34.4% equal-weight / 47.1% over windows 2–15 / 19.6% width-weighted. Report it; never gate on it. |
| `open_bench.py` `ray_fraction` | Moving denominator (`shifted = rays + band + unassigned`, `:148–150`), moved **selectively**: mandatory same-charge merging collapses the 93.2% concordant groups ~2.85→1 while the ~6.8% discordant groups — the population enriched for charge errors — survive unmerged, so the fraction rises with no new error created. Absolute ray counts are equally confounded (total PSMs fall ~2.85×). Quote **rays per accepted peptide** against `ray_fraction_chance` (`:160`), using the counts already exposed at `:231`. |
| per-spectrum means in `truthset.py` | Survivorship. `:137–140` appends one record per matched *spectrum* and `:157–163` reports `mean()`. Merge changes both the size and composition of that population; if it preferentially collapses the sparsest spectra the mean rises with nothing improved. |

### 4.2 The primary instrument: paired per-precursor spectrum quality

Run the pipeline **once**. Apply the merge as a post-processing pass to produce three files from one
emitted mzML:

- **arm A — OFF**: the emitted file unchanged
- **arm B — MERGE**: `merge:within_cycle` on
- **arm C — CONSOLIDATE (the null)**: the identical grouping with deletion instead of union

C is not optional. It is the same grouping predicate with the same tolerances differing only in the
reducer, it is already implemented, and it is the only arm that isolates *union* from *count
collapse*. If B does not beat C, the union half buys nothing and the correct action is to ship C.

Define the population **externally and once**: the set of DIA-NN `Precursor.Id` values from the
reference report. For each precursor, in each arm, take that precursor's single best-matching
spectrum (highest-scoring member in A and C, the merged spectrum in B). Compare **paired** over the
common set. Report at **both** depths — TOP-20 (`TOPN` at `truthset.py:23`) and full 500 — since
both are already emitted at `truthset.py:157`.

Per-precursor quantities:

1. **Annotated b/y fragment count**, at equal depth.
2. **Purity** — b/y fraction — **net of the +11.003 Da decoy-shift floor** (`collate.py:19`), which
   demonstrably has teeth on this data (`ms1-funnel.md:19–30`: floors 8.1% / 24.0% / 19.4% / 48.4%).
3. **Sage `matched_peaks` and score** for the PSM on that spectrum.

### 4.3 Why this is not gameable by the change itself

- **The denominator is fixed externally.** The precursor set comes from the DIA-NN reference, not
  from either arm's own identifications. Merge cannot change the population it is scored over. This
  is the correction README already applied at run level after the spurious dataset D "+4.1%"; the
  instrument has not received it.
- **Equal-depth truncation removes the density advantage.** Merge's most obvious way to look good is
  emitting more peaks; TOP-20 forbids it.
- **The decoy-shift floor removes the chance-annotation advantage**, which grows with n and is the
  same asymmetry that makes `min_matched_peaks: 4` rescue near-threshold targets (a true peptide at
  3 real matches crosses on one chance match, ∝ n; a decoy needs four, ∝ n⁴).
- **Pairing removes the survivorship artefact** that per-spectrum means suffer.
- **Zero variance.** All three arms derive from one emitted mzML by deterministic post-processing,
  and Sage is deterministic. There is no run-to-run noise to hide in — which matters because this
  project owns **no replicate data at all** (§4.4).
- **The consolidate null defeats the remaining confound.** Any effect that count collapse alone
  produces — FDR competition relief, PSM-pool shrinkage — appears in C as well as B. Only the
  B-minus-C difference is attributable to union.

### 4.4 On noise floors — what we do not have

The 0.17% figure (8,411 / 8,408 / 8,406 / 8,397) **is not replicate data.** Those are four different
command lines varying `charge:ambiguity_margin`, a parameter that is structurally inert under the
`charge:scoring = count` default: src:827 opens the `!envelope_scoring` branch and src:861 closes it
with `continue`, while the `ambiguity_margin` block sits at src:957–968 inside the envelope path.
The four numbers are also *monotonically decreasing*, which is not what run-to-run scatter looks
like. **Do not quote 0.17% as a noise floor of any kind.** The only variance estimate this project
owns is s ≈ 1.34 pp between-sample, which at n = 3 gives a ±3.34% non-inferiority interval — too
coarse to judge this change, and inapplicable to a paired within-sample comparison anyway. The §4.2
construction is chosen precisely because it needs neither.

### 4.5 Pre-registered kill criteria

The merge is **killed** — reverted, default stays OFF, and the result recorded — if **any** of:

- **K1.** Paired **decoy-net purity at full depth falls** on the common precursor set (arm B < arm
  A), by any amount outside the paired 95% CI. Rationale: this is the failure mode merge could
  plausibly cause, and the measurement already shows union *lowering* raw purity 0.0373 → 0.0351
  before the cap. This is the criterion most likely to fire; that is why it is first.
- **K2.** Paired **annotated b/y count at TOP-20 does not rise** (B ≤ A). Without an equal-depth
  gain, merge has added nothing but bookkeeping.
- **K3.** Arm **B does not beat arm C** on K1 and K2. Union bought nothing; ship consolidate instead.
- **K4.** **Open-search modified-peptide count at ≤1% entrapment FDR falls** (`open_bench.py:38`).
  Open/blind search is the stated purpose and is the instrument that does not reward redundancy the
  way closed `peptide_q` does. **Prerequisite:** `run_open.sh`'s envelope positive control must be
  fired and *recorded in `docs/benchmarks/`* first — it currently is not (results land in untracked
  `/path/to/scratch/bench/open/results.json`), while "2.82× PSMs for 1.05× modified peptides" is
  already published at `README.md:67` and `joint-2026-07-22.md:33` with no recorded provenance. The
  synthetic `ray_decomposition` validation at `bench/README.md:26` validates the *algorithm*, not
  the instrument on this data. Until that control is on record, K4 cannot be evaluated and the merge
  cannot be shipped.
- **K5.** **Post-merge within-cycle multiplicity ≈ 1.00.** The ~44%-of-pairs isotope-misassignment
  population (~1000 ppm apart) must *not* merge. A multiplicity at 1.00 means it did. The pass must
  emit a **declined-on-mass count**; a value near zero is itself the failure signal. (Caveat: the
  "~44%" is *inferred* from where `d_mz` p75 lands — `split_test.py:131` prints only median/p25/p75.
  Re-run reporting the bimodal split explicitly before pre-registering a number on it.)

**Not a kill criterion:** a fall in closed peptide count. It is reported and it may well fall;
`peptide_q` rewards redundancy, and the `consolidate:delta_rt` sweep already establishes the floor
for what count collapse costs there. A large fall (say >5%) triggers **investigation and a written
explanation**, not automatic reversion — but shipping over a fall requires K1–K4 all clean.

---

## 5. REJECTED — do not re-propose

1. **Across-cycle merge, any `rt_window` > 0.** §2. Four independent measured reasons. The
   population it targets shows no breakpoint, the benefit rises with lag, the benefit is present at
   45.7 s where the premise cannot hold, and the added intensities are copies.
2. **Summing intensities at matching m/z.** Wrong at *every* lag including within-cycle: shared
   peaks are the same trace's `getMaxIntensity()` scalar emitted twice (74.7% bitwise identical
   within-cycle, 62.9% at lag 1). Summing double-counts one number and re-ranks toward high-fan-out
   fragments before truncation.
3. **Averaging intensities.** Same defect, opposite sign — averaging identical numbers is a no-op on
   shared peaks and halves nothing, while making unique peaks and shared peaks incomparable.
4. **Apex assignment by highest TIC.** Selection on the quantity being scored. `consolidate` already
   does this (src:1758–1762) and it is a bug to be fixed, not a precedent.
5. **Using `ms1_traces[trace_idx].intensity` to choose among *cycles*.** It is a whole-trace apex
   scalar, constant across cycles. It is valid only for choosing among different traces in one
   frame, which is the only use §3.5 makes of it.
6. **A `merge:mz_tol` peak-matching tolerance.** Shared peaks are bitwise identical; anything else is
   two traces `MassTraceDetection` separated on purpose.
7. **Reusing the `consolidate` block's grouping code as-is.** Its IM gate has never fired (§3.3), and
   `cons_rt > 0.0` is required to enter the block (src:1728), so the exact-frame configuration the
   evidence supports is *unreachable* through that parameter.
8. **`merge:same_charge_only` as an optional flag.** Same charge is mandatory.
9. **Raising `assembly:max_fragments` to let the union through.** Separate experiment. n is already
   500 and open search pays n⁴.
10. **`ray_fraction` as a quoted quantity**, and **`truthset.py` RECALL as a kill gate.** §4.1.
11. **"Peptides flat while PSMs/peptide → 1.69" as the acceptance criterion.** §4.1: circular under
    `report_psms: 1`, and identical for merge and for deletion.
12. **The positional-counter scheme for recovering parent trace ids.**
    `ElutionPeakDetection::detectPeaks` (`ElutionPeakDetection.cpp:320–341`) is `#pragma omp parallel
    for` writing `single_mtraces` in thread-completion order under a critical section. Sub-traces
    from different parents interleave run-dependently; a counter would mis-assign parents *silently*.
    If provenance is ever needed, parse `MassTrace::getLabel()` — `MassTraceDetection.cpp:672` sets
    `"T" + trace_number` and `ElutionPeakDetection.cpp:534` appends `"." + String(min_idx + 1)`.
    Caveat: `trace_number` is local per `run()` call, so labels collide across MS2 bands.
13. **The baseline figures 8,524 and 6,940.** Zero hits anywhere in the repository, at HEAD or in any
    revision. The only committed record of that experiment is **8,925 → 6,628**
    (`docs/AlphaDIA-and-DIA-BERT.md:189`), consistent with src:531's "lost 26% of peptides"
    (−25.7%). Do not cite the others.
14. **`Δ(RT) = 0.000 s` as evidence of shared provenance.** Frame-grid quantisation (§1), proven from
    OpenMS source and confirmed by measurement.
15. **`union/larger = 1.355` as evidence of complementarity.** Poles inverted; 1.355 is near the
    *similar* end of a 1.0–1.95 scale (§1).
16. **"41.1% of peaks unique to one member."** Mislabelled statistic, sitting near its floor (§1).
17. **The 0.17% "replicate spread."** Not replicate data; the varied parameter was structurally
    inert (§4.4).
18. **Peaks/spectrum figures 228 and 23.** Measured: median 500, mean 435.8, 75.7% at the cap.

---

## 6. Repository corrections required before this ships

1. `docs/benchmarks/split-spectra.md` — correct the `union/larger` poles, remove the "41.1% unique"
   row, and annotate `Δ(RT) = 0.000 s` as frame-grid quantisation with the mechanism cited.
2. `harness/collate.py:124` ("we ~228") and `scripts/truthset.py:9` ("we emit ~23") — both wrong;
   correct to median 500 / mean 436, 75.7% at cap.
3. `src:1745`/`:1755` — the `consolidate` IM gate reads `MSSpectrum::getDriftTime()` (never set,
   default −1) instead of `getPrecursors()[0].getDriftTime()`. Fix, and re-run any consolidate sweep
   whose conclusion depended on IM.
4. `scripts/truthset.py:28` says "26 variable-width diaPASEF windows"; the list at `:29–33` holds
   **24**. Spectra in a missing window are counted `outside_window_scheme` and dropped from both OWN
   and CO-ISO (`:113`). Re-derive the scheme from the acquisition method.
5. `docs/benchmarks/joint-2026-07-22.md:19` and `README.md` — the 5.6× PSMs/peptide excess is a
   mixed-acceptance-set artefact; the single-set value is **4.67×**.
6. `README.md` OWN 5.79% / CO-ISO +7.50% — both lines date from the initial commit `a9a9b8a` and
   were never updated by `579ad70`, the commit that set the *current* defaults, and name no sample.
   **There is no valid pre-merge purity baseline for the shipped configuration.** Re-measure on a
   named sample before §4.2 runs; judging merge against a pre-defaults, sample-unnamed reference
   would be the mismatched-baseline error a fourth time.

---

## 7. What we do not know

Stated explicitly so these are not quietly assumed:

- **Whether count collapse helps anything.** Six redundancy-reducing levers have been falsified.
  The `consolidate:delta_rt` sweep gives a floor, not an answer, because it applies count collapse
  *and* deletion of the non-shared peaks simultaneously. Arm C in §4.2 is the same confound; only
  B−C isolates union.
- **Whether the 26% loss from frame aggregation was mechanistic or a metric artefact.**
  `REVIEW-2026-07-22.md:451` — "Under peptide_q, redundancy is free score" — predicts that every
  redundancy-reducing lever loses regardless of spectrum quality. That prediction is unfalsified and
  applies to this change too. It is the reason §4 does not gate on closed peptide count.
- **What fraction of emitted peaks are real.** ~91% of peaks are unannotated under b/y z≤2 / 20 ppm
  and ~96% under z=1 / 10 ppm; Sage's own `matched_peaks` median of 12 out of 500 agrees with the
  stricter figure. Some unannotated peaks are genuine (neutral losses, a-ions, isotopes, z>2), so
  every purity figure here is a **lower bound**. It is the *same* lower bound in every arm and every
  null, which is what the comparisons rest on — but the absolute level is unknown.
- **Full-group union gain.** `split_test.py:103–104` compares only `min(len(v), 3)` members
  pairwise, so 1.355 is a *pairwise* figure. With mean multiplicity 2.85 and 370 peptide-cycle pairs
  at 8+ members, the full-group figure is larger and unmeasured. This cuts in the design's favour and
  is exactly why it must not be quoted before it is measured.
- **The false-capture population's ΔIM.** Never measured. `merge:delta_im = 0.01` spends 52% of its
  budget on the target population's control (0.0052) with no measurement of what else it admits.
- **Whether the across-cycle result generalises.** One sample (dataset B), one arm (shipped defaults,
  `ms1_split_valleys` on), no replicate. The far control is n = 210. Conclusions §2.1, §2.2 and §2.4
  do not depend on that arm; §2.3 does.
