# Backlog

Ideas not yet tested, with the reasoning that makes them worth testing and the result that
would kill them. Falsified ideas move to the relevant docs/ file rather than being deleted —
knowing what does not work is most of what this project has established.

---

## 1. Unmix fragment XICs onto the co-eluting MS1 precursor basis (NNLS)

**The idea.** Today `scoreCandidates_` asks a *pairwise* question, one precursor at a time:
does fragment `f`'s elution profile correlate with precursor `p`'s? Every co-eluting precursor
asks independently, so under share-all a fragment is emitted into all of them at full intensity
(measured mean fan-out 6.45). Nothing ever asks the *joint* question.

The joint question is a regression. Let `P` be the matrix whose columns are the MS1 XICs of the
precursors co-eluting in this window, and `f` a fragment's XIC. Solve

    min || P w - f ||^2   subject to  w >= 0

The weights `w` are the fragment's intensity split across precursors, directly. Assign the
fragment to `argmax w`, or emit it into each precursor at `w_i` — an apportionment *derived from
the data* rather than a heuristic.

**Why this is better-posed than the decompositions already falsified.** Earlier analysis rejected
global SVD/NMF on this data because the rank is unknown and the matrix is not globally low-rank —
and rank selection *is* the answer you want, so assuming it is circular. That objection does not
apply here: **the number of components is not estimated, it is observed.** MS1 tells us which
precursors are present and gives their elution profiles. This turns an ill-posed blind
decomposition into a well-posed non-negative regression with a known design matrix. That is the
whole difference, and it is why this deserves a test where blind NMF did not.

**Relation to what exists.** `assembly:apportion` already splits intensity by `corr^p / sum(corr^p)`.
That is a crude, coordinate-wise stand-in for the same intent — it never solves a joint fit, so two
precursors with similar profiles both get high weight instead of competing. NNLS is the principled
version and subsumes it. `assembly:apportion` is currently **untested**, so measure it first: if the
heuristic already captures most of the gain, NNLS is not worth the complexity.

**Why it should help, in terms of measured defects.**
- Redundancy: we lose 14-16% of peptides at peptide-level FDR against the reference implementation's 5-7%, because we
  emit ~4.6 PSMs/peptide. Fragments landing in one precursor instead of six attacks that directly.
- Purity: only 5.79% (decoy-corrected, equal depth) of an emitted spectrum's intensity belongs to
  its own precursor; co-isolation contributes +7.50%.

**How to test.** The DIA-NN truth set already provides ground-truth fragment ownership: the
predicted library gives 405,336 true (precursor, fragment) pairs, and negatives are fragments of
other co-eluting precursors in the same window. So this is directly measurable *without* a search
engine — compute assignment precision/recall against that set for pairwise-Pearson vs `apportion`
vs NNLS, before running Sage at all.

**What would falsify it.**
- Assignment precision against the DIA-NN pairs does not beat the pairwise baseline. Then the
  limiting factor is the MS1 profiles themselves, not how they are combined.
- Precursor XICs in a window are too collinear to separate. **Measure this first** — it is cheap
  and decisive. Compute the condition number / pairwise correlation of `P` per window. Peptides
  co-eluting within a couple of seconds have near-identical profiles, and no solver recovers a
  split the data does not contain. If typical `P` is severely ill-conditioned, NNLS produces
  confident, arbitrary, wrong splits — worse than the current honest ambiguity.
- Peptides at `peptide_q` do not improve even if assignment precision does.

**Note on the ill-conditioning risk.** This is the real hazard, not compute cost. A regression that
silently picks one of many equally good answers is more dangerous than a heuristic that spreads the
intensity, because downstream it looks decisive. Regularisation (ridge, or an IM-similarity prior)
or an explicit ambiguity flag would be needed — and an ambiguity flag is arguably the more honest
output for a blind-search tool.

**Cost.** Windows hold ~50 co-eluting DIA-NN precursors and millions of fragment traces, so this is
a small dense NNLS (~50 unknowns, grid-length equations) per fragment. Warm-startable and trivially
parallel across fragments, but it sits in the inner loop, so runtime must be measured too.

**Prerequisite ordering.** Measure `assembly:apportion` -> measure conditioning of `P` -> only then
implement NNLS.

---

## N. Decompose the whole cycle jointly, using our own redundancy as the data

**Origin.** User proposal, 2026-07-22: *"use the large number of spectra to cluster the peaks in
the cycle (e.g. SVD/NNLS) to decompose everything"* — framed as **data cleanup**.

**The hook that makes this different from the blind decomposition already rejected.** Today's
cycle analysis ([benchmarks/cycle-redundancy.md](benchmarks/cycle-redundancy.md)) measured
within-cycle multiplicity of **2.85** — we emit ~3 pseudo-spectra per (peptide, cycle), where
the reference implementation emits 1.25. Those have been treated purely as a defect. But they are also **~3
independent noisy views of the same underlying cycle content**, and that is exactly the regime
where joint decomposition has something to work with. The redundancy stops being only a bug and
becomes the input.

Concretely: within one isolation window × one cycle, build `X` = (emitted spectra × fragment
m/z bins). We have many rows *because* of the over-emission. Factor `X ≈ W H` with `H` ≥ 0 the
component spectra and `W` ≥ 0 the per-emission weights. Components that recur across the
redundant emissions are real; components appearing in one emission only are extraction noise.
Output one consolidated spectrum per component — which is precisely the cleanup asked for, and
would attack the ~56% of within-cycle duplicates that `consolidate:delta_rt` addresses only by
nearest-neighbour merging.

**What the 2026-07-22 review established against blind decomposition, which still applies:**

* **Kruskal uniqueness needs component profiles to differ in at least one mode.** Co-eluting
  precursors inside one isolation window are exactly where profiles are collinear — PARAFAC is
  least identifiable precisely where we need it most.
* **Ion mobility is a partition, not a profile.** TIMS precedes fragmentation, so the mobility
  mode is a *label*, not a continuous third mode. The trilinearity argument for PARAFAC — the
  main reason it looked attractive — does not hold.
* **The quadrupole mask makes missingness structured**, not random: m/z visibility depends on IM
  and cycle position.
* Overall probability of success judged **<5%**.

**Why this variant partly escapes those.** It is not blind and not trilinear. It is a 2-D
non-negative factorisation over (emission × fragment) *within* one cycle-window cell, where the
row count is empirically large and the rank question is bounded by something observable — DIA-NN
reports ~50 co-isolated precursors per spectrum, and MS1 supplies candidate precursors for that
cell. It is closer to codex's recommendation of **anchored / coupled factorisation using
MS1-derived profiles as partial known factors** than to blind PARAFAC. It also overlaps
substantially with backlog item 1 (NNLS onto the MS1 basis) — **item 1 should be tried first**,
since it needs no rank estimate at all.

**Cheap falsification, before any implementation.** Take the ~31,621 duplicated (peptide, cycle)
groups already identified. For the ~56% whose members agree to ~15 ppm, ask offline: does a rank-1
non-negative fit reconstruct the members within noise? If yes, they are redundant views and joint
decomposition is justified. If the members need rank ≥ 2 to reconstruct, they carry distinct
information and merging them — by NMF *or* by `consolidate` — destroys signal. **That test is a
few hours on existing output and needs no new extraction run.**

**What would kill it.** Rank-1 reconstruction fails on the tight-m/z duplicates; or item 1 (NNLS)
already captures the gain, making the extra machinery unjustified; or the peptide count falls in
the `consolidate:delta_rt` sweep now running, which would show that within-cycle duplicates carry
information after all.

---

## N+1. Learned monoisotope/charge predictor from MS1 + fragment features

**Origin.** User proposal, 2026-07-22: *"correlate with MS1 precursor masses, perhaps in an ML
based predictor including these features."*

**Why this is now well-posed, when it was not this morning.** Three things changed today:

1. **The target is quantified.** Charge agreement is **74.6%** against a majority-class baseline
   of **69.6%** — inference beats "always answer z=2" by 5.0 points. That is the headroom, and
   it is large.
2. **Ion mobility is back in play.** The IM lever was retired on a statistics error: a
   *share-a-label* rate was compared against `max p_z` (69.6%) instead of `sum p_z^2` (55.8%).
   Corrected, (m/z, IM) carries **+19.2 points** of charge information over chance, not +5.4.
   See [charge-inference.md](charge-inference.md). A *monotonic* IM prior still fails — z=4 lies
   between z=2 and z=3 — but a **learned 2-D prior does not**, and that is precisely what a model
   can represent and a threshold cannot.
3. **The failure mode is localised.** ~44% of within-cycle duplicates are separated by ~1000 ppm
   ≈ one isotope spacing: the same peptide assigned to different monoisotopes. That is the thing
   to predict.

**Labels exist.** DIA-NN references give true charge and monoisotope for ~43,499 precursors on
dataset A, plus dataset B and dataset D — ~130k labelled precursors. Small for a transformer (DIA-BERT needed
276M) but ample for the model class the review actually recommends.

**Model class, per the 2026-07-22 review.** AlphaDIA — state of the art in this exact data type —
scores with **46 hand-crafted features into a 10,810-parameter MLP**, *not* a transformer. Verdict
there: transformer ~8% probability of beating the hand-crafted path; the small model is the
recommendation. Start with gradient-boosted trees or an MLP of that size.

**Feature set (all already computed or cheaply derivable):**

| group | features |
|---|---|
| MS1 isotope | intensity ratios M/M+1/M+2, observed-vs-averagine cosine, spacing residual in ppm, number of partners found |
| MS1–MS2 coupling | correlation of fragment XICs to the precursor XIC, fraction of fragments correlating, apex RT offset |
| mobility | 1/K0, deviation from the (m/z, IM) charge-2 line, and the tile's `ScanNumBegin/End` band — a hard transmission constraint the code currently discards |
| position | m/z, RT, precursor intensity, local co-isolation density |
| redundancy | within-cycle multiplicity for this precursor (measured 2.85 mean) |

**Target.** Correct charge, and correct monoisotope offset (−1 / 0 / +1 / +2), as two heads or
one joint class.

**Falsification, pre-registered.** The model must beat **both** baselines:
* the 69.6% majority-class constant, and
* the 74.6% shipped `count` heuristic,
by more than the 0.17% replicate spread, on a **held-out sample** — train on dataset A+dataset D, test on
dataset B. Training and testing on one sample would measure memorisation, and this project's
characteristic failure is exactly a metric that cannot tell the two apart.

**What would kill it.** Held-out charge accuracy does not beat 74.6%; or it does, but the peptide
count does not move — which would mean charge was never the binding constraint and the real limit
is elsewhere.

**Do first.** The `isotope_errors [0,0]` experiment now running measures how much the *search
engine* is already repairing our monoisotope errors. If Sage's correction fully absorbs them, a
predictor buys little for closed search and its value is confined to open search — which is the
stated purpose, but changes the cost/benefit.
