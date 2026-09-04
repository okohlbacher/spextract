# The within-cycle repeats are SPLIT SPECTRA, not duplicates (dataset B)

**Hypothesis (user, 2026-07-22):** either our spectra are misannotated (precursor mass, charge),
or *"the peaks of one spectrum split across two or more spectra — the opposite of chimaeras."*

The second is correct. Every prediction of the split model fires; every prediction of the
distinct-features model fails.

| observation | split predicts | distinct predicts | **measured** |
|---|---|---|---|
| **Δ(RT)** | ~0 | separated | **0.000 s** (control 1.4–4.2 s) |
| **Δ(IM)** | ~0 | separated | **0.0045** (control **0.0052**) |
| Jaccard peak overlap | low | incidental | **0.374** |
| peaks unique to one member | high | moderate | **41.1%** |
| **union / larger member** | **≫1** | ~1 | **1.355** |
| identical charge | — | — | 93.2% |
| Δ(precursor m/z) | small | large | 7.5 ppm median |

## The decisive numbers

**Δ(RT) = 0.000 s at every quartile.** Not small — zero. Pair members occupy the *same frame*.
That is one feature emitted twice, not two features.

**Δ(IM) = 0.0045 against a control of 0.0052.** Within-cycle pairs are *no more separated in
mobility* than the same peptide from different cycles. There is no mobility separation that would
justify treating them as distinct.

**Union/larger = 1.355** with **41% of peaks unique to one side.** Merging a pair gives a spectrum
36% richer than its largest member. The parts are **complementary**, not redundant.

## What this corrects

The cosine measurement earlier the same day gave median **0.845** and was read as *"partially
distinct spectra … merging destroys signal."* **That interpretation was wrong.** Splitness
predicts exactly that intermediate cosine — two halves of one spectrum share some peaks and each
carry unique ones. Cosine alone cannot separate the hypotheses; (IM, RT) geometry and peak
complementarity can, and they point the other way.

Consequences:

* **"Within-cycle duplication" was the wrong name** for this population, used for a full day.
* **`consolidate:delta_rt` is the wrong operation.** It keeps the richest member and *discards*
  the rest — throwing away the 41% of peaks unique to the others. The sweep running at the time
  of writing should lose peptides, and its pre-registered kill criterion should fire.
* **Merging, not pruning, is indicated.** Union recovers 36% more peaks per pair.
* **Hypothesis (a), misannotation, is largely excluded for this population**: 93.2% share charge
  and precursor m/z agrees to 7.5 ppm. It remains live for the *other* population — the ~44% of
  within-cycle pairs separated by ~1000 ppm ≈ one isotope spacing.

## Likely mechanism

Δ(RT)=0, 93% charge agreement and 7.5 ppm m/z agreement place the split **inside one frame, at
trace assembly**: one precursor's fragments partitioned across several emitted spectra instead of
assembled into one. Prime suspect is `trace:ms1_split_valleys`, made a default on 2026-07-21 on a
+46% closed-search gain, or the per-window fragment-to-precursor assignment.

Both are testable by re-running the split test on the `ms1_split_valleys 0` arm.
