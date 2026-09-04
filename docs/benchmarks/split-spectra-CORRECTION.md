# Correction to the split-spectra analysis (2026-07-22, same day)

An adversarial review found errors in [split-spectra.md](split-spectra.md). Three are confirmed
by independent re-measurement and are mine.

## 1. Δ(RT) = 0.000 s was TAUTOLOGICAL, not evidence

Reported as the decisive finding: *"Δ(RT) = 0.000 s at every quartile. Not small — zero. Pair
members occupy the same frame."*

**Sage RT is quantised to the frame grid.** Measured: 1,007 distinct RT values across 400,000
PSMs, modal spacing **1.385 s** — exactly the cycle period. Groups were built by binning
`rt/1.385`, so **every member of a group shares one RT by construction.** Reporting Δ(RT)=0
within those groups is circular. It carries no information about splitness.

The other split evidence (Jaccard 0.374, union/larger 1.355, 93.2% same charge, 7.5 ppm m/z)
does not depend on this, but the single most persuasive number is withdrawn.

## 2. `union/larger = 1.355` was read with its poles inverted

Reported as *"the parts are complementary, not redundant."*

Identical sets give 1.0; disjoint equal-sized sets give 2.0. **1.355 is much closer to identical
than to disjoint.** It means substantially overlapping with a modest unique fraction — not
complementary. The 36%-richer-union figure stands as arithmetic; the word "complementary" does
not.

## 3. Peaks per spectrum: the committed figure was from the wrong file

Reported median 228 (ours) vs 500 (the reference implementation). Re-measured on the actual `split_count` dataset B file:
**median 437, and 46.7% of spectra already at the 500-peak cap.** The 228 came from the older
dataset A `basepin2` output.

Consequence for the merge proposal: with ~47% of spectra at the cap, **merging cannot add peaks
to them.** Much of the union gain is unrealisable without raising `max_fragments`.

## 4. CODE BUG: consolidate's ion-mobility gate has never fired

`src/spextract.cpp:1137` stamps IM on the **Precursor** (`prec.setDriftTime(pc.im)`).
Consolidation at **:1745** and **:1755** reads it from the **MSSpectrum**
(`all_out[i].getDriftTime()`), which is never set and returns −1. The guard is

```cpp
if (im_i > 0 && im_j > 0 && fabs(im_j - im_i) > cons_im) continue;
```

so `im_i > 0` is false and the mobility check is skipped entirely.

**This confounds the consolidation KILL result.** The sweep reported at −12.7% peptides
(`drt=0.7`) merged on (RT, m/z, charge) with **no mobility constraint**. That is an
over-permissive merge, so the peptide loss is not clean evidence that pruning is wrong —
it is evidence that pruning *with a broken gate* is wrong. The sweep must be repeated after
the fix before the conclusion stands.

## 5. Partially disputed: the "0.17% replicate spread"

The review argues the four `ambiguity_margin` runs are not replicate data. The parameter being
inert is precisely what makes the four configs identical, so they *are* replicates — but the
runs were not all on one node: three on node-2 (8,411 / 8,408 / 8,406) and one on data (8,397).

**Same-node spread is 0.06%; the 0.17% figure mixes in cross-node variation.** The tighter number
should be used as the noise floor, and it is still a real floor rather than none.

## What survives

The split hypothesis is **weakened but not dead**. Surviving evidence: Jaccard 0.374, 41% of peaks
present in only one member, 93.2% charge agreement, 7.5 ppm precursor agreement, and
`ms1_split_valleys` quadrupling the affected groups (3,398 → 14,379). Withdrawn: Δ(RT)=0 as
evidence, and "complementary" as a description.
