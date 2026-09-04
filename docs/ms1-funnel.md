# MS1 funnel: where precursors are actually lost (2026-07-21, dataset A)

Anchored on the DIA-NN dataset A reference (43,499 precursors at q<=1%). Our tool dumped its
MS1 traces and inferred precursors via `diag:dump_ms1_tsv` (2,521,685 traces →
1,206,174 hypotheses); each DIA-NN precursor was then attributed to the stage where it
dies. A precursor counts as "detected" if a trace sits at the mono **or the +1/+2
isotope** m/z — a trace on the wrong isotope still proves the ion was seen, which is
what separates "invisible" from "mislabelled".

| stage | count | share | failing code |
|---|---|---|---|
| A — no MS1 trace at all | 977 | **2.2%** | `detectTraces_` sensitivity |
| B — trace, but no hypothesis at the mono | 4,585 | 10.5% | `inferPrecursors_` monoisotope walk |
| B — hypothesis, but **wrong charge** | 6,922 | **15.9%** | `inferPrecursors_` charge assignment |
| C — correct hypothesis | 31,015 | 71.3% | — |

> # ⚠️ RETRACTED (2026-07-21, same day): the numbers above are not supported
>
> A decoy floor was added after the fact — the same DIA-NN precursors shifted +11.003 Da, so
> they do not exist — and re-run through the identical pipeline:
>
> | stage | observed | **decoy floor** |
> |---|---|---|
> | no MS1 trace at all | 2.2% | **8.1%** |
> | trace but no hypothesis | 10.5% | **24.0%** |
> | hypothesis, wrong charge | 15.9% | **19.4%** |
> | correct hypothesis | 71.3% | **48.4%** |
>
> **A precursor that does not exist finds an MS1 trace 91.9% of the time and a
> charge-matching hypothesis 48.4% of the time.** With 2.5M traces and 1.2M hypotheses, a
> 20 ppm x 0.3 min x 0.05 1/K0 box — accepting a hit at the mono *or* +1 *or* +2 isotope —
> is so permissive that matching is nearly free. The assay was measuring **trace density,
> not correctness**.
>
> Chance-corrected: detection ceiling **72.3%**, correct-hypothesis **44.4%**.
> The wrong-charge stage (15.9%) sits *below* its own chance floor (19.4%) and therefore
> carries **no information at all**.
>
> **"97.8% MS1 detection ceiling" and "26.5% lost in inference, not detection" are
> withdrawn.** They were quoted in the README, a commit message, and repeatedly in
> discussion before the control existed. The qualitative direction (detection cheaper than
> inference) may still hold, but nothing here establishes it.
>
> A tolerance sweep is running to find an operating point where the floor is low enough for
> the observed-minus-floor gap to dominate. Until then, no stage attribution should be quoted.
>
> **Process note:** the spectrum-level recall measurement in the same session *did* get a
> decoy floor (35.6%, which cut 74.3% -> 60.1%). The funnel did not, and the caveat section
> below even said so — while the headline number was quoted unqualified anyway. Writing the
> caveat is not the same as heeding it.

**MS1 detection ceiling: 97.8%.** We see essentially every precursor DIA-NN finds.
**26.5% is destroyed downstream of detection, inside `inferPrecursors_`.**

## What this rules out

Extraction-sensitivity work targets a **2.2%** problem:
- way-#5 (replace `findBestPeak_` with joint m/z–IM–intensity scoring) — deprioritised
- AlphaDIA-style threshold-free MS1 extraction — not the MS1 bottleneck
  (it may still matter for MS2 fragments, which this funnel does not measure)

Anything touching noise thresholds, `chrom_peak_snr`, `min_length`, or `min_sample_rate`
at MS1 is bounded above by 2.2%. Stop tuning them.

## What this points at

Charge assignment (15.9%) and monoisotope identification (10.5%). Corroborated
independently by the truth-set charge measurement: **47.8% charge agreement vs
the reference implementation's 82.6%**, dominant confusion `4→2` (15,734 cases). At m/z 500 a z4→z2 error
fabricates ΔM = (2−4)(500 − 1.0073) ≈ **−998 Da**, which is fatal for open search.

## Mechanistically consistent with the way-2 result

way-2 **both** (`ms1_split_valleys` + `ms2_split_valleys` = 7.0) gave **+54.3% peptides**
at peptide-level FDR (5,131 → 7,917 on dataset B). It is the arm that splits **MS1** traces.

Merged MS1 traces → merged precursors → wrong monoisotope and wrong charge — exactly the
26.5% this funnel localises. Two independent measurements identify the same defect, and
the fix moves the number the funnel says is broken. See [[dataset A-results]].

## Caveats

- Matching tolerance 20 ppm / 0.30 min / 0.05 1/K0. No decoy floor was computed for the
  funnel itself; unlike the spectrum-level recall measurement (which needed one, and where
  the floor was 35.6%), stage attribution is a *conditional* breakdown among precursors
  DIA-NN already asserts exist, so chance matching inflates C at the expense of A/B rather
  than manufacturing precursors. The 97.8% ceiling is therefore an upper bound.
- DIA-NN's space is tryptic / z2–4 / 2 mods, so this says nothing about non-canonical
  precursors — the population the tool actually exists to find.
