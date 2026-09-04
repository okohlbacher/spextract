# Merging split spectra: works, but the ceiling is zero (dataset D-family, small file)

## The question

We emit 9.51 PSMs/peptide against the reference implementation's 1.69. Adjacent-cycle spectra of the same peptide
have median cosine **0.687** (within-cycle 0.845, unrelated 0.007), so they are similar enough
to merge in principle. Does merging help?

## Coordinate-gated merging loses, monotonically

Gate = (RT window, 20 ppm m/z, 0.02 1/K0, same charge). Intensities summed.

| window | absorbed | peptides | change |
|---|---|---|---|
| none | — | 1,991 | — |
| 1.4 s | 49,187 | 1,918 | −3.71% |
| 2.6 s | 59,625 | 1,841 | −7.58% |
| 5.0 s | 75,008 | 1,700 | −14.66% |

This was first measured with the ion-mobility gate **dead** (`consolidate`/`merge` read
`MSSpectrum::getDriftTime()`, never set, returning −1). Repairing it to read the **precursor**
drift time changed nothing (−4.4% → −3.71%), so the broken gate was not the cause.

## Two hypotheses tested and refuted

**"Summing distorts intensities."** Wrong, and the fix made it worse. Taking the max instead of
the sum costs a further **11 points** (−15.01% vs −3.92% at 1.4 s). Summation is not a problem
for identification.

**"The merge fuses unrelated precursors into chimeras."** Wrong. Measured over adjacent
identified pairs:

| | median Δm/z | median Δ1/K0 |
|---|---|---|
| SAME peptide | **9.0 ppm** | **0.0074** |
| DIFFERENT peptide | **205,065 ppm** | **0.0906** |

Four orders of magnitude apart. At the shipped gate the **chimera rate is 1.8%**, and tightening
to 5 ppm / 0.005 only reaches 1.3%. Coordinates already discriminate almost perfectly.

## The actual mechanism

Only **12,380 of 161,798** emitted spectra identify at all. A coordinate gate therefore merges
good spectra with *unidentifiable* ones, whose peaks dilute the matched fraction without
contributing matches. PSMs per **spectrum** confirms it: 0.432 → 0.398 (−8%) — merged spectra
identify *less* often, which no counting argument explains.

## A content gate fixes the harm — and gains nothing

Requiring a minimum cosine before merging:

| gate | absorbed | peptides | change |
|---|---|---|---|
| none | — | 1,993 | — |
| cos ≥ 0.5 | 31,739 | 1,964 | −1.41% |
| cos ≥ 0.7 | 14,985 | 1,976 | −0.80% |
| **cos ≥ 0.85** | 4,500 | **1,992** | **+0.00%** |
| cos ≥ 0.95 | 732 | 1,983 | −0.45% |

At cos ≥ 0.85 merging is **exactly free** where coordinate gating cost 7%. But the trade curve
is unfavourable: safe merging touches 2.8% of spectra, so redundancy barely moves; merging
enough to matter (20% of spectra) costs 1.41%.

## Conclusion

Merging is **neutral at best**. The redundancy is not recoverable by combining spectra that
already exist. If it is reducible at all, it must be upstream — by not emitting spectra that
were never going to identify.

A joint probabilistic gate (Δm/z, Δ1/K0, cosine → logistic regression) was considered and
deprioritised: coordinates already give a 1.8% chimera rate, so better discrimination cannot
change a ceiling that is zero. The same machinery has real headroom applied to charge and
monoisotope assignment instead (74.6% vs a 69.6% constant).

## Parameters added

`merge:rt_window`, `merge:min_cosine`, `merge:sum_intensity`, `merge:mass_ppm`,
`merge:delta_im`, `merge:any_charge`. All default OFF.
