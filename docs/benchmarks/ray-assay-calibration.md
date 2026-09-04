# Calibrating the delta-mass ray assay: it recovers ≤2.9% of charge errors

## What was done

`corrupt_charge.py` rewrote the reported charge of a known fraction of precursors, leaving m/z
untouched — exactly the error the tool makes in the field. The corrupted spectra were searched
through the symmetric −3000/+3000 Da open window and analysed by the ray decomposition.

| arm | rewrite | k | injected | rays detected |
|---|---|---|---|---|
| `4to2` | true z=4 reported as z=2 | −2 | 3,865 | **0** |
| `2to4` | true z=2 reported as z=4 | +2 | 99,984 | **0** |

## The bound, stated correctly

Naively 0 of 103,849 gives a 95% upper bound of 0.0029%. **That figure is wrong, and the error
is in the experiment design, not the arithmetic.**

Errors were injected into a **random** 10% of eligible spectra. But only **0.103%** of PSMs
survive 1% entrapment FDR in open search (1,512 of 1,468,765). So of 99,984 injected errors, only
**~103 landed in spectra that could ever have been identified** — the rest were injected into
spectra that produce no PSM regardless of their charge.

**Correct bound: 0 of ~103 detectable → 95% upper bound 2.9% recovery.**

Still decisive for the practical question: the assay recovers at most 1 charge error in 34.

## Consequence

**The ray assay cannot quantify charge errors.** The 26 slope-+2 events previously observed in
the `envelope` arm cannot be converted into a real count — at ≤2.9% recovery they would imply
≥900 true errors, but that inference rests on a bound derived from ~103 events and should not be
quoted as a number.

The earlier claim that "the control fired at p ≈ 2.4×10⁻⁶" remains arithmetically true about
excess over base rate, and remains uninformative about whether the instrument can *measure*
anything.

## Why recovery is so low

Not a defect in the geometry — the ray/band separation works, and was validated on synthetic
data at 300/300 recovery with zero noise leakage. The limit is upstream: **a spectrum with a
corrupted precursor mass rarely survives 1% entrapment FDR at all.** Open search on pseudo-DDA
spectra retains 0.1% of PSMs, so almost everything is filtered before the geometry can see it.

That is a property of open search on this data, not of the script.

## What would make this measurable

Inject into **spectra that are already identified** in an uncorrupted run, then measure recovery
on that population. That answers the question this experiment intended to answer — "can the assay
see a charge error in a spectrum that would otherwise identify?" — and gives a bound with useful
power instead of one resting on ~103 events.

Until then, the **DIA-NN-anchored charge measurement (74.6% agreement against a 69.6%
majority-class constant) remains the better instrument for charge quality**, despite its own
limitation of being blind to non-canonical species.
