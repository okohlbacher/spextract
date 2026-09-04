# MS1 feature FWHM vs retention time (dataset D)

Measured per MS1 trace from its own XIC — apex, then outward to half maximum with linear
interpolation. Traces that never fall to half max on both sides are reported as −1 and
**excluded**: a truncated peak is not a narrow one, and imputing it would bias the curve
downward exactly where peaks are broadest.

```
4,915,724 traces total
1,601,434 truncated  (excluded)
  317,505 too few XIC points (excluded)
2,996,785 used
```

## Result

**Median FWHM 3.61 s = 2.61 cycles.** Flat across the gradient: Theil–Sen slope
**0.0075 s/min**, ratio 8→28 min = **1.04×**. Slightly U-shaped (3.74 s at the start, 3.49 s
mid-gradient, 3.87 s at the end) but the whole excursion is ±5%.

| RT (min) | median FWHM | cycles |
|---|---|---|
| 0.0–9.8 | 3.74 s | 2.70 |
| 14.8–16.3 | 3.50 s | 2.52 |
| 26.6–31.0 | 3.87 s | 2.79 |

## What this settles

**1. The RT-dependence hypothesis is not supported.** Peak width does not grow with elution
time in this data, so the earlier "no breakpoint in the co-elution statistic" finding is *not*
explained by a smeared RT-dependent window. A single global `rt_window` is defensible.

**2. Our across-cycle spread is CORRECT, not excess.** We spread a peptide over **2.66 cycles**;
the chromatographic peak is **2.61 cycles** wide. Those are the same number. The across-cycle
term has been miscounted as redundancy all day — it is the one part of our behaviour that
matches the physics.

That splits the 9.51 PSMs/peptide cleanly:

| term | value | verdict |
|---|---|---|
| across-cycle spread | 2.66 | **correct** — equals FWHM |
| within-cycle multiplicity | 2.85 | **the defect** — one frame, ~3 spectra |

**3. It supplies the window the review said could not be derived.** One FWHM ≈ **3.6 s**, which
is **13× narrower** than the 45.7 s at which an earlier analysis still saw "benefit". Both can be
true: real structure inside one peak width, artifact beyond it.

**4. the reference implementation deliberately under-samples.** 1.29 cycles/peptide against a 2.61-cycle peak — it
takes roughly the apex only, and gets *more* unique peptides. So sampling the full peak is not
obviously better for identification: the apex has the best S/N, and one clean spectrum can beat
three noisy ones.

## Correction to an earlier claim

An earlier note put the elution window at ~7.5 cycles and concluded we might be *under*-covering
the peak. That came from DIA-NN's `RT.Start`/`RT.Stop`, which is a frame-bounded **integration
window**, not a fitted FWHM — and it is quantised to 1.385 s, so it cannot resolve a width change
below ~14%. **Real FWHM is 3.61 s. The "under-covering" claim is withdrawn.**
