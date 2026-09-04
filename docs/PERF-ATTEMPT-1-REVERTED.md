# Performance attempt 1 — REVERTED (2026-07-23)

Two independent audits each found a real defect and each proposed an **output-identical** fix.
Both claims were wrong, and the A/B caught it.

## The measurement

Small file (PXD027359, 361 MB), same binary tree, same config, OLD vs NEW:

| | wall | peptides | mzML bytes |
|---|---|---|---|
| OLD | 95 s | **1,992** | 1,041,375,208 |
| NEW | 94 s | **1,849** | 1,026,093,099 |

**−143 peptides (−7.2%) for −1 s (−1%).** The gate was identical peptide counts. It failed, and
there was no speedup to trade against it either.

## What was attempted

**Fix 1 — cap the nested OpenMP team (`EpdThreadCap`).** The underlying defect is real and
verified: `ElutionPeakDetection.cpp:328` is a bare `#pragma omp parallel for` with **no
`num_threads`**; `spextract.cpp:1545` enables 2 active levels; `epd.detectPeaks` runs
inside the window loop. So 24 windows × 120 threads ≈ **2,880 threads on 120 cores**, all
funnelling through the *program-global* critical at `ElutionPeakDetection.cpp:460/545`. That
explains the measured 43.5% system time.

*Why the output-identity claim was wrong:* the critical guards a **shared output vector**.
Changing team size changes the order in which traces are appended, and valley splitting is
order-dependent — so the trace *set* is not thread-count invariant, only the trace *population*
before splitting is.

**Fix 2 — parallelise the per-precursor scoring loop.** Also a real defect: `/proc` sampling
showed the loop **ends with 1,993 s at 1.000 user cores, 21 of 22 threads in `futex_wait`** —
the last window's scoring is serial, ≥45% of wall at 1/128 utilisation.

*Why the output-identity claim was wrong:* the mzML shrank by 15 MB, so **spectra are being
lost, not merely reordered**. Independent iterations were assumed; something in the chunked path
drops emissions.

## What stands

* **Both defects are real and verified in the source.** The 43.5% kernel time and the 1,993 s
  single-core tail are measured, not inferred. The *fixes* were wrong, not the diagnosis.
* The 9.7× decomposes as **3.67× CPU × 2.64× parallelism** against the reference implementation (457 s, 38.4 cores,
  6.8% system on the same node and `.d`). SpeXtract's *user* CPU alone would finish in ~5 min at
  120 cores — the algorithm is not 10× slow, it runs on a ninth of the machine.
* Window `[327.52, 468.80)` holds 489,986 of 2,145,053 precursors (**23.7%**), capping the loop
  at **4.2×** by Amdahl regardless of thread count.

## Method note

The small file made this cheap: **95 s per arm instead of 74 min**, so a failed attempt cost
minutes. But its phase split is *inverted* versus dataset D — MS1 51% / loop 33% here, versus MS1 7.9%
/ loop 91.3% on dataset D — because loop cost scales with co-eluting precursors and a 6-min gradient
has far fewer. **It is a valid output-identity gate and an MS1-phase benchmark; it is not a
window-loop benchmark.** Window-loop timings must still be confirmed on dataset D.
