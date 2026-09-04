# Charge inference: what the 15.9% is, and what will not fix it

The [[ms1-funnel]] localises **15.9%** of precursor loss to wrong charge assignment and
**10.5%** to failed monoisotope identification — against only 2.2% undetected. Truth-set
charge agreement is **47.8% vs the reference implementation's 82.6%**, dominant confusion `4→2` (15,734 cases).

## Why 4→2 happens

Isotope spacing is 1/z. At z=4 that is 0.25 Da; at z=2, 0.50 Da. A spurious peak halfway
between two real isotopes is indistinguishable from a genuine z=4 partner **on spacing
alone** — and the higher-charge hypothesis always has more chances to find one, so the
error is systematically biased toward overcalling charge.

Cost in open search: ΔM = (z′−z)(m/z − 1.00728). Measured at m/z 600:

| ours | true | phantom ΔM |
|---|---|---|
| z=4 | z=2 | **−1198.0 Da** |
| z=3 | z=2 | −599.0 Da |
| z=1 | z=2 | +599.0 Da |

## RETRACTED FALSIFICATION: ion mobility was retired on a wrong null

**This section previously read "FALSIFIED: ion mobility cannot break the tie". That conclusion
rested on a statistical error and is withdrawn.**

The nearest-neighbour test asks whether a precursor's nearest neighbour in (m/z, IM) **shares**
its charge. That is a coincidence rate between two draws, so the null is `sum(p_z^2)`. The
script compared it against `max(p_z)` instead — the null for *predicting* one label:

| | value |
|---|---|
| null used (`max p_z`) | 69.6% |
| **correct null (`sum p_z^2`)** | **55.8%** |
| observed | 75.0% |
| margin as published | +5.4 points |
| **margin corrected** | **+19.2 points** |

So the (m/z, IM) plane carries **substantially more** charge information than reported, and
"do not build an IM-based charge prior" was wrong.

What still stands, because it is a direct observation rather than a null comparison:

Hypothesis was that diaPASEF's mobility dimension constrains charge — the acquisition
tiles themselves track the charge-2 line, so charge ought to be separable in (m/z, 1/K0).
Tested against DIA-NN dataset A ground truth (43,499 precursors with m/z, IM **and** true charge):

| m/z 600–700 | z=2 | z=3 | z=4 |
|---|---|---|---|
| 1/K0 | 0.993 ± 0.031 | 0.937 ± 0.044 | **0.966 ± 0.035** |

**z=4 lies between z=2 and z=3 — mobility is not monotonic in charge.** At fixed m/z a
higher-charge peptide is proportionally heavier (mass = z × m/z), so its larger CCS
cancels the extra charge, non-monotonically.

Nearest-neighbour test (m/z ±2 Da, closest IM): **75.0% share charge vs a 69.6% chance
baseline — +5.4 points.** The plane carries almost no charge information beyond the base
rate. Do not build an IM-based charge prior.

This also retires the broader "use `ScanNumBegin`/`ScanNumEnd` to gate charge" idea. The
tile IM range is still worth using as a *transmission* filter (a precursor outside the
tile's mobility band could not have produced its fragments), but not as charge evidence.

## What is left

Isotope-envelope **intensity** matching: for a genuine z=4 precursor of mass ~2400 Da the
averagine model predicts specific M/M+1/M+2 abundance ratios. A spurious +0.25 Da peak
will not match them. Spacing is degenerate between z and 2z; **abundance is not**.

## CORRECTION (2026-07-21): envelope scoring is default-ON, not default-OFF

An earlier version of this file claimed `charge:scoring=envelope` was an untested lever
waiting to be enabled. It **was** the default until 2026-07-21; the code now ships
`registerStringOption_("charge:scoring", "<mode>", "count", ...)`. A benchmark arm passing
`-charge:scoring envelope`
returned byte-identical output to baseline (6509 spectrum_q / 5813 peptide_q), which is
what exposed it.

This matters for how every number here is read:

* The **15.9% wrong-charge rate and 47.8% charge agreement are what envelope scoring
  already produces.** They are not a "before" measurement.
* Averagine-shape cosine x isotope-XIC co-elution, with gapped envelopes allowed, is
  therefore **not sufficient** to break the z vs 2z degeneracy. The abundance argument
  above is correct in principle but the current implementation of it is not working.
* The informative experiment is the opposite direction: `-charge:scoring count` (the
  legacy partner-count heuristic, ties favouring LOW charge, break-on-first-miss) to
  establish whether envelope is better or worse than the alternative. Until that runs,
  we do not know that envelope scoring is helping at all.

The dominant confusion being `4->2` (overcalling charge) is consistent with a scorer that
rewards finding MORE envelope members: a z=4 hypothesis has twice as many candidate slots
as z=2 in the same m/z span, so gapped-envelope tolerance may actively favour it. Worth
checking whether the evidence weighting penalises hypothesis complexity enough.


---

# Moved here from the README (2026-09-04)

The README became a tool description; these measurements belong with the topic. Verbatim,
except that the comparison tool's name is written as `the reference tool`.

## Charge inference: the two defaults are coupled, not independent

Measured against a DIA-NN dataset B reference (`charge agreement` = our charge vs truth on matched
spectra; always answering z=2 would score **69.6%**):

| arm | charge agreement | `4->2` errors |
|---|---|---|
| ms1split + ambiguity_margin | **41.7%** | 34,554 |
| ms1split alone | 45.8% | 19,995 |
| baseline (`envelope`) | 49.7% | 10,199 |
| `count` | 71.8% | -- |
| **`ms1split` + `count`** (default) | **74.6%** | -- |

**MS1 splitting must not ship without `count`** — an *empirical* coupling, not an architectural
one. The two touch separate code paths (`trace:ms1_split_valleys` in `detectTraces_`,
`charge:scoring` in `inferPrecursors_`); nothing in the code enforces the pairing, so the
combination is a benchmark finding that a future edit could silently break. On its own MS1
splitting drops charge agreement below baseline and
more than triples the `4->2` confusion, which fabricates a ~-1198 Da phantom modification at
m/z 600 -- invisible to closed search (the engine enumerates charge) and fatal to open search.
Shipping `trace:ms1_split_valleys` without `charge:scoring count` would have been actively
harmful for this tool's actual purpose.

### `charge:ambiguity_margin` -- falsified, four tests

| margin | peptide_q | vs `split_count` 8,411 |
|---|---|---|
| 0.5 | 8,408 | -0.04% |
| 1.0 | 8,406 | -0.06% |
| 2.0 | 8,397 | -0.2% |

2.0 is the first value that can genuinely fire against an integer partner count (it admits any
hypothesis within 2 partners of the best) and it still changes nothing. Hedging charge does not
help; deciding it well does. Left at 0.
