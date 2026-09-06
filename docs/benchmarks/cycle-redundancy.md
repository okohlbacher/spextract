# Where the redundancy comes from (dataset B, Sage, common FDR)

SpeXtractor emits ~9.5 PSMs/peptide against the reference implementation's ~1.7. That ratio alone cannot separate two
very different situations, so it factorises exactly:

```
PSMs/peptide  =  (PSMs per (peptide, cycle))  ×  (cycles per peptide)
                  within-cycle multiplicity       elution spread
```

Only the first term is a defect: **one diaPASEF cycle contains one elution point for one
precursor.** Cycle period 1.385 s, read from `analysis.tdf` (1,342 cycles, 16,105 MS2 frames /
12 window groups), not assumed.

| | PSMs/pep | **within-cycle** | across-cycle | **% within-cycle dups** |
|---|---|---|---|---|
| **SpeXtractor** | 7.57 | **2.85** | 2.66 | **64.9%** |
| **the reference implementation** | 1.62 | **1.25** | 1.29 | **20.3%** |

**64.9% of all our PSMs are within-cycle repeats**, against the reference implementation's 20.3%.

## The histogram

PSMs per (peptide, cycle):

| | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| SpeXtractor | 12,318 | 5,221 | 2,915 | 1,889 | 1,206 | 838 | 567 | **370** |
| the reference implementation | 9,894 | 1,912 | 409 | 85 | 26 | 6 | 1 | 0 |

the reference implementation's tail dies at 4. Ours has **370 peptide-cycle pairs yielding 8+ PSMs** — one cycle,
one peptide, eight or more separate identifications.

## Why this is a different lever from everything falsified so far

Every redundancy knob tested to date — `assembly:competitive`, rank-pruning (`rp_max`),
`assembly:apportion`, frame aggregation — operates at the **fragment-to-precursor** level, and
every one lost peptides at `peptide_q`. This operates at the **precursor-within-cycle** level,
where duplication has *no physical justification at all*.

The across-cycle term (2.66 vs 1.29) is arguably legitimate: a peptide genuinely elutes over
several cycles, and the reference implementation may be under-reporting by collapsing them. **The within-cycle term
is not.**

## Bounded prediction

Collapsing within-cycle duplicates to one PSM should cut PSM count ~2.85× with **no expected
peptide loss** — PSMs/peptide 7.57 → ~2.66, against the reference implementation's 1.62.

Falsifier: if peptide count drops materially, the duplicates were carrying information and this
is just another count knob, like the four already falsified.
