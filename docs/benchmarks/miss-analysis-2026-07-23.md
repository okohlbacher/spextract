# Why we lose peptides: the labelled miss set (real-DDA matched pair)

1,590 peptides real DDA identifies and our pseudo-DDA does not. First **labelled** miss set in
the project — every miss has a known sequence, charge, and theoretical mass.

## First: "48.2% recall" was the wrong word (codex + vibe, confirmed)

Both external reviewers rejected the term, and they are right. DDA and DIA are **separate
injections sampling the peptide population differently by design**, so neither set is a superset.
48.2% is **DDA-conditioned sequence concordance**: P(we identify a sequence | DDA identifies it).
It is not recall against an exhaustive truth, and the headline in `real-dda-groundtruth.md`
overstated its own caveats. Corrected there.

## The A/B split cannot be cleanly resolved from summary TSVs

| RT window | class A (no spectrum) | class B (emitted, not ID'd) |
|---|---|---|
| 0.5 min | 99.4% | 0.6% |
| 1.0 min | 99.0% | 1.0% |
| 2.0 min | 98.6% | 1.4% |
| **any RT** | **14.8%** | **85.2%** |

The tight-RT numbers are confounded by **run-to-run RT drift** (two injections), the any-RT
numbers by **m/z coincidence** (161,805 emitted spectra over 600 Da, so almost every m/z has a
chance match — floor ~7 spectra per 20 ppm window). The truth is between and **is not
determinable without per-spectrum ID linkage the summary TSVs do not carry.** Reported honestly
rather than by picking a favourable window.

## What DOES resolve cleanly: the misses are weak, long, high-mass peptides

These crossings need no RT matching, so they are unconfounded.

**DDA ID strength — the decisive one.** Misses are DDA's *borderline* identifications:

| DDA hyperscore | MISS % | HIT % | enrichment |
|---|---|---|---|
| 0–15 (weak) | 40.2% | 13.5% | **2.97x** |
| 20–30 (strong) | 17.5% | 42.9% | 0.41x |
| 30+ (very strong) | 4.3% | 11.0% | 0.40x |

MISS median hyperscore **16.0** vs HIT **20.4**. **40% of what we "miss" are DDA's own weak IDs
near its FDR boundary** — not strong peptides we failed on. A large part of the 51.8% gap is DDA
scraping marginal identifications a shallower-sampling method would not confirm either.

**Peptide length and m/z — the genuine extraction gap.**

| length | enrichment | | m/z | enrichment |
|---|---|---|---|---|
| 20–45 | **6.07x** | | 900–1600 | **4.02x** |
| 15–20 | 1.15x | | 700–900 | 0.82x |
| 10–15 | 0.70x | | 500–700 | 0.77x |

Long, high-mass peptides are strongly over-represented in the misses. These carry more fragments
across more charge states and elute broader — harder to assemble into one clean pseudo-spectrum.
This is a real, characterised extraction weakness.

**Charge is NOT a factor.** z2 is 94.4% of misses vs 92.1% of hits — flat. Given how much charge
work this project did, the misses are notably *not* a charge problem.

## What we could NOT measure (stated, not faked)

* **Class C (identified as a different peptide)** — needs per-spectrum ID linkage absent from the
  summary TSVs.
* **True precursor intensity / the low-abundance tail** — no MS1 intensity column wired in; DDA
  hyperscore is only a proxy for ID strength.
* **Cycle density / co-isolation** — not crossed here; RT drift complicates cross-run density
  matching. This is the user's "high-density cycles" question and it remains open.

## Conclusion

The 51.8% "miss" decomposes into two very different populations:
1. **~40% are DDA's weak IDs** — borderline peptides near DDA's own FDR edge, which are not a
   SpeXtractor failure so much as a difference in what the two acquisitions confirm.
2. **A real gap on long, high-mass peptides** (6x at length 20+, 4x at m/z 900+) — the extraction
   genuinely struggles to build clean spectra for large peptides.

Charge, which absorbed most of this project's effort, is not where the peptides are lost.
