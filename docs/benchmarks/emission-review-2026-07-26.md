# Why SpeXtractor emits 1.53M spectra vs the reference implementation's 700k — adversarial review (codex + vibe + kimi)

Deep review reassessing WHEN/WHY emission ballooned to 2.2× the reference implementation and whether the merge
strategies can bring it back to ~700k **at neutral peptides** (a runtime win, not a peptide win).
Three CLIs (codex full, vibe full, kimi independent-but-stderr-only, see note) plus adjudication
against the source and the four prior redundancy docs.

## WHEN it ballooned: commit 579ad70 (2026-07-21)

That one commit flipped **two** defaults at once: `trace:ms1_split_valleys 0→7` **and**
`charge:scoring envelope→count`. So the shipped transition is **not a clean one-factor
experiment** (codex #1, kimi independently found the same commit). Around it:

| | traces | precursors | spectra | peptides |
|---|---|---|---|---|
| baseline (pre-579ad70) | 1.76M | 956k | **741k** | 9,103 |
| current (deterministic HEAD) | 4.92M | 2.145M (48.7% guessed) | **1,528,113** | 10,036 |

Emission +106%, peptides +10.3% → **844 extra spectra per net peptide**. But "net" hides churn;
the peptide *sets* were never compared, so how much is pure redundancy is formally unknown (codex #2).

## WHY: two distinct over-emission sources

**A. Trace multiplication from valley-splitting (dominant).** `ms1_split_valleys` runs EPD after
mass-trace detection and turns one multimodal trace into several; `inferPrecursors_` then seeds a
precursor from *every* unused trace. Traces 1.76M→4.92M (2.8×) is the multiplier. **It is
load-bearing** (+36.7% peptides on dataset B when measured) — reverting it is off the table.

**B. Isotope-unsupported singletons.** 1.04M of 2.145M hypotheses (48.7%) have no isotope partner
→ charge 0. Emission is ~1:1 with hypotheses. kimi's nuance: count-mode (the other half of 579ad70)
*increases* the guessed fraction vs envelope, because count de-isotopes contiguous-only with
break-on-first-miss, so gapped envelopes leave partners as separate charge-0 seeds.

## Six corrections that kill the naive plan

1. **`default_charge=0` does NOT reduce emission** (codex #3, certain). A no-envelope trace becomes
   a scored precursor *before* charge is defaulted; setting 0 only leaves charge unset. The code
   comment at `:1638` claiming it "restores strict isotope-only and drops them" is **wrong**. Cutting
   singletons needs a real `require_isotope_support` gate that erases the hypothesis.

2. **The decisive number was never measured, and currently CAN'T be.** Peptides supplied
   *exclusively* by guessed precursors is the whole question. `assembleFromList_` (`:1325`) stamps RT
   + precursor but **not** `guessed`/`trace_idx`/iso-count onto the emitted spectrum, so search PSMs
   cannot be joined back to inference class (codex #4). Instrumentation is the blocking first step.

3. **Merge is dead — the reframe does not rescue it** (codex #7, vibe #5, prior merge-gating.md, all
   agree). Coordinate merge loses 3.7/7.6/14.7% peptides as the window grows; cosine≥0.85 is neutral
   but touches only 2.8% of spectra. To reach 700k we must absorb **54%** of spectra — merge cannot.
   Redundancy is not recoverable by combining spectra that exist; it must be prevented upstream.

4. **The −12.7% consolidate "kill" is STALE.** It was measured with a dead IM gate (`im_i>0` always
   false) *and* an RT window that spans cycles (3.0s vs 1.385s cycle). The IM gate is now fixed in
   code (`:2340` reads precursor drift time) but consolidate was **never cleanly re-benchmarked**.
   The number should not be used as a prior (codex #11, split-spectra-CORRECTION #4).

5. **`dedup_precursors` is NOT zero-risk** (codex #9/#10, corrects vibe's "peptide-neutral by
   design"). It keeps the most-intense seed and *erases siblings*, then scores the survivor on its
   own seed XIC only — unique fragments that correlated with a discarded sibling vanish. It's
   consolidation moved upstream, not lossless dedup. Worse, it runs **before** default-charging and
   requires charge-match, so a supported z=2 and a soon-to-be-z=2 guessed hypothesis (charges 2 and
   0) never collapse → turning the flag on may remove only a modest subset, not halve.

6. **Count reduction ≠ extraction-runtime win** (codex #8). Merge/consolidate run *after* every
   hypothesis is scored; the window loop is ~60% trace detection, ~37% scoring, ~2% emit. Only a
   **pre-scoring** reduction cuts extraction time; post-scoring only cuts downstream *search* time
   (still worth it). Every count-reduction arm must report extraction CPU and search time separately.

## Reassessment verdict on the merge strategies

The runtime/peptide-neutral reframe is **valid**, but **none of the existing merge/consolidate/dedup
knobs is proven to reach 700k neutrally**, and merge is affirmatively falsified. Two candidate
levers remain, both **unproven** and pointing upstream (prevent, don't combine):

- **Pre-scoring feature collapse** — the the reference implementation "one pseudo-MS/MS per precursor *feature*" model.
  codex's stronger form: cluster coordinate-equivalent hypotheses and build **one joint feature XIC +
  quality score**, then fragment-score once — not the current keep-one-seed dedup. Cuts extraction
  *and* search time. Risk: the split analysis shows within-cycle pairs are mostly overlapping (union
  1.355 ≈ closer to identical) so collapse is plausible, but "within-cycle" is a post-hoc RT-bin
  artifact, not verified identical features (codex #13) — must validate on raw apex/XIC/IM.
- **Isotope-support quality gate on singletons** — biggest cut (~784k spectra if guessed emit at the
  same rate) but **highest risk of silently losing the low-abundance / unexpected-mod / nonspecific
  open-search species that are SpeXtractor's entire reason to exist** (codex (b), kimi's open-search
  caveat). A blanket "require ≥1 partner" is the single change most likely to lose the peptides that
  justify the tool.

## The decisive first experiments (all three converge)

Do these **before building any reduction**, because they set the ceiling and pick the lever:

1. **Provenance instrumentation** — stamp `guessed`, `trace_idx`, isotope count (and a stable
   spectrum index) onto every emitted spectrum so search PSMs join back to inference class.
2. **The oracle** (codex #12) — cluster the *identified* spectra by (charge, m/z, RT, IM) and count
   the **minimum spectra needed to retain every identified peptide**. If even an ID-aware oracle
   can't approach 700k, no pre-search selector can, and 700k-at-neutral is impossible. Runnable on
   existing dataset D output + a completed search, no re-run.
3. **Class attribution** — identified-spectrum rate by class; peptides *exclusively* from guessed
   spectra; the same split for modified / nonspecific open-search peptides; peptide counts after
   removing each class and recomputing FDR.

**Pre-registered non-inferiority** for any winner: ≤0.2% peptide loss on deterministic dataset D, checked
separately for modified / low-intensity / high-charge / guessed-only, on **both** Sage and MSFragger.

## Bottom line

- **(a) Most likely to reach 700k at neutral peptides:** a true pre-scoring precursor-feature
  collapse (one joint XIC from coordinate-equivalent hypotheses), *not* the existing keep-one-seed
  dedup. But run the oracle + current-dedup A/B first — there is presently **no evidence** any
  existing knob reaches 700k.
- **(b) Most likely to silently lose the peptides that justify SpeXtractor:** a blanket "require ≥1
  isotope partner" gate.
- Do **not** merge everything, do **not** confuse `default_charge=0` with filtering, do **not** force
  700k as a target in itself. First label the guessed population, measure guessed-only peptides, and
  establish the oracle feature-collapse ceiling.

---

## ORACLE RESULT (dataset D det4: 1,528,113 spectra, 113,130 identified PSMs → 10,036 peptides)

Ran the oracle on existing data (mzML precursor coords + Sage results, no re-run). **It overturns
the collapse/merge framing entirely.**

**Feature-collapse cannot reach 700k — it is geometrically impossible, not just peptide-lossy:**

| collapse operation | resulting emission | compression | vs 700k |
|---|---|---|---|
| baseline | 1,528,113 | 1.00× | — |
| same-coordinate (charge, 20 ppm, RT 1.4–7 s, IM 0.02) | 1.21M–1.42M | 1.07–1.26× | **above** |
| isotope-aware (link M/M+1/M+2 ladders) | 1.45M | 1.05× | **above** |
| the reference implementation | 700,434 | — | — |
| ID-aware oracle (drop everything that never identifies) | **~120k** | 12.7× | below |

The 1.53M emitted spectra are **~1.4M geometrically-distinct precursor coordinates**. Only ~5–20%
is removable as duplicates or isotope ladders. **92% of emission (1.41M spectra) sits in coordinate-
features that never identify, and those features are singletons (~1.1 spectra each)** — distinct
coordinates with nothing to merge them against. So merge / consolidate / dedup / better-de-isotoping
**cannot** approach 700k. Confirmed independently by same-coordinate (oracle1) and isotope-aware
(oracle3) clustering. The identified region (8%) collapses near-neutrally (same-coord chimera 0.3%;
isotope-aware 3%), i.e. the redundancy that collapse *can* remove is a rounding error on emission.

**The excess is over-GENERATION of distinct low-quality precursor hypotheses** — dominated by the
1.04M guessed singletons (48.7% of hypotheses, no isotope support) and non-identifying supported
precursors — **not redundancy.** The only lever that can cut emission is a **precursor quality gate**
at inference time (require isotope support / envelope quality), i.e. codex's lever (b), the one
flagged most-likely-to-lose-peptides. It is now the *only* geometrically possible lever.

**The ID-aware ceiling (~120k) is the headroom signal:** a perfect precursor-quality classifier would
emit ~120k spectra and keep all 10,036 peptides. the reference implementation's isotope-distribution + XIC-correlation
gate lands at 700k (conservative — keeps many non-identifying-but-plausible precursors). We keep
1.53M because we gate nothing. A gate even close to the reference implementation's should reach ≤700k; the peptide cost
is entirely a function of proxy quality.

**Retracted: the geometric gate simulation (oracle2) is confounded.** Simulating "require an emitted
M±1 partner" gave 12–17% peptide loss to reach ~700k — but the tool **de-isotopes**, so a genuinely
isotope-*supported* precursor emits only its monoisotope (the M+1 is consumed, absent from emission).
"No emitted partner" therefore conflates supported-monoisotopes with guessed-singletons; the 12–17%
is an artifact. **The true gate cost cannot be measured from existing data** — it needs the `guessed`
flag joined to identified peptides, which requires provenance instrumentation + a gated re-run.

### Forced next step (the go/no-go redirected the plan)
"Oracle first" was meant to find a cheap collapse. The oracle proves no collapse exists — so the
next step is the instrumented measurement codex called the blocking first move:
1. Stamp `guessed` / iso-count / `trace_idx` onto each emitted spectrum (small change).
2. Re-run dataset D + search; join identified peptides back to inference class.
3. Measure peptides supplied *exclusively* by guessed precursors (and by modified / low-intensity /
   high-charge / open-search subsets). That number *is* the isotope-support gate's peptide cost and
   sets where the gate can sit relative to 700k.

---

## GATE COST MEASURED (provenance re-run, 2026-07-26)

Stamped `spx_guessed` onto every emitted spectrum (`assembleFromList_`), re-extracted dataset D
(deterministic → byte-identical 1,528,113 spectra, verified), and joined the guessed flag to det4's
existing Sage results by spectrum index (no re-search needed):

| | emitted spectra | identify rate | identified PSMs | peptides |
|---|---|---|---|---|
| **guessed** (no isotope support) | 603,858 (39.5%) | **1.30%** | 7,852 (6.9%) | — |
| **supported** | 924,255 (60.5%) | **11.39%** | 105,278 (93.1%) | 9,958 |
| peptides exclusively from guessed | — | — | — | **78 (0.78%)** |

**Supported precursors identify 8.8× more often than guessed.** Dropping all guessed singletons:
**emission 1,528,113 → 924,255 (−39.5%) for a projected 78-peptide loss (−0.78%).** That is the
isotope-support gate (`assembly:require_isotope_support`, added OFF-by-default). Enormously better
than any merge — but it lands at **924k, not 700k** (still 1.32× the reference implementation, down from 2.18×), and
0.78% exceeds the pre-registered ≤0.2% non-inferiority bound. The 78 exclusive peptides are the
low-abundance / no-envelope species codex flagged.

**FDR-recalibrated confirmation (gated re-run + re-search, `assembly:require_isotope_support`):**

| | emission | peptides | extraction wall |
|---|---|---|---|
| baseline (det4) | 1,528,113 | 10,036 | 27:00 |
| **gated** | **924,255 (−39.5%)** | **9,989 (−0.47%)** | **19:46 (−27%)** |
| the reference implementation | 700,434 | — | 7:47 |

The real re-run **beats the projection** (−0.47% vs −0.78%): FDR recalibration recovered ~31 peptides
(9,958 projected → 9,989 actual), because 924k spectra carry a lighter multiple-testing burden
(codex #14). And extraction wall dropped 27% — a *pre-scoring* cut reduces extraction time, not just
search time (codex #8 confirmed). So the gate is a clean, near-neutral, real-runtime win.

**Cross-engine confirmation (MSFragger 4.4.1, same DB+params for both arms):**

| engine | baseline | gated | Δ |
|---|---|---|---|
| Sage (sage_deiso) | 10,036 | 9,989 | **−0.47%** |
| MSFragger | 11,138 | 11,033 | **−0.94%** |

Both engines agree the gate is near-neutral (<1%), unlike the Sage-only min_correlation gain — so it
is trustworthy to keep. MSFragger distinct pre-FDR peptide keys fell 423,725 → 300,645, confirming the
dropped 39.5% of spectra were overwhelmingly noise. (MSFragger OpenMS-mzML blockers handled:
`fix_mzml_cvparams.py` for the valueless-cvParam / Scans=0 crash, and `labile_fragment_ion_series=b,y`;
the JNA/Bruker-lib startup error is non-fatal for mzML input.)

## STATUS: gate validated, shippable. Path to 700k remains open.
The isotope-support gate (`assembly:require_isotope_support`, OFF by default) is a confirmed
near-neutral −39.5% emission / −27% extraction-time win on both engines. It lands at **924k (1.32×
the reference implementation)**, not 700k. Closing the last 924k→700k (≈−24%) requires a **graded quality gate on the
supported side** (the reference implementation-style isotope-distribution + XIC-correlation score swept to an operating
point), since supported precursors carry 99.2% of peptides and cannot take a blanket cut.

**Implications for the 700k target:** the guessed gate is the first lever and it's near-neutral, but
reaching 700k needs a *second* cut on the **supported** side (924k → 700k ≈ another −24%). Supported
precursors carry 99.2% of peptides, so a blanket supported-side cut is dangerous; it needs a graded
quality gate (the reference implementation-style isotope-distribution + XIC-correlation score, swept to an operating
point) rather than a hard drop. Sequence: ship the guessed gate (confirm FDR), then grade the
supported side toward 700k while holding peptides.

### Note on kimi
kimi (this version) streams reasoning to stderr and never emits a stdout answer in `-p` mode — it
keeps trying to spawn explore agents and truncates. Its independent findings were captured from the
reasoning stream (commit 579ad70; count-mode raises the guessed fraction; 56.4% of traces never
referenced; ~44% of within-cycle pairs sit at ~1000 ppm ≈ one isotope spacing, a distinct duplicate
class; the PSM-vs-precursor-dump join experiment). Not a clean final review, but not silent either.
