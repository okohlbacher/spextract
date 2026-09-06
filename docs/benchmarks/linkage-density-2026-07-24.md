# Per-spectrum ID linkage, density, and the long/high-mass gap (dataset D, 2026-07-24)

Follow-up to `miss-analysis-2026-07-23.md`. That analysis listed four next steps; this doc executes
them. Three of the four are **measurements**; the fourth (#3) is a decision. The headline is that
the linkage instrument **falsified the two extraction hypotheses the miss analysis proposed** — the
levers I named were wrong, and the data says where the peptides actually go.

## The instrument (#2): the join nobody had wired

Sage's `scannr` **is** the mzML spectrum id **is** the emission index. So one streaming pass over
`pseudo.mzML` yields the FULL emitted population (identified + not) and the Sage TSV labels the
identified subset by `scannr`. That single linked table is the A/B instrument the miss analysis said
was missing — and because it works INTERNALLY (our own run), it needs no cross-run RT matching, which
resolves the A/B confound the DDA comparison could not (two injections, RT drift).

* Script: `bench/linkage_density.py` (streams the 14.7 GB mzML once, caches the population to
  `pop_s30.tsv`, then all crosses are seconds).
* Population: **1,528,147 emitted spectra**; **113,495 identified** at `peptide_q<=0.01` (**7.4%**).
  Baseline Sage: **9,966 unique peptides**, 90,159 PSMs at 1% FDR.

**Join integrity** (checked, because the whole doc rests on it): of 113,495 identified spectra,
**0 have a `scannr` absent from the population** — the join is exact, no off-by-one or re-index.
Only **0.4%** of emitted spectra fall below Sage's `min_peaks=6` (structurally unscorable), and
removing them moves the ID rate 7.4%→7.5% — so "over-emission" is not an artefact of counting
unscorable junk. The emitted spectra are *full*, not empty (73.5% at ≥500 peaks; only 6,743 below
6 peaks); they carry peaks but do not identify.

Honest limits (stated up front — the MS1 funnel was retracted for skipping these):
* "identified" = a Sage PSM at `peptide_q<=0.01` on that spectrum id. A spectrum with only a
  sub-FDR PSM counts as NOT identified. Correct: we emitted it, we could not confirm it.
* This is internal to our dataset D run. It does **not** relabel the external DDA miss set.

## #1 — the long/high-mass hypothesis is FALSIFIED internally

The miss analysis proposed: *high-mass peptides are truncated by the 500-peak cap; raise
`max_fragments`.* The linkage kills it three ways.

**ID rate by precursor m/z** — 900-1200 is our BEST band, not a gap:

| m/z | emitted | id'd | ID rate |
|---|---|---|---|
| 300-500 | 445,307 | 8,262 | **1.9%** |
| 500-700 | 395,504 | 42,068 | 10.6% |
| 700-900 | 451,210 | 40,487 | 9.0% |
| 900-1200 | 203,083 | 21,655 | **10.7%** |
| 1200-2000 | 33,043 | 1,023 | **3.1%** |

The ID rate is flat across the core (500-1200 ≈ 9-11%) and drops only at the **extremes**
(300-500 and >1200), and the >1200 tail is only 2% of emitted spectra. There is no high-mass gap
internally; long peptides do land at high m/z (identified-peptide length rises 10→12→15→**19**→23 aa
across those bands) and they identify at the core rate.

**The 500-peak cap is saturated everywhere**, so high mass is not specially truncated:

| m/z band | npeaks median | at cap (≥499) |
|---|---|---|
| 500-700 | 500 | 69.9% |
| 900-1200 | 500 | 74.8% |

**Sage reads only the top-150 peaks** (`max_peaks=150`). So the 500-peak spectra are already
downsampled by the engine; raising `max_fragments` above 500 cannot add IDs, it only adds noise the
engine discards. The peak-count → ID-rate curve is monotone (0.0% at <10 peaks → 8.2% at the cap) but
**confounded**: rich spectra come from abundant, co-eluting peptides, so the correlation is not causal
and cannot be read as "more peaks help".

**Redirect:** the lever is not fragment *count* but the *cleanliness of the top-150* — whether they
are polluted by fragments belonging to a co-isolated precursor. That is the corrected experiment below.

## #4 — the "high-density cycles hurt us" hypothesis is FALSIFIED

Co-isolation density proxy: emitted precursors sharing an acquisition cell (RT cycle ±0.03 min × IM
band ±0.03 × ~25 Th isolation tile). 354,916 cells, max occupancy 49. **ID rate RISES with density**,
and it survives fixing the m/z band (which strips most of the elution confound):

| density (siblings) | all m/z | m/z 700-1200 | m/z 500-700 |
|---|---|---|---|
| 0 (alone) | 3.1% | 4.1% | 3.8% |
| 1-2 | 4.8% | 6.1% | 6.1% |
| 3-5 | 6.4% | 8.0% | 8.5% |
| 6-10 | 8.3% | 10.9% | 11.2% |
| 11-20 | 9.9% | 17.8% | 13.7% |
| 21+ | 12.0% | 26.0% | 17.8% |

The user's "high-density cycles" question is answered: **there is no density penalty** — the opposite.
Dense cells are mid-gradient, real-peptide-rich regions; the **sparse/isolated** spectra are the worst
(3.1%). Caveat: even within a fixed m/z band the proxy is partly an elution proxy, so this is not
proof co-isolation is harmless — but no penalty is visible on any axis measured.

## The unifying internal finding: over-emission

1.53M emitted spectra for 9,966 peptides (~150 emitted per identified peptide; ~13.5 per PSM), 74%
stuffed to a 500-peak cap the engine downsamples to 150, 92.6% never identified. The identifiable
signal concentrates in rich mid-gradient cells; the tool sprays low-quality spectra broadly. This is
the same redundancy `open_bench.py` pre-registered (4.6 PSMs/peptide) seen from the emission side.

## #3 — charge is deprioritised

Both the external miss set (charge FLAT) and this internal linkage (ID rate moves on peak-count and
density, not charge) say charge is not where peptides are lost. The joint charge model is shelved.
See memory `charge-not-the-bottleneck`.

## Corrected #1 experiment (de-chimerize the top-150): base done, treatment INTERRUPTED

`bench/plan_dechim.yaml` via the harness: `base` vs `assembly:rp_max=1` (keep each fragment in only
its single best-correlating precursor). If IDs rise, cross-precursor sharing polluted the top-150; if
they fall/flat, redundancy is load-bearing and de-chimerization is another falsified lever.

* **base arm: complete.** Extract 1863 s, 81 GB peak; Sage 90,162 PSMs — a faithful reproduction of
  the shipped-default baseline (reference 90,159), which validates the harness plan.
* **rpmax1 arm: interrupted at 99% of extraction** when node-1/06 were taken down for
  maintenance mid-run. Sage never ran on it. The number is pending the nodes' return; on resume,
  re-run the `rpmax1` arm against the SAME binary sha recorded in `base`'s manifest (the FAILURE-2
  one-invocation guarantee is preserved by asserting binary identity, since `base`'s mzML already
  exists).

<!-- RPMAX1_RESULT -->

## Adversarial review outcome (codex + vibe; kimi silent) — SUBSTANTIAL CORRECTIONS

Two independent reviews converged on the same headline flaw and I verified the sharpest claims
against the source. The **measurements stand; several interpretive claims do not.** Accepted:

1. **"Falsifies the external high-mass hypothesis" is a CATEGORY ERROR (both reviewers' #1 flaw).**
   The internal instrument measures P(ID | *emitted* spectrum); it is blind by construction to
   class-A "never emitted" peptides, which is exactly what the external 4.02× high-mass miss
   enrichment is about. RETRACTED wording: emitted high-mass spectra identify at the core rate, so
   IF a high-mass gap exists it is in EMISSION (class A) or in database censoring (Sage caps 5000 Da
   / length 35), NOT in per-spectrum ID of what we emit. The two analyses are not commensurable;
   neither overrides the other.
2. **Config mixing (codex #15).** The linkage ran on `rebase_s30` = `sage_closed_s30.json`
   (9,966 peptides / 90,159 PSMs). The FROZEN dataset D baseline is **10,072 peptides** under
   `sage_deiso.json` (`docs/dataset D-BASELINE.md`). Calling 90,159 the "reference" for a `sage_deiso`
   base arm was wrong. All linkage rates are `sage_closed_s30.json`; re-run on the base arm's
   `sage_deiso` pair when the cluster returns.
3. **`max_fragments` falsification not airtight (codex #6).** SpeXtractor ranks the top-500 by
   CORRELATION (`assembleFromList_`, `frag_scores`=Pearson `c`); Sage ranks its top-150 by
   INTENSITY. Different criteria, so raising `max_fragments` CAN change Sage's top-150. Corrected:
   the effect is untested and directionally unknown (added peaks are low-correlation = likely
   noise), not zero.
4. **Charge assertion unsupported (codex #16 / vibe #8).** `linkage_density.py` parses no charge;
   "internal linkage says not charge" is retracted. Charge deprioritisation rests ONLY on the
   external prevalence-flat signal, which has low power on the z2-dominated data — absence of
   evidence, not a measured null.
5. **Density proxy is endogenous (codex #9 / vibe #3).** Occupancy counts SpeXtractor's own emitted
   hypotheses, a collider on abundance; fixing m/z removes little of it. Corrected: "no density
   penalty" → "no penalty visible in a proxy that cannot isolate one." Real MS1 feature density +
   within-feature matching is needed to make the causal claim.
6. **peptide_q vs spectrum_q mixing (codex #3).** "Identified" uses `peptide_q<=0.01` (113,495);
   the 90,159 is `spectrum_q`. Rates should be reported under both. Arithmetic: "~13.5/PSM" mixed
   denominators — it is 1,528,147/113,495; per spectrum_q PSM it is 16.95.
7. **Length claim is survivorship-biased (codex #7)**; **"over-emission" needs a the reference implementation
   spectrum-level baseline to be a finding rather than a label (both #4)**; **no chance floors /
   cluster-bootstrap CIs on any rate (codex #18 / vibe #7)** — the exact MS1-funnel failure mode.

Where I push back: the "join is globally offset" BLOCKER (both #1) shows my *existence* check was
insufficient, not that the join is wrong — my parser reads the same `id="spectrum=N"` string Sage
reads from the same file. The realistic risk is Sage using a 1-based scan counter; I will assert
value-equality on sampled spectra when the cluster returns rather than claim it proven now. The
"rp_max=1 guaranteed no-op" (vibe #5) is too strong — codex's own #14 shows rp_max edits
fragment membership PRE-cap and can remove high-intensity top-150 peaks — but both are right that a
FLAT result would be ambiguous (no pollution vs lost-but-compensated), so the experiment needs a
top-150 Jaccard delta computed before searching.

**Net:** the linkage *instrument* is sound and the *measurements* (m/z ID-rate of emitted spectra,
cap saturation, density correlation) hold as descriptive facts about emitted spectra. The
*interpretations* — "falsifies the high-mass gap", "no density penalty", "charge not the axis" —
overreached and are downgraded to open questions requiring the controls above.

## The missing baseline, measured (answers the "over-emission" BLOCKER) — 2026-07-24, node-2

Both reviewers made this a BLOCKER: 7.4% has no meaning without a reference extractor's
spectrum-level ID rate on the same data. Measured on node-2 (node-1/06 in maintenance) by
running the SAME `sage_deiso.json` on the reference implementation's dataset D mzML from `/ceph`:

| | SpeXtractor | the reference implementation |
|---|---|---|
| emitted spectra | 1,528,147 | **700,434** |
| identified spectra (peptide_q ≤ .01) | 113,495 | 26,986 |
| **spectrum-level ID rate** | **7.4%** | **3.85%** (spectrum_q: 2.75%) |
| PSMs @1% FDR | ~90k | 19,458 |
| **peptides @1% FDR** | 10,072 | **11,517** |
| identified spectra per peptide | ~11.3 | **2.34** |

**"Over-emission = we spray low-quality spectra" is FALSIFIED by its own baseline.** Our
per-spectrum ID rate is roughly **twice** the reference implementation's; the emitted spectra are not junk. the reference implementation
also caps at 500 peaks (`defaultArrayLength="500"`), so cap saturation is not ours either.

**The real defect is redundancy vs coverage.** We emit 2.18× more spectra and identify **1,445
fewer peptides** (10,072 vs 11,517, −12.6%). Our identified spectra pile onto ~11.3 per peptide
against the reference implementation's 2.34 — the emission-side view of the documented 8.98 vs 1.66 PSMs/peptide.
So the deficit is **breadth**: peptides the reference implementation covers that we never do. That is precisely the
class-A/emission axis the reviewers said the internal instrument is blind to — which makes
"which peptides does the reference implementation get that we don't, and are they emitted at all" the single
highest-value next measurement, replacing both falsified levers.

Caveats: emitted-spectrum counts and the peptide comparison are config-matched and solid
(both `sage_deiso.json`; peptide reference is the frozen 10,072 baseline). The 7.4% figure is
still from the `sage_closed_s30.json` linkage run, so the ID-RATE ratio is preliminary until the
base arm re-run under `sage_deiso` lands; the direction (ours ≫ the reference implementation) is not in doubt at
this margin.

## What we actually lose: the coverage gap characterised (node-2, 2026-07-24)

The clean comparison the project never had: the reference implementation and SpeXtractor pseudo-spectra from the **same
dataset D acquisition**, searched with the **same `sage_deiso.json`** against the same FASTA. Extraction
is the only variable — no second injection, no RT drift, no DDA-conditioning. `bench/coverage.py`.

| | peptides |
|---|---|
| shared | 8,481 |
| **the reference implementation only (we lose)** | **3,036** |
| ours only (the reference implementation loses) | 1,545 |
| union | 13,062 |

We are **not a subset**: we find 1,545 peptides the reference implementation does not. Coverage of the union: ours
76.8%, the reference implementation 88.2%.

### The distinguishing property is ABUNDANCE, and nothing structural

MISS vs SHARED medians (the reference implementation's own best PSM for each of its peptides):

| property | MISS | SHARED | ratio |
|---|---|---|---|
| **precursor intensity** | **5,750** | **18,650** | **0.31** |
| ms2_intensity | 1.26e4 | 3.57e4 | 0.35 |
| matched_intensity_pct | 8.45% | 17.41% | 0.49 |
| matched_peaks | 9 | 12 | 0.75 |
| hyperscore (the reference implementation's own) | 26.9 | 35.1 | 0.77 |
| spectrum TIC | 2.76e5 | 3.65e5 | 0.76 |
| peptide length | 13 | 14 | 0.93 |
| precursor m/z | 650 | 685 | 0.95 |
| ion mobility | 0.962 | 0.985 | 0.98 |
| **charge** | 2 | 2 | **1.00** |

Precursor-intensity enrichment is monotonic: **9.4×** below 500, 5.7× at 1–2.5k, **0.40×** above
10k (70.3% of SHARED is >10k vs 28.4% of MISS). Everything structural is flat. Charge is exactly
1.00 — a third independent confirmation it is not the axis. Length and m/z lean SHORT and LOW
(7-10 aa 1.70×; m/z 300-500 1.43×; 900-1200 0.61×), the opposite of the DDA-conditioned analysis,
as expected for a non-DDA-conditioned comparison. Part of the gap is the reference implementation's own FDR edge
(6.2× enriched at hyperscore 15-25, 0.06× above 40).

### Class A/B: we DO emit them — the gap is identification, not emission

Well-posed here because both sides share one RT axis. Every number beside a decoy-coordinate
floor (same RT/IM, non-physical m/z offset), and a positive control that must be ~all-B:

| RT window | A never emitted | B emitted, not ID'd | A% | decoy floor for B |
|---|---|---|---|---|
| 0.2 min | 320 | 2,716 | **10.5%** | 16.0% |
| 0.5 min | 152 | 2,884 | 5.0% | 27.4% |
| 1.0 min | 74 | 2,962 | 2.4% | 37.3% |

**Positive control (SHARED, must be ~all B): 0.9% A / 99.1% B at 0.5 min — passes.** At the
tightest window **89.5% of missed peptides are class B** against a 16.0% coincidence floor.

So the "we never emit high-value peptides" hypothesis is ALSO falsified. We emit a spectrum at
~90% of the missed coordinates; we fail to *identify* it. Caveat: class B means an emission within
20 ppm / RT / IM of the coordinate, not proof that the spectrum is of that peptide — the decoy
floor bounds the coincidence, it does not identify the content.

**Synthesis.** We are not less accurate than the reference implementation — our spectra identify at ~2× its
per-spectrum rate. We are less **sensitive**: at faint precursors we emit a spectrum that is too
thin (8.5% matched intensity, 9 matched peaks) to clear the identification threshold, while
spending emission budget re-deriving spectra for peptides already covered (11.3 identified spectra
per peptide vs 2.34). The lever is the low-abundance sensitivity floor — trace-detection
thresholds and how weak precursor hypotheses are admitted and assembled — NOT emission coverage,
NOT charge, NOT fragment count, NOT de-chimerization.

### #1 experiment, settled: de-chimerization is strongly harmful

`base` vs `assembly:rp_max=1`, one binary, one invocation (node-2):

| | base | rp_max=1 | Δ |
|---|---|---|---|
| peptides | 10,124 | 4,819 | **−52.4%** |
| PSMs | 90,154 | 19,542 | −78.3% |
| emitted spectra | 1,528,083 | 1,427,890 | −6.6% |
| wall | 1,353 s | 3,270 s | 2.4× slower |

Spectra largely survive (−6.6%) but stop being identifiable (PSMs −78%). Cross-precursor fragment
sharing is **load-bearing**. Per the review's warning, the WRONG reading is "top-150 pollution is
absent"; the right one is that winner-take-all assignment strips real signal along with any
pollution. Lever closed.

Two bugs were caught by controls during this analysis and fixed before reporting: a native-id
parse that collapsed 700,434 the reference implementation spectra into 1, and FAILURE 8 again — our mzML writes
`scan start time` in SECONDS (`unitName="second"`) while Sage reports minutes, so the first A/B
attempt compared seconds to minutes and the positive control read a physically impossible 100%
class A. `coverage.py` now reads the unit from the CV accessor and aborts if the RT ranges fail to
overlap.

## Sensitivity-floor sweep: the co-elution gate is NOT the floor (5th lever closed)

The coverage analysis said the missed peptides are faint, already emitted (~90% class B), and
carry thin spectra (9 matched peaks / 8.45% matched intensity vs 12 / 17.41%). Prime suspect:
`gate:min_correlation=0.3` admits a fragment only if its XIC correlates >=0.3 with the precursor's,
and Pearson correlation is attenuated toward 0 as S/N falls — a gate structurally biased against
low-abundance precursors. Pre-registered in `bench/plan_sens.yaml`: **if peptides do not rise at
0.20 or 0.12, the gate is not the sensitivity floor.**

| `gate:min_correlation` | peptides | Δ vs base | PSMs | emitted spectra | wall |
|---|---|---|---|---|---|
| **0.30 (base)** | **10,138** | — | 90,206 | 1,528,093 | 1456 s |
| 0.20 | 9,829 | **−3.0%** | 90,257 | 1,532,398 | 1463 s |
| 0.12 | 9,696 | **−4.4%** | 90,102 | 1,533,107 | 1471 s |

Monotone, three points, all far outside the 0.17% replicate spread. The base arm reproduced the
independent `plan_dechim` base to **+0.14%** (10,138 vs 10,124), so the reference is sound.

**Criterion applied: peptides FELL. The co-elution gate joins the closed levers** (`max_fragments`,
density, charge, de-chimerization, co-elution gate = five).

Mechanism, and it is consistent with everything else measured: PSMs are flat to ±0.2% and emitted
spectra move only +0.3%, so loosening the gate does not fail to find things — it **converts unique
peptides into redundant re-identifications of peptides already covered**. Admitted low-correlation
fragments are noise, and because Sage scores only the **top-150 by intensity**, intense noise peaks
displace true fragments from the window that actually gets scored. Loosening any fragment-admission
knob therefore costs peptides, which is the same lesson `rp_max=1` taught from the opposite side.

**Useful residue: the derivative points the other way.** Peptides rise monotonically WITH gate
strictness across the whole tested range (9,696 → 9,829 → 10,138), so 0.3 may itself be too loose.
`min_correlation` 0.40 / 0.50 is the evidence-motivated next test — grounded in three measured
points, not a hunch. Untested; do not quote a benefit until measured.

### STRICT direction: the shipped default 0.3 is TOO LOOSE — first real gain of the session

The loosening sweep's derivative pointed the other way, so the same sweep was extended past the
default. It **wins, and has not turned over**:

| `gate:min_correlation` | peptides | Δ vs base | PSMs | emitted spectra |
|---|---|---|---|---|
| 0.12 | 9,696 | −4.0% | 90,102 | 1,533,107 |
| 0.20 | 9,829 | −2.7% | 90,257 | 1,532,398 |
| **0.30 (shipped)** | **~10,122** (3 reps: 10,124/10,138/10,104) | — | 90,176 | 1,528,014 |
| 0.40 | 10,380 | **+2.7%** | 89,819 | 1,514,995 |
| **0.50** | **10,657** | **+5.5%** | 88,130 | 1,485,454 |

Base replicate spread is **0.34%**, so +2.7% and +5.5% are 8–16× noise. **PSMs FALL as peptides
RISE** (90,176 → 88,130) — redundancy converting into coverage, which is the direction that matters
given we over-produce PSMs ~8.9/peptide against the reference implementation's 1.69. PSMs/peptide improves 8.92 → 8.27.

Against the reference implementation's 11,517 on the same acquisition, this narrows the peptide gap from **−12.0% to
−7.5%** — roughly 37% of the gap closed by one parameter.

Extended (`bench/plan_strict2.yaml`) — **the curve TURNS OVER, so there is a real interior optimum
at 0.60**, and the pre-registered "unbounded → pick nothing" clause does not fire:

| `gate:min_correlation` | peptides | Δ vs base | PSMs | PSMs/peptide | emitted |
|---|---|---|---|---|---|
| 0.12 | 9,696 | −4.2% | 90,102 | 9.29 | 1,533,107 |
| 0.20 | 9,829 | −2.9% | 90,257 | 9.18 | 1,532,398 |
| **0.30 (shipped)** | **10,125** (4 reps: 10,124/10,138/10,104/10,135) | — | 90,214 | 8.90 | 1,528,014 |
| 0.40 | 10,380 | +2.5% | 89,819 | 8.65 | 1,514,995 |
| 0.50 | 10,657 | +5.3% | 88,130 | 8.27 | 1,485,454 |
| **0.60** | **10,903** | **+7.7%** ← peak | 84,288 | **7.73** | 1,423,257 |
| 0.70 | 10,675 | +5.4% | 76,813 | 7.20 | 1,302,939 |

Unimodal, peak at 0.60, +7.7% over the shipped default against a **0.34% replicate spread** (4
independent base runs across two harness invocations). Emitted spectra fall monotonically with
strictness (1.528M → 1.303M), so past 0.60 real spectra start dropping below `min_fragments` —
a plausible reason the curve turns.

Against the reference implementation's 11,517 on the same acquisition this narrows the peptide gap from **−12.0% to
−5.3%** — ~56% of the gap closed by one parameter, while PSMs FALL 90,214 → 84,288.

### RETRACTED: the +7.7% is an FDR RECALIBRATION ARTEFACT, not a sensitivity gain

All three reviewers (codex, vibe, kimi) independently made this their #1 BLOCKER, and a direct
test on data already on disk **confirms it**. Sage derives `peptide_q` from its OWN decoys per run,
so `q<=0.01` is a different score cut in every arm. Applying ONE FIXED discriminant threshold across
arms (`bench/fdrcheck.py`) removes that confound:

| arm | q-based Δ | fixed-score Δ (s≥1.0 / 1.5 / 2.0 / 2.5) | decoys | own q≤.01 threshold |
|---|---|---|---|---|
| base | — | — | 217,058 | 2.1711 |
| corr040 | +2.4% | **+1.3 / +3.3 / +14.3 / +13.6%** | 219,649 | 2.2034 |
| corr050 | +5.1% | −0.2 / −7.0 / −7.7 / −11.5% | 217,351 | 2.1218 |
| **corr060** | **+7.6%** | **−5.6 / −14.4 / −10.2 / −17.7%** | 203,513 | **2.1031** |
| corr070 | +5.3% | −19.3 / −30.3 / −13.5 / −17.1% | 171,204 | 2.1010 |

At **every** fixed score cut corr060 has FEWER peptides than base. The apparent gain comes entirely
from the 1%-FDR threshold sliding down (2.1711 → 2.1031) as decoys drop (217k → 203k): tightening
the gate deletes spectra, shifts the target-decoy balance, and lets more peptides clear q≤0.01 with
*less* real signal. **The "+7.7%" and "56% of the gap closed" claims are withdrawn.**

Only **corr040** survives both criteria (+2.4% q-based AND positive at all four fixed cuts) and is
the one arm still worth pursuing — modestly, and only under entrapment FDR.

Caveat on the test itself: `sage_discriminant_score` is a per-run LDA output, so a "fixed score" is
not perfectly commensurable either. The direction is unambiguous (decoys down, threshold down,
q-count up, fixed-score count down) but the decisive experiment remains **entrapment FDR**, which
all three reviewers named as the single highest-value missing measurement — and which this project
has never run despite quoting "1% FDR" throughout.

Compounding methodological faults the reviewers caught, all valid:
* **Winner's curse.** Seven swept points, optimum chosen on the same sample used to evaluate it, and
  the 0.60/0.70/0.80 pre-registration was written *after* watching 0.12→0.50 rise. Only `base` has
  replicates (n=4); every swept arm is n=1.
* **Wrong engine for the headline.** `docs/dataset D-BASELINE.md` already records the MSFragger gap as
  **−19.5%** (10,470 vs 13,014), not the −12% I quoted from Sage. Project policy says Sage flatters
  SpeXtractor; I published the flattered number.
* **`rp_max=1` re-tested a closed lever.** `dataset D-BASELINE.md` lists `competitive` and `rp_max` 2/4/8
  as already falsified. ~1.5 h spent re-deriving a known result.

### MS1 sensitivity arm: BLOCKED BY IMPLEMENTATION, not answered

`ms1sens` (`trace:noise_threshold_int` 100→10, `trace:ms1_chrom_peak_snr` 3.0→1.0) ran **40 min
with zero log progress at 1.5 of 128 cores** (RSS 79 GB) while a whole base run takes 24 min. Cause
is architectural and already known in this project: OpenMS `MassTraceDetection` has **no OpenMP** —
a serial greedy pass over an intensity-sorted apex list — and lowering the MS1 thresholds explodes
that list, so runtime blows up superlinearly. `perf:trace_bands=12` did not rescue it (we would see
~12 cores, not 1.5). Killed; **no metrics were recorded, so nothing garbage entered the results.**

This is a result about the implementation, NOT about MS1 sensitivity: whether a lower MS1 floor
would recover faint precursors remains **unmeasured**. To test it at all, either MassTraceDetection
must be parallelised or a milder setting (e.g. noise 50 / snr 2.0) must be used. Do not record the
MS1 direction as falsified — it was never measured.

## UN-RETRACTED: entrapment FDR says the gate gain is REAL (2026-07-24, node-2)

The retraction above was **premature**, and the instrument that settled it is the one all three
reviewers demanded. Entrapment search (20,416 human + 16,343 Arabidopsis `ENTRAP_` proteins, same
`sage_deiso.json`), nominal `peptide_q<=0.01` in every row:

| arm | accepted | entrapment hits | raw entrap% | est. true FDR% |
|---|---|---|---|---|
| base | 9,502 | 91 | 0.96% | **2.39%** |
| corr040 | 9,811 | 89 | 0.91% | 2.27% |
| **corr060** | **10,342** | **86** | 0.83% | **2.08%** |
| the reference implementation | 10,729 | 81 | 0.75% | **1.89%** |

Poisson test — if the EXTRA peptides an arm accepts were as error-prone as base's, entrapment hits
would grow in proportion:

| arm | peptides vs base | entrapment observed | expected at base's rate | |
|---|---|---|---|---|
| corr040 | **+309** | 89 | 94 | −0.5σ |
| **corr060** | **+840** | **86** | **99** | **−1.3σ** |
| the reference implementation | +1,227 | 81 | 103 | −2.1σ |

**corr060 accepts 840 more peptides with NO measurable increase in errors** (86 observed against 99
expected). That is the signature of real signal, not recalibration — an artefact would have shown
entrapment hits growing with the accepted count.

**Why `fdrcheck.py` misled me.** It compared arms at a fixed `sage_discriminant_score`. That score
is a **per-run trained LDA output**, so its scale is NOT commensurable across runs — I flagged the
caveat when running it and then let the retraction headline stand anyway. Fixed-score comparison is
invalid here; direct error counting is not. Evidence ranking, weakest to strongest: nominal q-count
(per-run calibrated) < fixed-score (invalid across per-run-trained scores) < **entrapment (counts
errors directly)**. `bench/fdrcheck.py` is kept for the decoy/threshold diagnostics it prints, but
its cross-arm delta must not be read as evidence.

**The finding that outranks all of this: nominal 1% FDR is really ~1.9–2.4%.** Every "1% FDR"
number in this project's history — ours and the reference implementation's — is roughly **2× optimistic**. This is the
first time the project has measured its own error rate rather than modelling it.

**And the gap to the reference implementation is much smaller than reported**: 10,342 vs 10,729 under one
entrapment-controlled search = **−3.6%**, against the −12.0% (Sage) and −19.5% (MSFragger) nominal
figures. the reference implementation is also the cleanest arm (1.89%), so part of its lead is a stricter effective
threshold, not extra signal.

Caveat, stated because the counts are small: entrapment hits are 81–91, so Poisson noise is ±9–10
and the FDR *differences between arms* are NOT individually significant. The defensible claim is
"+840 peptides at no measurable error cost", not "corr060 has a significantly lower error rate".
Replication on dataset B/dataset A and under MSFragger is still required before the default moves.

## The four literature-motivated changes, benchmarked — three fail, and the failures are informative

Five arms, one binary, one invocation (node-1). Nominal `peptide_q<=0.01`:

| arm | change | peptides | vs base | verdict |
|---|---|---|---|---|
| base | — | 10,104 | — | — |
| wl025 | wavelet 0.25×FWHM | 10,069 | −0.35% | flat (within 0.34% noise) |
| wl050 | wavelet 0.50×FWHM | 10,101 | −0.03% | **flat — clean null** |
| rankint | intensity ranking | 9,855 | **−2.5%** | **falsified** |
| ms1sens | MS1 floor lowered | 5,605 | **−44.5%** | **falsified** |

**ms1sens — the physics, not a bug.** Parallelization worked (2.4 h vs never), and the diagnostics
explain the collapse: lowering the MS1 floor produced **66.6M MS1 traces vs 4.9M (13.6×)** and
**9.5M emitted spectra vs 1.5M (6.2×)**. The tool now emits a pseudo-spectrum at every noise blip;
those flood the target-decoy competition and bury real IDs. The literature premise — "a
pseudo-spectrum is created only if there is MS1 signal" — cuts both ways: relaxing it does not
surface faint real precursors, it manufactures spurious ones. **The MS1 floor is load-bearing at
its current setting.** The MS1 sensitivity direction is falsified; the parallelization that made
the test possible is kept (it is correct and free at default thresholds).

**rankint — falsified, against the literature's own logic.** Ranking the 500-fragment cap by
intensity (as the reference implementation does) LOST 2.5%. The reasoning "the engine re-ranks by intensity, so keep
intense peaks" was incomplete: our correlation-ranked cut is already a *quality* filter, and
replacing it with raw intensity lets co-eluting-neighbour peaks into the kept set. The correlation
ranking is doing real de-noising work that the reference implementation gets from its smoothed-XIC feature detection
instead. Closed.

**wavelet — a clean flat null.** Smoothing engaged correctly (`0.5 × 3.61 s = 1.81 s`, confirmed in
the log) and moved nothing at nominal FDR. Two honest readings, not yet separated: (a) the
grid-Pearson correlation is already robust enough that pre-smoothing the precursor is redundant, or
(b) the benefit is real but only visible under entrapment / on the faint tail, which the nominal
count cannot see. Given the gate sweep already showed a real entrapment-confirmed gain from
*tightening the correlation constraint*, (a) is the more likely explanation — smoothing and a
higher threshold are two routes to the same "trust the correlation more" effect, and the gate got
there first. Not adopted; not worth an entrapment run on a flat count.

**Net:** of the four post-review changes, one is a pure infrastructure win (MS1 banding, 6.4×),
one is a measured null (wavelet), two are falsified (intensity ranking, MS1 sensitivity). The
genuine peptide gain this session remains the entrapment-confirmed `min_correlation` tightening.

## Cross-dataset validation of the min_correlation winner (dataset B + dataset A, 2026-07-25)

The dataset D gate result was one sample; policy requires a second before the default moves. Tested
base/0.5/0.6 on dataset B and dataset A, one binary, one invocation each, then entrapment-scored.

**Nominal peptides (`peptide_q<=0.01`, human-only `sage_deiso`) — replicates on all three, and at
corr060 reaches near-parity with the reference implementation:**

| sample | base | corr050 | corr060 | Δ @0.60 | the reference implementation | corr060 vs dt |
|---|---|---|---|---|---|---|
| dataset D | ~10,120 (4 reps) | 10,657 | 10,903 | **+7.7%** | 11,517 | −5.3% |
| dataset B | 8,406 | 8,833 | 8,938 | **+6.3%** | 8,948 | **−0.1%** |
| dataset A | 8,983 | 9,870 | 10,059 | **+12.0%** | 10,242 | −1.8% |

**Entrapment true FDR — the win is unambiguous: peptides UP and error DOWN, on all three:**

| sample | base (acc / FDR) | corr060 (acc / FDR) | the reference implementation (acc / FDR) |
|---|---|---|---|
| dataset D | 9,502 / 2.39% | 10,342 / **2.08%** | 10,729 / 1.89% |
| dataset B | 8,016 / 2.40% | 8,459 / **1.95%** | 8,652 / 2.05% |
| dataset A | 8,452 / 1.71% | 9,622 / **1.53%** | 10,037 / 2.09% |

Tightening the gate does not just add error-free peptides — it **lowers the measured error rate
while adding peptides** on every sample. That is the opposite of the recalibration signature that
retracted the first reading. The lever is confirmed real and portable.

**Gap to the reference implementation, entrapment-matched (corr060):** dataset D −3.6%, dataset B −2.2%, dataset A −4.1% (mean ~−3.3%),
against the −12% (nominal Sage) and −19.5% (nominal MSFragger) figures. On dataset B and dataset A our corr060
error rate is *below* the reference implementation's, so at equal FDR the gap narrows further. The true standing vs
the reference implementation is low-single-digit percent, not double digit.

### MSFragger cross-engine (2026-07-25): the peptide GAIN is Sage-specific, the FDR gain is not

MSFragger now runs on our mzML (two blockers fixed: empty `labile_fragment_ion_series`, and
valueless cvParams → `Scans=0`, patched by `bench/fix_mzml_cvparams.py`). Same entrapment DB,
target-decoy FDR on hyperscore applied identically to both engines:

| engine | base (acc / FDR) | corr060 (acc / FDR) | Δ peptides |
|---|---|---|---|
| Sage | 9,502 / 2.39% | 10,342 / 2.08% | **+8.8%** |
| MSFragger | 10,032 / 1.67% | 10,026 / 1.42% | **−0.1% (flat)** |

**The +7–14% peptide gain does NOT replicate under MSFragger.** The FDR improvement DOES (both
engines get cleaner: Sage 2.39→2.08%, MSFragger 1.67→1.42%). Coherent reading: Sage's base
UNDER-performs MSFragger's base (9,502 vs 10,032 under entrapment); tightening removes
low-correlation fragments that confuse Sage's linear discriminant, closing that gap — but
MSFragger's hyperscore is already robust to those fragments, so it has nothing to gain. The
tightening compensates for a Sage-specific scoring weakness, it is not a universal extraction
improvement. This is exactly the engine-dependence the search-engine policy warns about, and it is
why Sage-only would have shipped a misleading +7% default.

**Decision: do NOT change the shipped default on this evidence.** Tightening is defensible as a
Sage-development convenience (neutral-to-positive on MSFragger via lower FDR, positive on Sage
counts), but the "+7% peptides" headline is retracted as a universal claim — it is a Sage artefact
the entrapment control on Sage alone could not expose, because entrapment measures error, not
cross-engine reproducibility.

Also gated: 0.5 vs 0.6 is within Poisson noise of the entrapment counts, so "tighten to ~0.5-0.6"
was always a range, not a value.

## Implemented after the review (2026-07-24)

Four changes, each traceable to a review finding or to the the reference implementation paper (Nat Commun 16:95),
none to parameter search — which matters, because parameter search on this file already produced
one retracted result.

**1. Band-parallel MS1 tracing (`perf:ms1_trace_bands`, default 12) — MEASURED, 6.4×.**
The cause was not "MassTraceDetection is serial" but simpler: **the MS1 call never passed a band
count**, so it defaulted to 1 while MS2 was banded across the machine. On the exact config that
previously ran 40 min with zero log progress:

| | unbanded | banded |
|---|---|---|
| peak CPU in MS1 phase | **152%** | **967%** |
| MS1 trace phase | no progress in 36 min | completes |

This unblocks the MS1 sensitivity floor, which the literature names as the root cause of
low-abundance loss for this entire tool class — DIA-Umpire lineage: *"a pseudospectrum is created
only if a peptide generates detectable signal in the MS1 data, which can limit sensitivity for
low-abundance species."* That is our miss profile exactly (MISS precursors 3.2× dimmer).

**2. Wavelet smoothing of the precursor XIC (`trace:wavelet_smooth`, default 0 = off).**
Stationary (à trous) B3-spline, scale = multiple of the **measured** MS1 FWHM. The literature
cross-check found we had copied the reference implementation's **0.3 correlation threshold without copying its
smoothing** (2D Gaussian + Savitzky-Golay on the precursor profile). Pearson r between two noisy
profiles is attenuated toward 0 as S/N falls, so an unsmoothed 0.3 is a systematically weaker
constraint — and worst at faint precursors. Raising the threshold to compensate was tried and
produced only the FDR artefact above; this attacks the cause. À trous (undecimated) is required,
not preferred: output stays sample-aligned, so the touched grid points — and therefore what
`gate:min_correlation_points` means — are unchanged. Fragment XICs stay RAW, as in the reference implementation.
Verified by `-diag:selftest_wavelet` against the SHIPPED functions (DC preservation, noise
reduction, peak retention, mirror edges, level-vs-sampling): 7/7 pass.

**3. Intensity ranking of the fragment cap (`assembly:rank_by`, default `correlation`).**
the reference implementation keeps the "top N highest intensity peaks"; we keep the top N best-correlating. The
engine re-ranks **by intensity** downstream (Sage `max_peaks=150`), so a correlation-ranked cap can
discard exactly the peaks that would have been scored. 70–75% of our spectra sit at the cap.

**4. Entrapment FDR (`bench/entrapment.py`).**
The measurement this project has never made. Counts errors directly against a foreign proteome
(20,416 human + 16,343 Arabidopsis) instead of trusting a per-run decoy calibration — the axis
that produced the retraction above. Reports raw entrapment fraction AND size-corrected FDR, never
the corrected figure alone.

## Harness fix landed in passing

`resolved_defaults()` parsed defaults from `--helphelp` "(default: '...')" text, which OpenMS 3.6
builds no longer emit — the parser silently returned `{}` and the FAILURE-5 no-op guard aborted every
run. Repointed at `-write_ini` (authoritative across versions; INI item paths minus the tool/instance
wrapper NODEs are exactly the CLI parameter names). Verified: 58 defaults parsed, correct values.
