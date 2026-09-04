# Next project: close the recall gap upstream (feature detection), not at the fragment level

## Why this project exists
The recall gap vs the reference implementation on dataset D is ~1,500 peptides (9,989 vs 11,517, ~87%), entirely the
**low-abundance tail** (MISS precursors are ~3.2× dimmer; nothing else discriminates). Seven
fragment-level levers are now falsified (deisotope, adjacent completion, sharing-subtraction,
correlation-competitive, apex-competitive, logoverlap gate-statistic, greedy peel-bright
deconvolution) — all fail because fragments are shared across 11–50 co-eluting precursors, so any
redistribute/reduce/remove/re-score depletes the faint precursors. **The lever is not in fragment
assignment.** It is upstream, in how the precursor (and its XIC) is detected and defined.

## The the reference implementation difference (from the paper, Nat Commun 16:95 / PMC11160675)
the reference implementation's FRAGMENT side is identical to ours (Delta Apex IM 0.02, RT 3, Corr 0.3, RFmax 500,
top-N by intensity — we copied these). Its advantage is entirely UPSTREAM:
1. **Neighbouring-scan aggregation BEFORE centroiding** — "data in the same isolation window of
   adjacent neighbour RT frames are aggregated to form a composite 3D matrix by summing intensity...
   extracts weak signals and filters out noise." Enhances faint SIGNAL (not a lower threshold).
2. **2D (m/z–1/K0) adaptive Gaussian feature detection** — 2D Gaussian smoothing of the m/z×IM
   matrix, local maxima as seeds, **Gaussian fitting** to determine peaks; Savitzky-Golay + Z-score
   for elongated peaks. Produces smooth, accurate precursor XICs BY CONSTRUCTION.
3. Correlates the **smoothed monoisotopic precursor XIC** vs raw fragment XICs at the 0.3 gate.

Our own doc already isolated the mechanism: *"the correlation ranking is doing real de-noising work
that the reference implementation gets from its smoothed-XIC feature detection instead."* We use OpenMS greedy
`MassTraceDetection` (m/z tracing, valley-split) → raw, noisy XICs → at low S/N the fragment
correlation is a lottery → true faint fragments gated out / interferents admitted → faint peptide's
spectrum is unidentifiable at the engine's intensity top-150 (measured: our top-150 b/y coverage
0.083 vs the reference implementation 0.208).

## What we ALREADY tried (do not repeat)
- **Wavelet post-smoothing of our precursor XIC:** benchmarked FLAT (−0.03%, engaged correctly).
  Post-smoothing a greedy XIC ≠ a Gaussian-fit XIC. Not the fix.
- **Lower MS1 noise floor:** −44.5% (13.6× more traces, spurious). Threshold relaxation is wrong.
- **Intensity fragment ranking:** −2.5% (correlation ranking is load-bearing de-noising).
- **MS2 frame aggregation AFTER centroiding:** −26%, but UNRESOLVED — mechanistic or a peptide_q
  redundancy artefact ("under peptide_q, redundancy is free score").
- All 7 fragment-level levers (above).

## What we MISSED / have NOT done
1. **Native pre-centroiding neighbour-frame aggregation** (the reference implementation's actual weak-signal method).
   Implemented as `trace:native_ms1_neighbors` / `native_ms2_neighbors` (+ `perf:stream_load`, sums
   RAW frames in the reader before centroiding) but **NEVER BENCHMARKED**. We only tested the
   after-centroiding version (−26%). Pre- vs post-centroiding is a different operation.
2. **2D (m/z–1/K0) adaptive Gaussian feature detection** — the reference implementation's core upstream method, the
   source of its smoothed-XIC de-noising. We use greedy `MassTraceDetection`. NOT implemented.
3. **Isotope-distribution charge determination + precursor-cluster de-isotoping** (the reference implementation removes
   confident isotopes from further clustering). We default to partner-count charge; envelope mode
   exists but is off. Partially explored.

## Proposed plan (hypothesis-ranked, entrapment-gated)
Every arm is judged on ENTRAPMENT FDR (not nominal peptide_q), because the "redundancy is free
score" confound makes nominal counts unreliable for anything that changes spectrum/precursor count.
Both Sage AND MSFragger (standing policy).

- **Tier 1 (days): native pre-centroiding aggregation.** Benchmark `native_ms1_neighbors` and
  `native_ms2_neighbors` (1 = 3-frame sum, 2 = 5-frame sum) with `stream_load`, on the gated
  baseline, under entrapment FDR. This is the single cheapest test of the reference implementation's actual weak-signal
  mechanism and it is already coded. Resolve the 26% question at the same time (native vs
  after-centroiding, entrapment-gated). Decision gate: if faint-tail recall rises at neutral true
  FDR, this is a large part of the gap for near-zero engineering.
- **Tier 2 (weeks): 2D m/z–1/K0 adaptive Gaussian feature detection** as an alternative MS1 (and
  optionally MS2) feature detector to greedy `MassTraceDetection`: neighbour-aggregate → 2D Gaussian
  smooth the m/z×IM matrix → local-maxima seeds → 2D Gaussian fit → RT tracing with SG/Z-score
  segmentation → Gaussian-fit (smooth) XIC used for fragment correlation. This replicates the
  smoothed-XIC de-noising directly. Validate the XIC-quality proxy first (does a Gaussian-fit XIC
  raise the top-150 fragment coverage toward 0.208 on the MISS set?) before full integration.
- **Tier 3 (optional): precursor-cluster de-isotoping + isotope-distribution charge** to sharpen
  precursor definition, only if Tiers 1–2 leave a residual.

## Risks / falsifiers to pre-register
- Native aggregation may just manufacture more spectra (the −44.5% MS1-floor failure mode) — gate on
  entrapment FDR and on faint-tail recall specifically, not total count.
- The 2D Gaussian detector is a large surface; validate the XIC-quality proxy on the MISS set BEFORE
  committing to full integration.
- Aggregation autocorrelates adjacent RT points → inflates Pearson r → must re-tune the correlation
  gate and re-verify under entrapment.

---

## ADVERSARIAL REVIEW VERDICT (codex + vibe + kimi, 2026-07-26): scope REVISED

**All three reviewers converge: do NOT build the Tier-2 2D-Gaussian detector yet — it is likely the
8th falsification.** Its stated rationale (a smoother correlation-reference XIC) was already tested by
the wavelet arm and returned flat; "post-smoothing ≠ Gaussian fit" is a shape quibble, not an
information argument. Three critical redirections:

1. **The gap may be SEARCH-SIDE, not extraction (all three).** the reference implementation's validated pipeline is
   FragPipe (MSFragger + MSBooster + Percolator rescoring); we benchmarked on Sage with
   `chimera:false, report_psms:1`. A faint peptide buried under a bright co-isolate loses the
   single-PSM competition and gets no MS1-informed rescoring — a SEARCH failure. **Free pre-step
   confirms it:** of the 3,032 MISS peptides, 1,025 (34%) already have a near-threshold PSM (q
   0.01–0.10) in our output and ~70% have SOME PSM; only 918 (30%) have no PSM at all. So a large
   fraction is rescuable at search time, not extraction.
2. **Charge/monoisotope is a tail-specific factor (kimi)** — funnel attributes 26% of precursor loss
   to metadata; the `count`-vs-`envelope` charge A/B history should be verified (current default IS
   count = 74.6% agreement vs the old envelope 47.8%; the reference implementation 82.6%).
3. **IM as a continuous per-fragment weight (codex + kimi)**, `exp(-dIM^2/2sigma^2)`, down-weights
   interferents WITHIN a spectrum without redistributing fragments ACROSS precursors — escapes all 7
   falsifications (which were cross-precursor) and directly un-buries faint fragments under the
   engine's top-N-by-intensity. Hours to test. NOT in the falsified list.

**The decisive cheap experiments (do these before any detector build):**
- (a) Re-search existing spectra with Sage `chimera:true, report_psms:5`, and via
  MSFragger+MSBooster+Percolator — recovery of the fixed MISS set (RUNNING).
- (b) IM continuous-weight fragment scoring — coverage proxy on the MISS set, then search.
- (c) Oracle decomposition (kimi Finding 7): for each MISS peptide emit the TRUE fragments at true
  charge/mono, search — factors the gap into metadata vs content vs search, and tests whether the
  ceiling is information-theoretic at all. Also fix the entrapment estimator (codex #12: use
  peptide/precursor hypothesis ratio, not protein ratio).

**Tier 2 (2D Gaussian) is deferred** pending a demonstrated CONTENT gap after (a)-(c). Tier 1 (native
aggregation) remains a cheap valid test for the 13%/30% not-emitted, but is not the main lever.
