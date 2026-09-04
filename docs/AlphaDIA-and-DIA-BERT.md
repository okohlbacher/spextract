# AlphaDIA & DIA-BERT — architecture review for spextract

Reviewed 2026-07-21. Both proposed as candidate model architectures for cycle
decomposition / annotation. Verdict: **steal AlphaDIA's stage-1 signal
extraction; do not build a transformer yet.**

---

## AlphaDIA (Wallmann et al., Nat Biotech 2025, 10.1038/s41587-025-02791-w)

Apache-2.0, Python + Numba/NumPy + PyTorch. `github.com/MannLabs/alphadia`.

### It cannot do what we do — our niche is unoccupied

Grep across the whole repo for `open search`, `mass tolerant`, `delta_mass`,
`blind search`, `de novo`, `denovo`: **zero hits**. No delta-mass dimension in
scoring or FDR. No pseudo-DDA/MGF export (`pseudo` matches only
`pseudo_reverse`, a decoy strategy).

Library-free only in the DIA-NN sense: in-silico FASTA digest + full PeptDeep
prediction, with the search space **closed and fully enumerated before search
begins** (`library_prediction:` — enzyme, `max_var_mod_num: 2`,
`missed_cleavages: 1`, `precursor_charge: [2,4]`, `precursor_mz: [400,1200]`).
The closed library is load-bearing for both candidate generation *and*
target-decoy FDR.

Its answer to unexpected PTMs is "predict the PTM into the library, then
transfer-learn onto it" — you must name the modification up front. It cannot
find a variant you did not enumerate.

**=> the reference implementation's pseudo-MS/MS-for-open-search niche is genuinely unoccupied by
the state of the art. This project is not redundant.**

### Stage 1 — "feature-free" — is what we should steal

Replaces discrete trace detection entirely. Peptide-centric, per candidate
(`search/selection/selection.py::_select_candidates_pjit`):

1. take top-k library fragments (default 12) + precursor isotope ladder
2. RT window (`get_frame_indices_tolerance`) + IM window (`get_scan_indices_tolerance`)
3. **cut a dense tensor straight from raw TOF-indexed data** (`get_dense_intensity`),
   gated by the quadrupole cycle mask. No centroiding, no trace linking, no threshold.
4. FFT-convolve with a 2D Gaussian (`fft.convolve_fourier`)
5. collapse: `log_fragment + log_precursor`, each `np.sum(np.log(smoothed + 1), axis=0)`
   over fragments — **log-sum, not pairwise correlation**
6. peak-pick local maxima (`find_peaks_2d`), keep top 3-5, merge within 3 scans / 3 cycles
7. integration bounds via `symetric_limits_2d`

Co-elution is **enforced by construction** (all fragments summed over the same
RT×IM grid) rather than measured post-hoc by correlation. Correlation features
appear later, as scoring features on already-bounded candidates.

Why this matters for us:
- **No detection threshold.** Weak fragments our `chrom_peak_snr × noise` apex
  gate silently drops still contribute in log space. That gate cost 401 peptides
  when merely misconfigured — this removes the bug class.
- **Log-sum > our zero-padded Pearson.** Likelihood under independence: one
  missing fragment doesn't zero the score (unlike a product); one huge fragment
  can't dominate (unlike a sum). Ours has neither property. This is review
  item #3, solved by a published benchmarked method.

### Paper oversells two things — check code, not text

- *"learned convolution kernels"* → shipping code is a **single plain Gaussian**
  (`selection/kernel.py::GaussianKernel.gaussian_kernel_2d`) with `n_features = 1`
  in `_build_features`. The "learning" is EM-style re-estimation of FWHM from
  confident IDs each round (`recalibration_handler.py:84-85` sets `fwhm_rt` /
  `fwhm_mobility` to the median `cycle_fwhm` / `mobility_fwhm`). Cheap, copyable,
  no ML infrastructure.
- fragment competition: paper says `k_max = 1` shared fragment; code says
  `if fragment_overlap >= 3` (`fragcomp/fragcomp.py:141`).

### Tensors

Selection (`get_dense_intensity`), 4D — intensity only:
`(1, n_fragment_mz_slices, n_scans, n_precursor_cycles)` float32

Scoring (`get_dense`, `absolute_masses=True`), 5D:
`(2, n_mz_slices, n_observations, n_scans, n_precursor_cycles)` float32
- dim0: `[0]` intensity, `[1]` observed m/z (or ppm deviation)
- dim2 `n_observations`: distinct quadrupole observations within one DIA cycle
  — what makes synchro-PASEF / midia-PASEF work
- dim3 scans = ion mobility, dim4 = RT in DIA cycles

Discretization happens at exactly one place: `_find_peaks` + `symetric_limits_2d`.
Everything upstream is continuous-grid.

(The 5D comment at `selection.py:163` is stale — the call is `get_dense_intensity`, 4D.)

### Ion mobility is first-class

Dedicated Bruker TIMS reader (`bruker_jit.py`) over TOF-indexed push data.
IM used at three levels: window bounding, 2D Gaussian + 2D peak-picking, and
5 of 46 features (`base_width_mobility`, `mobility_observed`,
`fragment_scan_correlation`, `template_scan_correlation`, `mobility_fwhm`).
Gated on `has_mobility` so Orbitrap drops them.

Cost, stated in the paper: explicit IM modelling gives **>1 h per file** for
large libraries, "will need improvement in future versions."

### Scoring: hand-crafted features + tiny MLP

`NUM_FEATURES = 46` (`constants/settings.py:5`), enumerated in
`search/scoring/scoring.py::DEFAULT_FEATURE_COLUMNS`. Groups: shape/location 5,
MS1 10, isotope 2, fragment agreement 8, ion series 4, profile correlations 9,
peak shape 3, mass error 2, **interference 3** (`n_overlapping`,
`mean_overlapping_intensity`, `mean_overlapping_mass_error` — we compute none of these).

Classifier: MLP hidden `[100, 50, 20, 5]`, **10,810 parameters**
(`fdr/classifiers.py:210`). `BatchNorm1d → (Linear→ReLU→Dropout)×4 → Linear → Softmax`.
Per-run SGD, 10 epochs, batch 5000, lr 0.001. Output P(decoy), count-based
target-decoy competition. Architecturally a Percolator-style rescorer.

### Chimeras: suppression, not deconvolution

`fragcomp/fragcomp.py` — candidates at 5% FDR bucketed by quadrupole window;
pairs within `rt_tol_seconds` compared at `mass_tol_ppm`; if shared fragments
>= 3, the worse-scoring one is invalidated (`valid_window[j] = False`).
**Winner-take-all at the identification level; shared intensity is never
apportioned.** Co-eluting precursors sharing few fragments are both reported.

=> **Neither AlphaDIA nor DIA-BERT deconvolves. That space is open.**

### Transfer learning

Fine-tunes four PeptDeep models (MS2 transformer, RT, charge, CCS) on **its own
first-pass output**: 1% precursor + 1% protein FDR, 3 best observations per
modified precursor, high-quality subset = median fragment correlation > 0.5.
Self-training loop. GPU **optional** — `get_torch_device` defaults to `"cpu"`,
upgrades only if `use_gpu: true` and backend available.

### Performance vs DIA-NN

Head-to-head numbers live only in figure panels (Fig 4b-c, 5b-d); no Source
Data file, so DIA-NN's counts are not extractable. Quotable AlphaDIA-only
numbers: HeLa timsTOF Ultra 21 min >73,000 precursors / ~6,800 protein groups,
median CV 7.7%. Lou mouse-in-yeast timsTOF: up to 81,500 mouse peptides,
"matching and even exceeding" DIA-NN/Spectronaut/MaxDIA. Dimethylated HeLa with
transfer learning: 96,000 precursors vs 65,000 without = **+48% precursors,
+25% protein groups**; RT median error 317 s → 11 s.

**Their strongest defensible claim is FDR calibration, not raw ID count** —
Arabidopsis entrapment holds 1% protein FDR while "some of the other tested
tools reported up to threefold more false-positive identifications than intended."

### Portability

Numba `@njit` over NumPy = effectively C. Directly portable (~1-2k lines):
`selection.py`, `kernel.py`, `fft.py` (swap pocketfft for FFTW/Eigen),
`utils.py` (`find_peaks_1d/2d`, `symetric_limits_2d`), `fragcomp.py` (~300 lines).
The 46 feature functions are pure NumPy reductions — mechanical.

**Real cost is `bruker_jit.py` (645 lines)** — `_assemble_push`,
`_get_push_indices`, `_cycle_mask` depend on alphatims' TOF-indexed sparse
layout (`tof_indptr`, `push_indices`, cycle model). We'd need the equivalent
index over our own representation. That's where the work goes.

Apache-2.0 is compatible with OpenMS BSD-3 for reimplementation (attribution + NOTICE).

---

## DIA-BERT (Nat Commun 2025, 10.1038/s41467-025-58866-4)

`github.com/guomics-lab/DIA-BERT`. Encoder-only transformer, 8 self-attention
blocks, 2 conv/avgpool/layernorm layers before it, BCE head.

**Wrong task, wrong data type.**

- **Library-dependent scoring, not decomposition.** Input is a **330 × 16 peak
  group matrix**: 330 rows = precursor m/z + 5 isotopes + up to ~325 *library*
  fragment ions; 16 columns = RT points (ref, -7, +8). The library fragments
  *are the query* — exactly the information we lack and are trying to infer.
  Cannot run without a library ⇒ cannot do blind search.
- **No ion mobility.** Paper explicit: timsTOF CCS values "are not included in
  the current version." v1.0 Orbitrap-only.
- Chimeric spectra not addressed.
- **GPU required.**
- Trained on **276M** precursors from 952 files. Positives = DIA-NN 1.8.1 ∩
  Spectronaut 14.6 at q<0.01; negatives = DIA-NN 1.7.12 decoys. Quant model:
  34M from 360 simulated files (6.87 TB, modified Synthedia).
- Reports +51% proteins / +22% precursors vs DIA-NN.

### The one valuable lesson: keep the RT axis, don't sum it

The 330×16 input keeps **16 RT points as separate columns** and lets conv +
attention decide how to combine them.

**Hypothesis (untested):** our N=5 cross-frame aggregation failed (fragments
211→412 but peptides 8,925→6,628) not because cross-frame information is
useless, but because we aggregated by **summation**, collapsing the per-cycle
elution profile into a scalar and destroying the shape that carries the
co-elution evidence. Test: retain per-frame intensities as an N-vector and
score on the profile rather than the sum. Contained change to `aggregateFrames_`.
This re-opens a direction previously recorded as falsified.

Second lesson: their +22% over DIA-NN is a *learned* scorer beating DIA-NN's
handcrafted co-elution features, which are considerably more sophisticated than
our zero-padded Pearson. External evidence that review item #3 leaves real
signal unclaimed.

Label-protocol note: they define positives as the **intersection of two
independent tools**. We can only partly copy this — Sage-on-our-pseudospectra
is not independent of our extraction, so it cannot serve as the second vote.
A limitation of our label set, not a protocol we have.

---

## Convergent verdict on "should we build a transformer?"

|  | DIA-BERT | AlphaDIA |
|---|---|---|
| scoring model | 8-block transformer | 46 hand-crafted features → MLP [100,50,20,5] |
| parameters | undisclosed | **10,810** |
| labels needed | **276M** precursors | per-run self-training |
| GPU | **required** | optional (CPU fallback) |
| ion mobility | **none** | first-class |
| our data type | no | yes |

The state of the art *in our exact data type* deliberately does **not** use a
transformer for scoring. It uses a ~10k-parameter MLP over hand-crafted
features and matches/exceeds DIA-NN, with FDR calibration as its strongest claim.

We have ~130k labels and **no GPU** (node-1: 128 CPU cores, no NVIDIA
device) against DIA-BERT's 276M + GPU requirement.

**Order of work:**
1. AlphaDIA-style stage 1: dense tensor → Gaussian → log-sum → peak-pick.
   Threshold-free, no library needed for the extraction half.
2. Compute AlphaDIA's feature set (esp. the 3 interference features we lack).
3. Fit a small MLP on the DIA-NN truth set.
4. Transformer only if (3) demonstrably saturates.

### Adapting stage 1 to our untargeted setting

AlphaDIA's stage 1 is **precursor-hypothesis-driven** — it convolves around a
known library m/z/RT/IM. We don't have that. But the halves differ:

- **Convolve + peak-pick needs no library.** Cut the full dense tensor per
  isolation window (all fragment m/z bins × IM × RT), smooth, peak-pick in 3D.
  Threshold-free fragment features with precise (RT, IM) apexes — strictly
  better input than our gated traces. Untargeted.
- **Log-sum needs a hypothesis, not a library.** We already generate precursor
  hypotheses from MS1 traces. Score each by the log-sum of smoothed fragment
  intensity at its (RT, IM) apex across candidate fragments. Drop-in replacement
  for the Pearson scorer.

This is a design problem for the grouping half, not a transcription job.

---

## Cross-tool file-format traps (found while building the truth set)

Both would have silently faked a comparison:

- **the reference implementation writes the precursor as its isolation window** (±0.01 Da offsets);
  we write the real acquisition window (26 variable-width windows, 19–141 Da).
  Trusting either file's window empties the reference implementation's co-isolation set by
  construction and makes it look perfectly pure.
- **RT units differ**: our mzML is in seconds, the reference implementation's in minutes, DIA-NN
  reports minutes. Read the CV `unit_info`, never the bare number.
- Precursor m/z key differs: we write `selected ion m/z`; the reference implementation omits it and
  only has `isolation window target m/z`. IM is in `selectedIon` for us, in
  `scan` for the reference implementation.
