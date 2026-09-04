# Self-calibration options for emitted spectra (deep-research + 3-CLI adversarial review, 2026-09-01)

Basis: deep-research fan-out (104 agents, all claims 3-0 verified; sources: COFI/JPR2005,
delta-histogram recal/AnalChem2025, DOLM/AnalChem2011, MaxQuant software-lock-mass, MSFragger
calibrate_mass, DIA-NN, OpenMS InternalCalibration, Waters patent US9991103B2) + adversarial review
by codex, kimi, vibe of options A-H.

## THE CONCISE LIST (integrated verdict)

1. **B+C as ONE joint 2-parameter fit — the v0.2.0 open-fallback fix.** (all 3 reviewers, kimi's
   framing: they are complementary constraints, not competitors.)
   - B = search-free ABSOLUTE anchors: MS2 immonium ions (70.065/86.096/110.071/120.081/136.076) +
     tryptic y1 K/R (147.1128/175.1190); MS1 ambient siloxanes/phthalates, presence-gated.
     Anchors pin the OFFSET (t0) exactly where its m^-1/2 error is largest (low m/z).
   - C = residue-mass-delta histograms (top ~15 amino-acid deltas + isotope spacings; NO H2O/NH3
     neutral losses -- patent-vocabulary exclusion) pin the SLOPE across the full range, covering
     B's anchor hole above m/z 175.
   - Fit the two sqrt-chord parameters (bounded model; fail closed on inadequate anchor support);
     per-RT-chunk (fixed 60-300 s); deterministic robust fitting (Theil-Sen / exhaustive subset --
     NOT RANSAC, even seeded: breaks byte-determinism across library versions); centroid in TOF
     space BEFORE m/z-dependent gates (no circularity); two-pass (collect -> freeze model -> emit);
     collapse to distinct traces (no pseudoreplication); decoy masses + held-out anchors as QC.
   - Consensus accuracy bet: median 8.3 -> **1.5-4 ppm central** (bow -> ~+-2 ppm), 4-8 ppm tails.
     NOT sub-ppm -- that territory belongs to downstream search-dependent tools (0.3 ppm ceiling).
   - Effort: ~1 week MVP (B) + C as its shape constraint; hardening 2-4 more weeks.
2. **Option I (was missing from A-H): centroiding/mono-reporting forensics — the ONLY fix for the
   vendor-path +1.7 ppm residual + mass-degraded tail** (unanimous: no axis calibration touches it).
   Step 1 (one afternoon): residual-vs-INTENSITY diagnostic on raw frames (saturation-pull
   hypothesis, kimi) + codex's stage-by-stage mass ledger (raw TOF -> frame centroid -> IM
   aggregation -> trace representative -> isotope/charge selection -> emitted m/z); fix the first
   divergent stage (candidates: TOF-space profile-fit centroids, saturation censoring,
   envelope-fit mono instead of apex). Realistic post-fix target: 0.5-1 ppm central.
3. **E (isotope-spacing) = per-file QC DIAGNOSTIC only, never a corrector** (5 uDa signal at 5 ppm;
   derivative-only -- cannot fix the gauge). Ship as a bow-shape QC readout.
4. **H (open port of the vendor polynomial) = quarantined cross-check during development only**
   (upstream author disabled it; bimodal failure risk). Never a release gate.

## KILLED / REJECTED
- **D (Waters neutral-loss regression, US9991103B2): legally dead** (in force to ~2035; claim 1 is
  broad -- not limited to regression, DIA expressly covered; BSD gives no patent immunity;
  inducement exposure). Also technically weaker than C (17/18 Da losses = poor ppm leverage).
  Architectural rule: NO precursor->fragment neutral-loss pairing anywhere in calibration code.
- **F (2D RT x m/z surface): premature** -- a model shape in search of a signal; only after B+C
  density is demonstrated. **A/G: release-gate violations, not options.**

## Resolved disagreements
- Gauge: vibe claimed C/D/E can fix the absolute gauge alone; kimi/codex's worked math wins --
  deltas see multiplicative error at full strength but the offset only at ~Delta/m (10x attenuated;
  a pure additive offset cancels exactly). Hence B+C joint, not C alone.
- C standalone accuracy: 2-8 ppm with false-lock risk (false pairs grow N^2 vs true N; mode
  survives only with charge/window/mobility coherence filters + local background model).
- COFI's 2.49 ppm is roughly the 40k-resolution CENTROID FLOOR, not an aspiration (kimi).

## Interactions with live work
- The deterministic-fit requirement dovetails with the det8 FAIL (threads-8 vs threads-100 mzML
  differ by 3 bytes at identical spectrum counts -- last-ulp thread-order sensitivity; no-SDK
  control running to attribute).
- Option 1 dissolves the two-tier product conflict (license-blocked users get ~2-4 ppm open-path
  masses); the vendor SDK remains the benchmark/best path.

## REASSESSMENT (user directive "polynomials are not dead", 2026-09-01) — the polynomial IS the answer

Empirical results that invalidate the ranking above for TDF inputs:
1. **Model-family check vs the vendor oracle:** the 2-param sqrt-chord (which the review said to
   refit) CANNOT represent vendor truth (max residual 13.0 ppm); a QUADRATIC in sqrt-space fits it
   to 0.0001 ppm, frame-stable to ~1e-7 across the run -> ONE per-run polynomial suffices.
2. **The closed form was derived and verified EXACT (3e-5 ppm max, all 3 files):**
   t_ns = C0 + (1e6/sqrt(C1_eff))*sqrt(m) + C2*m,  C1_eff = C1*(1+dC1*(T1_ref-T1_frame)/1e6),
   t_ns = tof*DigitizerTimebase + DigitizerDelay -- ALL constants from the TDF MzCalibration table.
   Clean-room numeric derivation against locally-licensed vendor output; no vendor code consulted.
   use_recalibrated_state made zero difference on these files; temperature term = the 0.6 ppm
   residual of the C2-less form.
3. **Upstream bug found:** the reference open implementation (timsrust-calibration; adapted by
   mzdata then disabled as "not consistently better") DROPS the C2*m term -> -11..-40 ppm vs vendor
   on this data. That omission, not translation subtlety, is why the port "wasn't better".
   Report upstream.
4. **Implemented:** TdfTableTof2MzConverter in the loader (patches/openms-brukertims-mz-calibration.patch,
   supersedes+includes the sdk-mz patch): m/z priority = vendor SDK (if env set) > TDF-table exact
   model > legacy chord (warn). License-free, deterministic, no fitting.

**Reassessed ranking:**
1. **H1 table-driven exact model (implemented)** -- vendor-exact open path for ALL TDF users incl.
   license-blocked; dissolves the two-tier gate entirely. 3-file comparison arms running.
2. **Option I (centroiding/mono forensics)** -- unchanged; still the only fix for the +1.7 ppm
   residual + mass-degraded tail (which are downstream of ANY axis calibration).
3. **B+C statistical fit** -- demoted to (a) QC/verification layer on top of H1, (b) fallback for
   non-TDF inputs or ModelType != 1 files. Its 1.5-4 ppm bet is now the fallback-of-the-fallback.
4. **E isotope-spacing QC** -- unchanged (diagnostic).
D remains legally dead; F remains premature; the vendor SDK remains as cross-check/benchmark oracle.
detctrl determinism control VOIDED (self-kill + mid-run rebuild) -- re-run on the settled binary.
