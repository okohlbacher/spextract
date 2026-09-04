# Project review — 2026-09-01 (8-agent fan-out: 4 mappers, 3 judges, 1 red-team)

Full rendered report: claude.ai/code/artifact/258035d0-0fb0-4782-be7f-12d4eca6c7a7
State reviewed: branch coverage-analysis-2026-07-24 @ d4686f0.

## Verdict
The science phase succeeded and is over-documented; the product phase has not started and is
under-shipped. Highest-leverage week = turning finished, honest science into something a third
party can build, run, and verify — not another quality lever.

## Scoreboard (dataset D; engines deterministic sigma=0; the reference implementation denominators NOT yet reconciled)
| metric | SpeXtract | the reference implementation | ratio |
|---|---|---|---|
| closed Sage (corr_power=2 default) | 10,802 | 11,517 | 93.8% (+im_weight opt: 11,112 = 96.5%) |
| closed MSFragger | 11,560 | 13,932 | 83.0% |
| OPEN search targets @ ~1% entrap FDR | 8,887 @1.29% [1.02-1.58] | 9,989 @1.07% [0.85-1.33] | 89% — UNEXPLAINED |
| isotope-shifted PSMs (phantom-mod class) | 8.20% | 3.46% | 2.4x — PTM-interpretation defect, NOT recall |
| emission | 924,255 | 700,434 | 1.32x |
| extraction wall | ~19:30 | 7:47 | 2.5x (4.1x compute: ~1x intrinsic, ~1.5x waste R2-R4, ~0.5x constants) |
| FDR calibration @ nominal 1% | 1.0-1.45% (proven, CIs) | ~1.1% | the one axis we are AHEAD (provable) |

## Red-team headline
**The flagship open-search metric has never been measured in the right currency, for either
tool.** open_bench counts target peptides@FDR, but the raison d'etre is PTM-SITE discovery;
the -1.003 Da work proved the currencies diverge (a shifted precursor still counts as a target
while being a wrong PTM assignment), and we carry 2.4x the reference implementation's rate of that error. Scored
delta-mass-bin-level (known-PTM bins vs artifact bins, entrapment-controlled per bin), the 89%
could move either direction. One script against existing output; must precede any new open-search
experiment. Second: two unreconciled the reference implementation reference sets (Sage 11,552/11,517; MSFragger
13,014/13,932 — 7% apart) + the 9,989-vs-10,072 baseline discrepancy → no ratio is citable yet.

## Key gaps (judges)
- Science: open-search gap unexplained (live candidates: emission competition, faint-tail
  quality; ppm is interpretation-class — audit free, don't build). Generalisation is technical
  (same-series holdouts); patient data unreproducible in principle → public-data replication
  required. Point comparisons at unequal actual FDR → need full FDR curves at matched FDR.
- Engineering: zero tests/CI (fatal for a determinism-headline tool); README 27 commits stale;
  main 27 behind; v0.1.0 tag burned (day-one commit); EPD determinism patch documented nowhere;
  evidence chain on purgeable /scratch (verified INTACT 2026-09-01 after the 34-day gap — luck,
  not policy); nothing verifies repo-src == cluster-binary; entrap_apply.py imports untracked
  /scratch modules. Code landmines: apportion/rp_max silently bypass corr_power/im_weight
  (3 hand-synced emit paths); falsified flags' help still RECOMMENDS them (rp_max "useful
  range", consolidate "Try 5-10", merge "principled window"; merge help says "summed" but code
  does MAX); charge:iso_im_tolerance default 0.05 contradicts its own help; macOS silently
  serialises the window loop; only the mono peak is excluded from fragment lists (M+1/M+2 can
  emit as fragments); mono can cross isolation-window routing boundary.
- Product: adoption case is structural (license-blocked users, pipeline builders, calibration
  trustworthiness), not benchmark parity. Never claim parity. Binary must be renamed SpeXtract.
  Closed-search gap chasing = non-goal (corrupt metric).

## THE PLAN (implementing in this order; kimi+codex adversarial review after each step)
00. TODAY — archive the evidence chain: untracked /scratch analysis modules + engine configs +
    MSFragger params into git; cmp outputs/logs into git + /ceph tarball; FASTA SHA256s;
    provenance stamping (src sha recorded at deploy/build; harness refuses on mismatch later).
01. Re-denominate the flagship metric (days, no new runs): delta-mass-bin-level entrapment
    scoring of BOTH tools' existing open output; reconcile the reference implementation denominators + 9,989/10,072.
02. ONE pre-registered emission-controlled open-search arm (1 week, time-boxed): quality-gate
    924k→~700k + the reference implementation-precursor-matched sub-arm, identical search, full FDR curves, scored
    under the 01 metric. Gap closes → emission competition → principled gate. Flat → faint-tail
    quality → roadmap. Ship regardless.
03. Public-PXD replication + v0.2.0 as ONE unit (1-2 wks): the reference implementation paper's own public
    diaPASEF data end-to-end (both tools, both engines, corrected estimator); then main FF,
    README rewrite (dataset D-BASELINE + public numbers), EPD patch + OpenMS pin docs, apportion/
    rp_max guard, falsified-flag help pass, binary rename to SpeXtract, build CI + determinism
    smoke test on public-derived fixture, tag v0.2.0 (v0.1.0 burned) + CHANGELOG + CITATION.cff.
04. R2+R3+R4 runtime batch (3-5 d, low risk): cull span-0 MS2 traces, two-lower_bound band
    construction (byte-identical), double-score dedup → −25-40% wall → headline <2x.
ROADMAP (in README): Q5 rescoring sidecar (FIRST recover the chimeric+MSBooster arm launched
07-26, result never recorded — nearly free); Q4 mono-offset stamp (mandatory once 01 exists);
R5 task graph (only behind CI determinism gate); upstream OpenMS PR (separate product, gated on
EPD patch PR + public fixture).
NON-GOALS: closed-search parity; merging/collapse; runtime parity pre-R5; im_weight default;
-1 Da fixes as defaults.

## History lessons (mappers)
Wins are within-spectrum/upstream, losses cross-precursor ("share-all at full intensity is
robustly optimal"); cheapest decisive arm beats best-reasoned design; the metric is the project;
adversarial review is a control, not an oracle — run the decisive arm anyway; determinism is a
force multiplier.
