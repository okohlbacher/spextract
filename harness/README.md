# SpeXtractor benchmark harness v2

`bench2.py` produces artifacts and provenance. `collate2.py` produces numbers. Nothing else may
write into `docs/benchmarks/`.

```
bench2.py   plan-check --plan P              # validate a plan; no side effects, no dirs created
bench2.py   pin        --sample dataset A          # write pinned content digests into samples.yaml
bench2.py   run        --plan P              # execute; prints a runset id
bench2.py   verify     --runset R [--deep]   # rehash sealed artifacts; writes nothing
bench2.py   show       --runset R
collate2.py truth      --runset R            # compute metrics; writes runset revision r2
collate2.py collate    --runset R.r2 [--publish]
bench2.py   selftest                         # 64 guard cases
collate2.py selftest                         # 88 guard cases, incl. two end-to-end runs
```

Exit codes: `0` ok/published · `2` guard abort · `3` refusal to publish.

---

## Why v2 exists

Ad-hoc benchmarking produced eleven data-handling failures in one session. v1 added eleven
guards. Two independent adversarial reviews found that **six were theatre** — paired baselines,
scale checking, binary identity, reference existence, complete provenance, process exclusivity —
the claimed invariant had no implementing code. And the largest hole was one no guard addressed:

> **No immutable run identity.** A nonzero converter exit left a stale nonempty mzML that was
> accepted; Sage failed leaving an old `results.sage.tsv` that was rescored; the manifest was
> overwritten with fresh metadata. One vault row could be assembled from three different attempts.

The root cause of that whole class is **location addressing**. `outroot/{sample}__{arm}/pseudo.mzML`
is a *name* that means different things at different times, so a write that fails leaves the
previous meaning in place and the next reader cannot tell which write it is looking at.

The fix is not a check. **The tool writes into a directory that did not exist one second earlier,
whose name embeds a nonce no other process knows.** There is no pre-existing `pseudo.mzML` to
survive a nonzero exit, because there is no pre-existing directory. Everything else in v2 is
scaffolding around that one property, plus digest-chaining so it survives across stages.

Second architectural change: **`truth` and `collate` no longer glob.** v1's
`root.glob("*__*")` under a directory derived from the *plan filename* collated whatever happened
to be present — including a copied directory carrying a foreign manifest, and including nothing at
all (an empty root produced an apparently valid empty vault report). In v2 **no glob discovers
work**: the sealed runset contract enumerates every cell and set equality with the plan is exact.
(One glob remains, in `load_runset`, and it only resolves which *revision* of an explicitly named
runset id to open — it cannot find a runset you did not name.)

---

## Tiers

Used consistently below and in the source. The reviewers' central complaint was that v1's
docstrings claimed more than the artifacts delivered, so the tier is always stated.

| tier | meaning |
|---|---|
| **STRUCTURAL** | the wrong thing cannot be expressed; there is no code path that accepts it |
| **VERIFIED** | checked against a digest committed before the data existed |
| **CHECKED** | a runtime abort on an observable condition |
| **RECORDED** | measured and stored, **not** prevented |

---

## Guards — `bench2.py`

Every row aborts. Every row has at least one `--selftest` case that deliberately triggers it and
asserts *that specific guard id* fired.

| id | tier | prevents | cases |
|---|---|---|---|
| **G-B01** | STRUCTURAL | Unknown plan key aborts. The `basline:` typo class: it parses, `.get("baseline")` ignores it, and the arm publishes as if unpaired by design. Schema version mismatch aborts — no forward-compat guessing. | 3 |
| **G-B02** | STRUCTURAL | Arm/sample name regex + duplicate detection. **FAILURE: arm names were path components** — `../escape` left the tree, a duplicate silently overwrote another arm's directory. In v2 names are never path components at all. | 2 |
| **G-B03** | STRUCTURAL | Every treatment resolves to a baseline **in this plan**. **FAILURE 2**: v1 did `if not arm.baseline: continue`, skipping all validation for an arm that omitted the field — then ran and published it. Baseline chains deeper than 1 abort. | 3 |
| **G-B04** | STRUCTURAL | A parameter listed in registry `tool.scaled_params` may **only** enter via a ≥3-point sweep. **FAILURE 6**: `ambiguity_margin: 0.1` compared at a single point against an *integer* partner count. The registry declares it, not the plan, so a hand-written arm cannot opt out. | 3 |
| **G-B05** | VERIFIED | Sample must be pinned; `raw_content_id` and `acquisition_id` must match the pin. **FAILURE 1**: v1 checked that the string `"dataset A-A-3"` appeared in the path. A symlink retargeted between preflight and execution defeated that completely. | 3 |
| **G-B06** | VERIFIED | References, Sage binary, Sage config and FASTA must exist **and hash to their pinned digest**; a bare path with no digest is refused; the Sage config is parsed and its FASTA compared to the registry. **FAILURE 9**: v1 printed `MISSING` and continued. Also: `reg["search"]["fasta"]` was declared and never used by anything. | 5 |
| **G-B07** | CHECKED | `flock` host lock held for the whole run, with write-readback. Replaces v1's `ps` substring scan, which was race-prone, BSD-fragile, defeated by a renamed binary, and had a real bug (`str(pid) not in line` false-matches 123 against 5123). **Scope is honest — see UNCOVERED.** | 1 |
| **G-B08** | CHECKED | `staging/`, `stages/` and `failed/` must share a device. Put staging on node-local disk and `os.rename` raises `EXDEV`; the reflex fix is `shutil.move`, a non-atomic 12 GB copy, and the entire commit protocol collapses silently. | 1 |
| **G-B09** | VERIFIED | The binary is hashed immediately **before and after every exec**. **FAILURE 7**: v1 hashed once at preflight and stamped that cached value into every manifest, so two different binaries could execute while both manifests carried binary A's hash. v1's comment claimed identity was "asserted equal"; no assertion existed. | 1 |
| **G-B10** | STRUCTURAL | A nonzero return code aborts the attempt. No "but the file exists" branch, for the converter *or* Sage. v1's `if not mzml.exists() or size == 0` accepted a failed run whenever a previous attempt's output was present; Sage's return code was discarded outright. | 4 |
| **G-B11** | CHECKED | Declared outputs exist, are nonempty, and are **complete**: mzML ends with a closing tag, `<spectrumList count="N">` equals the observed spectrum count, TSV ends with a newline, every row has the header's field count, required columns present. | 11 |
| **G-B12** | STRUCTURAL | Fresh attempt dir → fsync → atomic `os.rename` commit. Death before the rename leaves work in `staging/`, which no reader enumerates; death during is atomic; death after cannot expose unflushed bytes. Seals detect post-commit edits. | 8 |
| **G-B13** | CHECKED | The raw `.d` is content-digested before the run and its `(dev, ino, mtime_ns, size)` plus `analysis.tdf` digest re-verified **after**. A retarget during the run aborts the attempt, so nothing seals and nothing publishes. **Not prevention — detection.** See UNCOVERED. | 1 |
| **G-B14** | STRUCTURAL | No-op detection by **output-digest equality** with the baseline, plus resolved-INI equality. **FAILURE 5**: three arms passed values that were already defaults and burned ~2.5 h. v1 regexed `--helphelp` and diffed strings, missing `0.1` vs `0.10`, boolean aliases, aliased flags, tool-side clamping, and parameters ignored under another mode. Output-digest equality is exact. | 5 |
| **G-B15** | STRUCTURAL | Stages are reachable only through the sealed runset contract — `load_stage` takes a stage id, and stage ids exist only inside a sealed runset. No code path discovers work by listing a directory. | 1 |
| **G-B16** | STRUCTURAL | Converter argv is restricted to `-ini/-in/-out`, and a parameter absent from the tool's own resolved INI aborts. A command-line parameter would override the INI and be invisible to `ini_sha256`, making no-op detection blind to it. | 3 |
| **G-B17** | STRUCTURAL | Exact set equality between the plan's cells and the runset's cells — no missing, no extra, **no duplicates**. Two completed runs for one `(sample, arm)` is a refusal, not a coin flip. | 4 |
| **G-B18** | CHECKED | Acquisition windows are extracted **per sample** from that sample's `DiaFrameMsMsWindows`; a missing table aborts. **`WINDOWS_S23` is deleted and there is no default table.** v1 applied dataset B's 24 one-dimensional tiles to dataset A and dataset D with nothing recording that they shared a method. | 4 |
| **G-B19** | CHECKED | Free-space preflight. On this cluster quota exhaustion silently kills jobs; a convert that dies at 95% is indistinguishable after the fact from a tool bug. | 2 |
| **G-B20** | STRUCTURAL | One canonical serializer with `allow_nan=False`, and `parse_constant` rejection on read. `json.dumps` emits bare `NaN` by default, which is not valid JSON and silently poisons every downstream reader. | 2 |

## Guards — `collate2.py`

| id | tier | prevents | cases |
|---|---|---|---|
| **G-C01** | CHECKED | RT unit read from a **closed** accession map; unknown or absent unit aborts; the unit must be constant within a file. **FAILURE 8**: RT in seconds read as minutes matched 31 of ~300,000 spectra. v1 read the CV unit — correct — but then `getattr(v, "unit_info", "minute")` *guessed minutes* exactly where it needed to refuse. | 5 |
| **G-C02** | STRUCTURAL | Precursor m/z resolved by two explicit branches, and the source (`selected_ion` vs `isolation_target`) is counted and folded into `estimand_id`. v1's `.get(a, .get(b, 0))` default fires only on **absence**, so a present-but-zero m/z skipped the fallback and the spectrum vanished. | 2 |
| **G-C03** | STRUCTURAL | Missing charge → `None`; missing ion mobility → `None` and a named disposition. **A defect neither review found:** v1's `im ... or -1` combined with a filter reading `im < 0 or abs(...) <= D_IM` **silently disabled the ion-mobility gate** for that spectrum. The effective tolerance changed mid-file. | 2 |
| **G-C04** | CHECKED | Sage schema contract: required columns, duplicate-column check, pinned header digest, per-row field count, typed coercion with domains, `NA`/blank/NaN rejection, and every row's `filename` must name this cell's mzML. **FAILURE 3**: `x.get("peptide_q", 1)` turned a renamed column into a published bold zero. | 7 |
| **G-C05** | STRUCTURAL | `peptide_q` is the published metric; `spectrum_q` lives in a diagnostics sub-object and **never appears in a vault table row** (asserted against the rendered report). v1 printed it in the headline row where it was quotable; it inflated counts by 9–11 points against the reference implementation. | 2 |
| **G-C06** | STRUCTURAL | One division in the file. Zero denominator → `ZERO_DENOMINATOR`; below `n_min` → `INSUFFICIENT_N`. Replaces `max(len(sp),1)`, `max(tot_z,1)`, `max(len(present),1)`, `max(1-recd,1e-9)`, `i2.sum() or 1.0`. `Value` rejects NaN and `n=0`; `Value / Value` is a `TypeError`. | 5 |
| **G-C07** | STRUCTURAL | Empty observation list → `EMPTY_POPULATION`. Replaces `np.array(own or [0])`, which manufactured a synthetic zero observation. | 1 |
| **G-C08** | CHECKED | A missing or thin fragment library **aborts**. v1's `if lib and lib.exists()` left every fragment set empty, so matched spectra got `own=0, coiso=0, own_decoy=0, coiso_decoy=0` — four invented measurements from absent reference data. `Path("")` → `"."` is gone with it. | 3 |
| **G-C09** | CHECKED | The DIA-NN report's own `Run`/`File.Name` column must name **this sample's** raw file; NaN `Q.Value` aborts. Multi-run reports are **filtered, not rejected** — a rule that refused them would fire on normal DIA-NN output and push users to hand-make copies outside the harness. Pinning a digest makes a wrong reference *stable*; this is what makes it *detected*. | 5 |
| **G-C10** | CHECKED | 2-D (m/z × 1/K0) window assignment from the per-sample table. Zero matches → `OUT_OF_WINDOW_TABLE`; multiple matches → `AMBIGUOUS_WINDOW`. **No first-match fallback** (v1 resolved overlaps by taking the first tile). Above 0.5% of MS2 spectra, the affected metrics become `METHOD_MISMATCH` and the report refuses. | 6 |
| **G-C11** | CHECKED | Every spectrum leaves the loop through a named disposition; the counters must partition exactly, **cross-checked against the mzML's own `spectrumList count`** — a self-check alone would pass on a truncated file because both sides shrink together. | 5 |
| **G-C12** | CHECKED | One `attribute()` returning exactly one category, for targets and decoys alike, with `own + coiso + unexplained == top-K total` asserted per spectrum. v1 used `if/elif` for targets and two independent `if`s for decoys, so a peak matching both sets counted once for targets and twice for decoys. | 2 |
| **G-C13** | STRUCTURAL | One `select_match()` (single best) used for the target pass **and** the control pass. v1 selected one best candidate for targets but added *every* candidate within tolerance for decoys — different counting rules for the value and the floor that corrects it. | 1 |
| **G-C14** | CHECKED | The decoy control's **opportunity count** is recorded. Zero matches across the whole run, or too few opportunities, → `CONTROL_DEGENERATE`, which blocks publication. Also refuses when >20% of decoy queries shift outside their own tile (positional non-exchangeability). This is the check that would have caught the retracted 97.8% number. | 4 |
| **G-C15** | CHECKED | Abbott correction is **never clamped**: below-floor → `OUT_OF_DOMAIN`, floor >50% → `CONTROL_TOO_HIGH`, result outside [0,100] → `OUT_OF_DOMAIN`. v1's `max(1-recd, 1e-9)` published negative recalls and 500% recalls without complaint. | 4 |
| **G-C16** | STRUCTURAL | Numerator and denominator populations are named and checked against a declared subset relation. v1 computed charge agreement over *matched spectra* and its majority baseline over *all truth precursors*, then subtracted. | 3 |
| **G-C18** | STRUCTURAL | `Delta` is constructible only from a `Pair`; `Pair` carries one sample and refuses to span two. **FAILURE 1** — the worst published error in this project — was a cross-sample ratio (7,430 was dataset D, divided by dataset A's the reference implementation reference). A percent metric yields percentage **points** only; relative change of a percentage is `NOT_COMPUTED`. | 3 |
| **G-C19** | CHECKED | Pair validation: same runset, host, converter digest, Sage binary/config/FASTA, `method_id`, `estimand_id`; differing INIs; both members complete; baseline present. **FAILURE 2 / FAILURE 7**: v1 recorded binary and host in every manifest and never compared them, then printed only the first row's values for the whole sample. | 8 |
| **G-C20** | STRUCTURAL | `render()` has **no branch that returns a string for an `Undefined`** — it raises. Any refusal blocks the whole report; the incomplete report goes to `benchmarks-incomplete/`, never the vault, with `UNDEFINED(CODE)` in the cell. v1 rendered a dash next to real numbers. | 9 |
| **G-C21** | CHECKED | `estimand_id` — a digest over the **metrics source code** plus every tolerance, `n_min`, `RT_BIN`, decoy model and `method_id` — must be equal across a column, and a delta across differing estimands is refused. The `mz_source_profile` is bucketed, not carried as a raw fraction, so the gate cannot mismatch forever on 0.998 vs 0.997. | 1 |
| **G-C22** | VERIFIED | Chain binding at read time: search's input digest must equal this cell's mzML digest; metrics' input digests must equal both. Plus a size sweep over every sealed artifact. Bytes from three attempts cannot satisfy three digest equalities in one cell. | 2 |
| **G-C23** | STRUCTURAL | Vault write is `O_EXCL` and the filename carries the runset id. Two plans sharing a filename stem overwrote each other's report in v1. | 3 |
| **G-C30** | CHECKED | Token-aware lint over the measurement region bans `.get(` with a default, `or [0]`, `or 1.0`, `max(x,1)`, `max(x,1e-9)`, 3-argument `getattr`, `json.dumps`, and bare `except:`. **Scope is the region only — see UNCOVERED.** | 2 |

---

## The selftest is the deliverable

> The v1 harness was committed without a single guard ever having been demonstrated to fire.
> Six turned out to be theatre. **A guard that has never been observed to abort is a comment.**

`--selftest` builds a synthetic world — a real sqlite `analysis.tdf` with a
`DiaFrameMsMsWindows` table, real parquet references, real mzML files parsed by pyteomics, a
stand-in converter and a stand-in Sage whose behaviour is driven by environment variables — and
then deliberately triggers each guard.

Two properties make it more than a smoke test:

1. **It asserts the guard id.** `Abort` carries the id, so "the right guard fired" is
   distinguishable from "something crashed and happened to prevent the bad outcome". A crash
   where a guard was expected is reported as `FAIL`, not `PASS`.
2. **Every abort path is paired with a positive control.** A harness that refuses everything gets
   bypassed by hand, and then it guarantees nothing. `expect_ok` cases assert the guard does *not*
   fire on valid data — including the case v2 was warned about, a multi-run DIA-NN report.

`collate2.py selftest` also runs **two complete pipelines** (`bench2 run` → `collate2 truth` →
`collate2 collate --publish`):

- one whose decoy control is inert — must exit `3`, must leave the vault **empty**, must write an
  incomplete report with `UNDEFINED(CONTROL_DEGENERATE)` in the cell;
- one whose control fires — must exit `0`, must publish, must contain no `UNDEFINED` cell, must
  not mention `spectrum_q` in any table row, must state `n = 1`; then re-collating must refuse to
  overwrite, `verify --deep` must pass, and appending one byte to a sealed artifact must make both
  `verify --deep` and `collate` abort.

Current status: **64/64 (bench2) and 88/88 (collate2) guard cases pass.** Writing the selftest
found two real bugs in this rewrite (a commit-order `chmod` that made `os.rename` fail, and an INI
artifact sealed with its canonical digest instead of its file digest) that code review had not.

---

## UNCOVERED

Failures this harness does **not** prevent. They are listed here rather than described in a
docstring as if they were handled.

### Measurement validity (the harness enforces bookkeeping, not estimands)

1. **The decoy model is not calibrated.** `+11.003 Da` is versioned, digested, and its opportunity
   count published. G-C14 refuses when the control never fires. It **cannot certify the control
   when it does** — whether the shifted query is an exchangeable null remains unestablished. The
   two reviews disagreed on whether the shift can ever match at 20 ppm; v2 refuses to settle that
   by argument and instruments it instead.
2. **The decoy is positionally non-exchangeable within a tile.** A precursor in the top ~11 Da of
   its window has a shifted query that leaves the window. G-C14 *measures* this fraction and
   refuses above 20%; it does not correct for it below that.
3. **The Abbott correction's model is unverified.** `observed = true + (1-true)·floor` requires an
   independent, comparable control probability. Target and decoy here are unions of precursor IDs
   across correlated spectra sharing RT/IM/window geometry. Bounds-checked, not justified.
4. **top-20 purity is a top-20 estimand.** Named honestly, `depth_basis` recorded. It is still not
   total search utility — Sage searches every peak, and a tool can emit unlimited low-ranked
   contamination without penalty.
5. **`q_metric_count_gap_pct` is not a redundancy measure.** Renamed from v1's `fdr_loss`, which
   invited that reading. It is the gap between two q-value filters, confounded by score
   calibration, PSM multiplicity, peptide aggregation and search-space composition.
6. **No replication, no uncertainty.** Every cell is `n = 1`. Arm order is randomized within a
   sample and the seed recorded, which removes *systematic* order bias but **cannot estimate
   variance**. The `n` on a `Value` is a count of spectra, and spectra are not independent
   observations of a precursor — it is not a basis for a confidence interval, and the report says
   so.
7. **1/K0 tile bounds use a linear scan-number interpolation**, not the Bruker SDK conversion. The
   model name is in `method_id` and `estimand_id`, so numbers computed under different models
   cannot share a column — but the interpolation itself is unvalidated.
8. **`DiaFrameMsMsWindows` extraction is itself unvalidated code** sitting on the critical path.
   It is unit-tested against synthetic tables, not against a real Bruker acquisition.

### Bypasses that remain

9. **Seals are accident-evident, not tamper-evident.** The digest is unkeyed and `seal_digest()`
   ships in this repository. Anyone who can write `$BENCH_ROOT` can edit an artifact and recompute
   every digest up the chain. This resists accidental mixing and silent drift. It does **not**
   resist a determined operator, and no claim of SHA-256 hardness is made anywhere.
10. **Raw-file retarget is detected, not prevented.** G-B13 re-verifies inode and `analysis.tdf`
    digest after the run and aborts, so nothing seals. A retarget-and-retarget-back inside the run
    window still evades. There is no fd-pinned `/proc/self/fd/N` exec — that is Linux-only and this
    harness must also run on darwin.
11. **`.tdf_bin` is sampled, not hashed.** Files above 8 MiB are head+tail 1 MiB. A modified middle
    — i.e. the actual spectra — is invisible to `raw_content_id`. The coverage string travels in
    every recipe; it is never called "content identity" without it.
12. **`flock` semantics are assumed on shared storage.** Over NFSv3 without lockd, or Lustre
    without `-o flock`, the lock may be a silent no-op. The write-readback detects a broken
    lockfile, not broken exclusion.
13. **The G-C30 lint covers one region of one file** and is defeated by `d[k] if k in d else v`,
    `defaultdict`, `setdefault`, `try/except KeyError`, or string concatenation. It is a regression
    lock on known-bad idioms, not a proof.
14. **`Value._v` is private by convention**, enforced by that lint. Python cannot prevent attribute
    access; an author who wants to divide two payloads can. It stops the accident.
15. **`collate` does not rehash large artifacts.** It size-checks them (catching truncation and
    appending) and chain-checks digests recorded at truth time. A same-size edit to a sealed mzML
    is caught only by `bench2.py verify --deep`, which the report footer prints the command for.
    Nothing forces an operator to run it.

### Not prevented at all

16. **Third-party CPU, memory and I/O contention.** G-B07 excludes other *bench2* invocations. It
    does not exclude another user's job, another Sage, NUMA placement, cgroup limits or thermal
    state. Load average and available memory are RECORDED in each seal so a suspicious result can
    be interrogated afterwards.
17. **Effective thread count is not verified.** The requested count goes into the INI; nothing
    confirms the tool used it, and Sage's threading is not controlled.
18. **Output determinism is never tested.** No mode runs the same recipe twice and compares output
    digests. With 120 threads the converter may well be nondeterministic, and G-B14's
    output-digest no-op detection assumes it is not — a nondeterministic tool would make a genuine
    no-op invisible to the byte-equality half of that guard (the INI half still fires).
19. **No garbage collection or retention policy.** `stages/` is append-only forever and `failed/`
    retains partial multi-GB mzMLs. G-B19 checks free space before a run; nothing reclaims it.
20. **No retraction mechanism.** A published vault file that later turns out to be wrong cannot be
    marked retracted from inside the harness. Given this project's history, that is a real gap.
21. **Cost.** Default `run` recomputes every stage; there is no `--reuse-convert`. A 3-sample ×
    2-arm plan is a multi-hour serial job on one host. **The eleven original failures happened
    because people went ad hoc.** If the correct path is much more expensive than the incorrect
    one, it will be bypassed — and then the harness guarantees nothing. This is the most likely
    way v2 fails in practice.

---

## Migration

v1 output under `/path/to/scratch/harness/<plan>/` has no seals and no chain. `collate2.py`
cannot read it and **there is deliberately no import path**. Move those directories aside.

Every number currently in `docs/benchmarks/` and in the top-level `README.md` was produced by v1
or by ad-hoc scripts. Each must either be re-derived under v2 or carry an explicit
*no provenance chain* marker. That is a consequence of this design, and it is the correct one.

`samples.yaml` needs a v2 shape before anything runs: `bench_root`, per-sample `pinned` digests
(`bench2.py pin --sample dataset A`), references as `{path, sha256}` mappings rather than bare paths,
`search.*` as `{path, sha256}` mappings, and `tool.scaled_params`. The `fingerprint:` field is
removed by `pin` — a substring is not an identity, and believing otherwise is FAILURE 1.
