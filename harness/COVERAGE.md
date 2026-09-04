# Measured guard coverage (mutation testing)

**Not a claim — a measurement.** Produced by `mutate.py`, which disables each `die("G-…")`
call in turn and reruns the module's selftest.

```
KILLED    selftest failed  -> the guard is genuinely exercised
SURVIVED  selftest passed  -> the guard can be deleted and nothing notices
```

Reproduce:

```bash
python3 mutate.py bench2.py   --json /tmp/mut_bench.json
python3 mutate.py collate2.py --json /tmp/mut_collate.json
```

## Why this file exists

Two harness generations shipped guards that were documented, "tested", and inert.

v1 had six such guards; an external review classified them as theatre. v2 replaced them,
added a 64/88-case selftest, reported **64 passed, 0 FAILED** — and reproduced the same
failure. Disabling the content-digest comparison at `bench2.py:1387` (the guard for
**FAILURE 1**, the cross-sample comparison that was actually published as a wrong result)
left the suite fully green.

The selftest case for that guard *reimplements the comparison in its own body and calls
`die()` itself*. The production code is never reached.

**A passing selftest is not evidence a guard works.** The only acceptance criterion is:
delete the guard, and the suite must fail.

## Measured, 2026-07-22

| module | sites | KILLED | coverage | pure theatre |
|---|---|---|---|---|
| `bench2.py` | 119 | **43** | **36%** (was 28%) | G-B00 only (was G-B00 + G-B13) |
| `collate2.py` | 65 | 30 | 46% | G-C00, G-C10 |

**Critical gate: PASSES for `bench2.py`** (4 entries). `G-B05` (cross-sample, FAILURE 1) and
`G-B13` (raw TOCTOU, FAILURE 7) are now killed at every gated site by tests that call the
PRODUCTION path. `G-B13` went 0/6 -> 6/6.

Still open: `collate2.py G-C10` (the `WINDOWS_S23` replacement) has no test, and `collate2.py`
has no `match=` support so its 46% carries the same-id ambiguity below.

**109 of 184 guard sites remain unverified claims.** The harness is still NOT adoptable as
the sole sanctioned benchmarking path: the critical guards are verified, but the bulk is not,
and `collate2.py` has had no work at all. `bench.py` v1 remains in use with its limitations
known rather than disguised.

### Guards whose every site survived

| guard | what it claims | assessment |
|---|---|---|
| `G-B00` | pyyaml missing → abort | **acceptable untested.** Fires at import; the real failure mode is a loud `ImportError`. Not domain logic. |
| `G-C00` | pyarrow/pyteomics missing → abort | same. |
| `G-B13` | raw `.d` missing / retargeted | **FIXED: 6/6**, every site message-pinned against the production `verify_raw_unchanged`. |
| `G-C10` | produce stage carries no acquisition window table | **STILL A REAL GAP.** The guard replacing v1's hardcoded `WINDOWS_S23` (FAILURE 8). Lines 1382, 1530, 2126 all survive. Next item. |

## The subtler defect: same-id fallthrough

`expect_abort` originally matched on **guard id alone**. Mutation showed that disabling
`raw not found` at `bench2.py:426` still left its assertion green — a *sibling*
`die("G-B13")` further down fired instead, reporting an unrelated condition
(`raw member analysis.tdf changed`). The suite could not tell the two apart.

So a SURVIVED site is ambiguous:

* genuinely untested, **or**
* tested, but masked by a sibling sharing its id — in which case that specific check can be
  deleted, detection of a distinct condition is lost, and the suite stays green.

**Same-id fallthrough is itself a theatre mechanism.** `expect_abort` now takes `match=`,
which pins an assertion to the message and makes each site independently killable. It
immediately caught a wrong assumption in the first test written against it: a plain file
does not abort with "not a directory" (no such check exists) but falls through to
`contains no files` at `:447`.

**Every new guard test must pass `match=`.** Without it the test cannot distinguish which
check fired, and coverage numbers overstate protection.

## Gate

```bash
python3 mutate.py bench2.py --critical critical_guards.txt
```

Exits 3 if any listed guard survives. The critical set is guards covering failures that
have **already happened in this project**, held to 100% mutation coverage; everything else
is reported but not gated. Claiming protection everywhere and verifying it nowhere is what
produced two generations of theatre.

**Status: PASSES for `bench2.py` (4 entries), and the `collate2.py G-C10` entry is still
unmet.** Running the gate against `collate2.py` reports `CRITICAL LIST IS STALE` because that
guard has no test — which is the intended behaviour, not a bug.

## Instrument limitations

* `mutate.py` replaces the `die()` **call**, not the detection logic. A guard whose
  condition is wrong but whose `die()` is reachable still shows KILLED.
* One mutant at a time; interacting guards are not tested jointly.
* Requires a green baseline — otherwise every mutant appears KILLED and the sweep reports
  perfect coverage. `mutate.py` aborts rather than reporting that.
* The AST extractor must match this interpreter (`python3` here is **3.7.4**, which has
  neither `ast.Constant` for strings nor `Call.end_lineno`). Matching only `ast.Constant`
  silently found **zero** sites — an instrument reporting perfect coverage because it
  measured nothing. `mutate.py` now cross-checks against a raw `die("G-` count and refuses
  to report if the two disagree.
* `collate2.py` is not yet patched with `match=` support; its coverage figure carries the
  same-id ambiguity described above.


## Gate design: two bugs found in the gate itself

The first version reported **"CRITICAL GATE PASSED (6 entries verified)"** while checking two.

1. **Line-pinned entries go stale silently.** Adding tests above a guard shifted every line
   below it; `G-B05:1387` then matched no site and passed. Entries now pin to a **message
   fragment**, which survives edits.
2. **Cross-module entries were skipped.** `G-C10` lives in `collate2.py`; running the gate on
   `bench2.py` found no such site and passed it. An entry matching nothing is now an **ERROR**
   (`CRITICAL LIST IS STALE`), not a pass.

Third instance of this pattern in one session — in the tool built to detect the pattern.
