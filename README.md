# SpeXtractor

Extract pseudo-DDA ("pseudo-MS/MS") spectra from Bruker timsTOF **diaPASEF** data.

SpeXtractor turns a data-independent acquisition into a set of pseudo-tandem spectra that any
ordinary DDA search engine can read. Each emitted spectrum pairs one precursor hypothesis with the
fragment traces that co-elute with it in both retention time and ion mobility. The output is plain
mzML, so it feeds Sage, MSFragger, Comet or anything else that reads mzML — no library, no
spectral prediction, no enumeration of the search space in advance.

That last point is the reason the tool exists: because the search space is not enumerated before
the search runs, an **open or blind search** over the output can report variants and unexpected
modifications that a library-based DIA workflow cannot represent.

SpeXtractor is BSD-3-Clause and is a standalone application. **OpenMS is a prerequisite, not a
host** — the tool links against an OpenMS installation but lives outside the OpenMS source tree.

## Requirements

| | |
|---|---|
| OpenMS | 3.6.0 (verified against `3.6.0-pre-exported-20260717`), built `WITH_OPENTIMS` to read Bruker `.d` |
| Compiler | C++20 — GCC 13+ or Clang 16+ |
| OpenMP | strongly recommended; without it the tool runs single-threaded |
| CMake | 3.16+ |
| **Memory** | **80–125 GB peak** for a typical 30–60 min gradient, and it scales with the acquisition |
| Disk | output mzML is roughly 6–9 GB per run |

The memory figure is not a suggestion. A full diaPASEF acquisition holds millions of mass traces in
flight; on a 32 GB machine the tool will not complete. Peak usage is reported at the end of every
run, so you can size a node from your own data.

## Install

```bash
git clone https://github.com/okohlbacher/spextractor.git
cd spextractor
scripts/apply_openms_patches.sh /path/to/OpenMS    # see "OpenMS patches" below
cmake -B build -DOpenMS_DIR=/path/to/OpenMS/build
cmake --build build -j
```

The binary is `build/spextractor`. `cmake --install build` puts it in `bin/`.

### OpenMS patches

`scripts/apply_openms_patches.sh` installs `src/TdfMzCalibration.h` into the OpenMS tree and applies
two patches. **Both matter, and the build will not tell you if they are missing:**

- **Vendor m/z calibration.** Without it the Bruker reader converts flight time to m/z with a
  two-point linear-in-sqrt chord that is **−5 to −11 ppm biased** (m/z dependent) on every file
  measured, costing roughly **6–11% of closed-search peptide identifications**. Every emitted mzML
  records which calibration was used in the `spx:mz_calibration` userParam
  (`tdf_table_modeltype1` | `bruker_sdk` | `legacy_chord_APPROXIMATE`). An unsupported calibration
  table **fails closed** rather than silently producing biased masses.
- **Lock-free elution-peak detection.** OpenMS guards a shared vector with a program-global critical
  section. Called from inside SpeXtractor's parallel window loop, that one lock serialises the tool.

Optionally, `OPENMS_BRUKER_SDK_PATH=/path/to/libtimsdata.so` uses Bruker's own library for the
conversion instead. It is an independent cross-check, not a requirement, and is not redistributed.

### Optional: `.mzpeak` input

Build the mzPeak C++ library once against the same Arrow/Parquet/libzip that OpenMS uses
(`scripts/build_mzpeak_lib.sh`), then configure with `-DMZPEAK_ROOT=<checkout>`. `.d` remains the
primary input and is faster.

## Run

```bash
spextractor -in sample.d -out pseudo.mzML -threads 64
```

That is the complete command. Every default is the configuration the project benchmarks; a run that
needed extra flags to reproduce a published figure would be a bug, and there is a test for exactly
that.

Then search the output like any DDA file:

```bash
sage sage.json -o results pseudo.mzML
```

### Options worth knowing

| option | default | what it does |
|---|---|---|
| `-threads` | 1 | worker threads; the window loop scales to ~70 concurrent windows |
| `-trace:detector` | `integer` | `integer` works on the instrument's flight-time bin; `openms` uses OpenMS `MassTraceDetection`. Different algorithms — see below |
| `-charge:min_charge` | 2 | lowest precursor charge to emit. Singly-charged hypotheses are ~30% of emission and ~1.7% of peptides |
| `-assembly:require_isotope_support` | `true` | drop precursor hypotheses with no isotope partner. `false` roughly doubles emission, is ~7× slower, and identifies fewer peptides |
| `-perf:malloc_trim` | `true` | return the allocator's free pages at the two phase boundaries (glibc). The window loop otherwise stacks on top of everything the earlier phases freed but the allocator kept: measured **-26.6% peak** on a 2-hour acquisition, byte-identical output |
| `-perf:stream_load` | `true` | read the `.d` frame by frame. `false` holds the whole run in memory (90 GB floor) **and changes the output** |
| `-trace:max_span_sec` | 120 | trim a mass trace to this many seconds around its apex |

`spextractor --help` lists every option with the measurement behind each default.

### Reading the output

Emitted spectra are MS2 with a synthetic precursor. Provenance is recorded as userParams on the
run — `spx:detector`, `spx:require_isotope_support`, `spx:mz_calibration` — so a file can always be
attributed to the configuration that produced it.

## Two detectors

`-trace:detector` selects between genuinely different algorithms, not two implementations of one.

- **`integer`** (default) finds a peak's candidate traces by arithmetic on the instrument's own
  flight-time index, and never converts its compact store back to double m/z. It needs the vendor
  calibration and falls back to `openms` — loudly — without it.
- **`openms`** runs OpenMS `MassTraceDetection` on a materialised peak map.

They agree on about 85% of the union of identified peptides. Which one identifies more depends on
the search engine and on the file, so **the reason `integer` is the default is memory** — roughly
40% less — not a peptide gain. If you are chasing identifications on a particular dataset, it is
worth trying both.

## Performance

Measured on a 30-minute gradient (33,553 frames), 100 threads on a 128-core node, shipped defaults:

| | |
|---|---|
| wall time | ~5:15 |
| peak RSS | ~80 GB |
| window-loop occupancy | 56–70× |
| pseudo-spectra emitted | ~0.66 M |

Preloading tcmalloc is worth more than any structural memory work in this repository (measured
−18.4% peak RSS).

Runtime is dominated by the window loop, which is parallel across isolation windows and, within a
window, across flight-time bands. The loader is a few percent. Memory, not CPU, sets how many
windows can be in flight: an admission gate sizes concurrency from free RAM at startup.

## Determinism

The spectrum data is **byte-identical at any thread count and any loader batch size** for a fixed
binary (verified on three independent pairs: 8 vs 100 threads twice, loader batch 64 vs 256).
`bench/semantic_digest.py` hashes `<spectrumList>`..`</spectrumList>`, excluding the wall-clock
stamp, the recorded parameters and the trailing byte-offset index — none of which can match across
invocations.

This does **not** survive a code change: a last-ulp difference in one arithmetic path cascades to
~2% in peptide counts. Compare two builds by the accepted peptide **set** and the ppm median, never
by raw counts. The two statements are about different things and both have been measured.

## Tests

```bash
python3 test/test_spextractor.py build/spextractor
```

Thirteen end-to-end checks against a synthetic acquisition: isotope ownership, the charge floor,
parameter validation, thread-count determinism, mass-trace extension across missing cycles, and
that the shipped detector is the one that actually runs. No fixtures, no framework, no network.

`tests/` additionally holds the calibration unit tests and the entrapment-FDR scoring used for the
benchmark record.

## Layout

```
src/        the tool (one translation unit) plus the calibration model
scripts/    build helpers, OpenMS patch application, cluster deployment
test/       the end-to-end suite
tests/      calibration and entrapment tests, golden calibration values
bench/      benchmark drivers and the semantic digest
docs/       measurements, design notes and review records
evidence/   frozen run records with checksums
```

## Documentation

- [docs/dataset D-BASELINE.md](docs/dataset D-BASELINE.md) — the decision record: every default, what was
  measured to justify it, and what was falsified
- [docs/MZ-AXIS-DESIGN.md](docs/MZ-AXIS-DESIGN.md) — the flight-time index and the integer detector
- [docs/charge-inference.md](docs/charge-inference.md) — charge assignment, its coupling to MS1
  splitting, and the levers that did not work
- [CHANGELOG.md](CHANGELOG.md)

## Citing

See [CITATION.cff](CITATION.cff).

## License

BSD-3-Clause. Derived from OpenMS (BSD-3-Clause) — see [NOTICE](NOTICE).
