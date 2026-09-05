<p align="center">
  <img src="assets/logo.svg" alt="SpeXtract" width="132">
</p>

<h1 align="center">SpeXtract</h1>

<p align="center">
  <em>Pseudo-DDA spectra from Bruker timsTOF diaPASEF data, for open and blind searching.</em>
</p>

<p align="center">
  <a href="https://github.com/okohlbacher/spextract/actions/workflows/tests.yml"><img alt="tests" src="https://github.com/okohlbacher/spextract/actions/workflows/tests.yml/badge.svg"></a>
  <a href="LICENSE"><img alt="license" src="https://img.shields.io/badge/license-BSD--3--Clause-blue.svg"></a>
  <a href="https://github.com/okohlbacher/spextract/releases/latest"><img alt="release" src="https://img.shields.io/github/v/release/okohlbacher/spextract"></a>
</p>

Extract pseudo-DDA ("pseudo-MS/MS") spectra from Bruker timsTOF **diaPASEF** data.

SpeXtract turns a data-independent acquisition into a set of pseudo-tandem spectra that any
ordinary DDA search engine can read. Each emitted spectrum pairs one precursor hypothesis with the
fragment traces that co-elute with it in both retention time and ion mobility. The output is plain
mzML, so it feeds Sage, MSFragger, Comet or anything else that reads mzML — no library, no
spectral prediction, no enumeration of the search space in advance.

That last point is the reason the tool exists: because the search space is not enumerated before
the search runs, an **open or blind search** over the output can report variants and unexpected
modifications that a library-based DIA workflow cannot represent.

SpeXtract is BSD-3-Clause and is a standalone application. **OpenMS is a prerequisite, not a
host** — the tool links against an OpenMS installation but lives outside the OpenMS source tree.

## Requirements

| | |
|---|---|
| OpenMS | 3.6.0 (verified against `3.6.0-pre-exported-20260717`), built `WITH_OPENTIMS` to read Bruker `.d` |
| Compiler | C++20 — GCC 13+ or Clang 16+ |
| OpenMP | strongly recommended; without it the tool runs single-threaded |
| CMake | 3.16+ |
| Memory | a large-memory machine: a full acquisition holds millions of mass traces in flight |

This is a server-class tool, not a laptop one. Peak memory is reported at the end of every run, so
size a node from your own data rather than from a number quoted here.

## Install

```bash
git clone https://github.com/okohlbacher/spextract.git
cd spextract
scripts/apply_openms_patches.sh /path/to/OpenMS    # see "OpenMS patches" below
cmake -B build -DOpenMS_DIR=/path/to/OpenMS/build
cmake --build build -j
```

The OpenMS tree must be checked out at the commit named in `patches/openms.lock` before patching —
the patches are cut against exactly that commit, and are verified to reproduce the benchmarked tree
from it byte for byte:

```bash
git clone https://github.com/OpenMS/OpenMS.git
git -C OpenMS checkout $(sed -n 's/^OPENMS_BASE=//p' patches/openms.lock)
```

The binary is `build/spextract`. `cmake --install build` puts it in `bin/`.

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
  section. Called from inside SpeXtract's parallel window loop, that one lock serialises the tool.

Optionally, `OPENMS_BRUKER_SDK_PATH=/path/to/libtimsdata.so` uses Bruker's own library for the
conversion instead. It is an independent cross-check, not a requirement, and is not redistributed.

### Optional: `.mzpeak` input

Build the mzPeak C++ library once against the same Arrow/Parquet/libzip that OpenMS uses
(`scripts/build_mzpeak_lib.sh`), then configure with `-DMZPEAK_ROOT=<checkout>`. `.d` remains the
primary input and is faster.

## Run

```bash
spextract -in sample.d -out pseudo.mzML
```

That is the complete command. Threads default to every core on the machine, and every other default
is the configuration the project benchmarks; a run that needed extra flags to reproduce a published
figure would be a bug, and there is a test for exactly that.

Then search the output like any DDA file:

```bash
sage sage.json -o results pseudo.mzML
```

### Output format

**The extension decides.** `.mzpeak` is the default when the name does not say: it is a columnar
archive (Parquet in a zip) and comes out at roughly a third the size of the equivalent mzML.

**But DDA search engines read mzML, not mzPeak**, and OpenMS' own `FileConverter` does not accept
mzPeak as an input type either. So write `.mzML` — or pass `-out_type mzML` — whenever the next step
is a search, which in practice is most of the time. The tool reports the format it wrote and warns
when that format is mzPeak.

### Options worth knowing

| option | default | what it does |
|---|---|---|
| `-threads` | all cores | worker threads. Set it explicitly on a shared machine |
| `-trace:detector` | `integer` | `integer` works on the instrument's flight-time bin; `openms` uses OpenMS `MassTraceDetection`. Different algorithms — see below |
| `-charge:min_charge` | 2 | lowest precursor charge to emit. Singly-charged hypotheses are ~30% of emission and ~1.7% of peptides |
| `-assembly:require_isotope_support` | `true` | drop precursor hypotheses with no isotope partner. `false` roughly doubles emission, is ~7× slower, and identifies fewer peptides |
| `-perf:stream_load` | `true` | read the `.d` frame by frame. `false` holds the whole run in memory, **changes the output, and on a large acquisition finds ~2% more peptides at the same measured FDR** — the more sensitive setting if you can afford ~1.75× the memory |
| `-trace:max_span_sec` | 120 | trim a mass trace to this many seconds around its apex |
| `-out_type` | from extension | force `mzpeak` or `mzML` instead of taking it from the `-out` name |

`spextract --help` lists every option with the measurement behind each default.

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

Window-loop occupancy reaches 56–70× on 100 threads, and a 30-minute gradient extracts in a few
minutes. Preloading tcmalloc is worth more than any structural memory work in this repository
(measured −18.4% peak RSS).

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
python3 test/test_spextract.py build/spextract
```

Eleven end-to-end checks against a synthetic acquisition: isotope ownership, the charge floor,
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
