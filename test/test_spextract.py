#!/usr/bin/env python3
"""End-to-end checks on a synthetic input. No framework, no fixtures: one file, plain asserts.

Verification on this project has been whole-file cluster benchmarks, which catch regressions in
aggregate but cannot say WHY a number moved, and cannot run before a commit. These are the smallest
checks that fail if the things most likely to break do break:

  1. isotope ownership   one z=2 envelope (M, M+1, M+2) must yield ONE precursor at the mono, not
                         three, and the M+1/M+2 peaks must not reseed precursors of their own
  2. min_charge          a supported z=1 envelope is dropped at the default and kept at 1, while a
                         z=2 envelope is unaffected either way
  3. parameter bounds    a negative count must be REJECTED, not wrapped through Size into an
                         empty-but-successful run
  4. determinism         the same input at 1 and 4 threads must give byte-identical spectrum data
  5. trace extension     one missing cycle must leave ONE trace and its gapped fragment in the
                         spectrum; seven missing cycles must give TWO traces. Every version of the
                         integer detector that lost 92% of the peptides passed checks 1-4.

These exercise the default (OpenMS) detector only. `trace:detector=integer` refuses to run
without the vendor flight-time calibration, which a synthetic mzML cannot carry, so it is judged on
real data in the benchmark harness and not here.

Usage: test_spextract.py /path/to/spextract [workdir]
Exit status is 0 only if every check passes.
"""
import base64, hashlib, os, re, struct, subprocess, sys, tempfile

ISO = 1.0033548
PROTON = 1.007276


def b64f(vals, double=True):
    fmt = "<%d%s" % (len(vals), "d" if double else "f")
    return base64.b64encode(struct.pack(fmt, *vals)).decode()


def spectrum(idx, ms_level, rt, peaks, im, prec=None):
    """One mzML spectrum. `peaks` is [(mz, intensity)], `im` one 1/K0 per peak."""
    mz = b64f([p[0] for p in peaks])
    inten = b64f([p[1] for p in peaks], double=False)
    imarr = b64f(im)
    pre = ""
    if prec:
        target, lo, hi = prec
        pre = f"""<precursorList count="1"><precursor><isolationWindow>
<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{target}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{target-lo}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{hi-target}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
</isolationWindow><activation><cvParam cvRef="MS" accession="MS:1000044" name="dissociation method"/></activation></precursor></precursorList>"""
    # "frame=<N>" is the vendor key frameIdOf() parses; without it every frame is unmappable and
    # the integer detector falls back, which is how all 8 checks silently ran on the OTHER detector.
    return f"""<spectrum id="frame={idx + 1} scan=0" index="{idx}" defaultArrayLength="{len(peaks)}">
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{ms_level}"/>
<cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum"/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination"/>
<scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt}" unitAccession="UO:0000010" unitName="second" unitCvRef="UO"/></scan></scanList>
{pre}
<binaryDataArrayList count="3">
<binaryDataArray encodedLength="{len(mz)}"><cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{mz}</binary></binaryDataArray>
<binaryDataArray encodedLength="{len(inten)}"><cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{inten}</binary></binaryDataArray>
<binaryDataArray arrayLength="{len(im)}" encodedLength="{len(imarr)}"><cvParam cvRef="MS" accession="MS:1002816" name="mean inverse reduced ion mobility array" unitAccession="MS:1002814" unitName="volt-second per square centimeter" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{imarr}</binary></binaryDataArray>
</binaryDataArrayList></spectrum>"""


def synth(path, n_cycles=14, cycle_s=1.4, gap_cycles=()):
    """Two precursors that co-elute in one isolation window: a z=2 envelope and a z=1 envelope.

    Both are given a 3-peak isotope envelope so `require_isotope_support` keeps them, matching
    fragments that co-elute with the precursor, and a distinct ion mobility.

    `gap_cycles`: cycles in which the z=2 precursor AND its first fragment are ABSENT. This is
    the input that separates a correct trace-extension rule from a wrong one: OpenMS tolerates up
    to five consecutive missing frames before it ends a trace, so one gap must still give ONE
    trace, and seven consecutive gaps must give TWO. Every version of the integer detector that
    lost 92% of peptides passed the other checks in this file; none would have passed this.
    """
    win = (600.0, 590.0, 620.0)                    # target, lo, hi
    z2_mono, z2_im = 601.3, 0.90                   # a doubly-charged precursor
    z1_mono, z1_im = 610.6, 1.20                   # a singly-charged one
    z2_frags = [(233.11, 0.7), (348.19, 1.0), (461.27, 0.55)]
    z1_frags = [(288.14, 0.8), (401.22, 0.9)]
    specs, idx = [], 0
    for c in range(n_cycles):
        rt = 1.0 + c * cycle_s
        shape = 1.0 - abs(c - n_cycles / 2.0) / (n_cycles / 2.0 + 1.0)   # a triangular elution
        amp = max(shape, 0.05)
        ms1, ms1_im = [], []
        gap = c in gap_cycles
        for k, rel in enumerate((1.0, 0.55, 0.20)):                      # averagine-ish envelope
            if not gap:
                ms1.append((z2_mono + k * ISO / 2.0, 6.0e5 * amp * rel)); ms1_im.append(z2_im)
            ms1.append((z1_mono + k * ISO / 1.0, 4.0e5 * amp * rel)); ms1_im.append(z1_im)
        order = sorted(range(len(ms1)), key=lambda i: ms1[i][0])
        specs.append(spectrum(idx, 1, rt, [ms1[i] for i in order], [ms1_im[i] for i in order])); idx += 1
        ms2, ms2_im = [], []
        for fi, (mz, rel) in enumerate(z2_frags):
            if gap and fi == 0: continue                                   # first fragment gapped too
            ms2.append((mz, 2.0e5 * amp * rel)); ms2_im.append(z2_im)
        for mz, rel in z1_frags: ms2.append((mz, 1.5e5 * amp * rel)); ms2_im.append(z1_im)
        order = sorted(range(len(ms2)), key=lambda i: ms2[i][0])
        specs.append(spectrum(idx, 2, rt, [ms2[i] for i in order], [ms2_im[i] for i in order],
                              prec=win)); idx += 1
    body = "\n".join(specs)
    open(path, "w").write(f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml"><mzML version="1.1.0" id="spextract_test">
<cvList count="1"><cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/></cvList>
<fileDescription><fileContent><cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum"/></fileContent></fileDescription>
<softwareList count="1"><software id="so_test" version="0"><cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="spextract-test"/></software></softwareList>
<instrumentConfigurationList count="1"><instrumentConfiguration id="ic_test"><cvParam cvRef="MS" accession="MS:1000031" name="instrument model"/></instrumentConfiguration></instrumentConfigurationList>
<dataProcessingList count="1"><dataProcessing id="dp_test"><processingMethod order="0" softwareRef="so_test"><cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML"/></processingMethod></dataProcessing></dataProcessingList>
<run id="run_test" defaultInstrumentConfigurationRef="ic_test">
<spectrumList count="{len(specs)}" defaultDataProcessingRef="dp_test">
{body}
</spectrumList></run></mzML></indexedmzML>
""")
    return z2_mono, z1_mono


def synth_tdf(path, n_frames):
    """A minimal analysis.tdf: the real dataset D MzCalibration row plus a Frames table.

    The calibration constants are the measured dataset D ones from tests/calibration_golden.json, not
    invented numbers -- an invented row would either fail isSupported() or silently define a
    different mass scale. T1 is given a small per-frame spread so the per-frame factor is actually
    exercised rather than collapsing to the reference.
    """
    import sqlite3, json, os
    g = json.load(open(os.path.join(os.path.dirname(__file__), "..", "tests",
                                    "calibration_golden.json")))[0]
    if os.path.exists(path):
        os.remove(path)
    db = sqlite3.connect(path)
    db.execute("CREATE TABLE MzCalibration (Id INTEGER, ModelType INTEGER, DigitizerTimebase REAL,"
               " DigitizerDelay REAL, C0 REAL, C1 REAL, C2 REAL, T1 REAL, dC1 REAL, dC2 REAL,"
               " C3 REAL, C4 REAL)")
    db.execute("INSERT INTO MzCalibration VALUES (1,?,?,?,?,?,?,?,?,0,0,0)",
               (g["model_type"], g["timebase"], g["delay"], g["C0"], g["C1"], g["C2"],
                g["T1_ref"], g["dC1"]))
    # Real tdfs reference the calibration row per frame (Frames.MzCalibration); mirror that, since
    # the loader selects the row the frames reference rather than "the only row".
    db.execute("CREATE TABLE Frames (Id INTEGER, T1 REAL, MzCalibration INTEGER)")
    db.executemany("INSERT INTO Frames VALUES (?,?,1)",
                   [(i + 1, g["T1_ref"] + 0.03 * (i % 3)) for i in range(n_frames)])
    db.commit(); db.close()
    return path


def run(binary, inp, out, extra=(), threads=1, expect_fail=False, env=None):
    cmd = [binary, "-in", inp, "-out", out, "-threads", str(threads),
           "-assembly:min_fragments", "2", "-assembly:require_isotope_support", "true",
           "-trace:ms2_min_length_sec", "0", "-trace:ms1_split_valleys", "0",
           "-trace:ms2_split_valleys", "0"] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if expect_fail:
        assert r.returncode != 0, "expected a non-zero exit, got success:\n" + r.stdout[-800:]
        return r
    assert r.returncode == 0, "run failed:\n" + (r.stdout + r.stderr)[-2000:]
    return r


def precursors(mzml):
    """(mz, charge) of every emitted spectrum."""
    out, mz = [], None
    for line in open(mzml, errors="replace"):
        if "selected ion m/z" in line: mz = float(re.search(r'value="([^"]+)"', line).group(1))
        elif "charge state" in line and mz is not None:
            out.append((mz, int(re.search(r'value="([^"]+)"', line).group(1)))); mz = None
    return out


def peaks_of(mzml, prec_mz, tol=0.02):
    """m/z list of the first spectrum whose selected ion is within tol of prec_mz."""
    want, arr, prec = False, None, 8
    for line in open(mzml, errors="replace"):
        if "<spectrum " in line: want = False
        elif "selected ion m/z" in line:
            v = float(re.search(r'value="([^"]+)"', line).group(1)); want = abs(v - prec_mz) < tol
        elif want and 'name="m/z array"' in line: arr = "mz"
        elif want and 'name="intensity array"' in line: arr = "it"
        elif want and 'name="64-bit float"' in line: prec = 8
        elif want and 'name="32-bit float"' in line: prec = 4
        elif want and arr == "mz" and "<binary>" in line:
            b = line[line.find("<binary>") + 8: line.find("</binary>")]
            raw = base64.b64decode(b)
            return list(struct.unpack(f"<{len(raw)//prec}{'d' if prec == 8 else 'f'}", raw))
    return []


def ms1_traces_near(tsv, mz, tol=0.02):
    """(n_xic, ...) of every MS1 trace in a diag:dump_ms1_tsv dump near an m/z."""
    out = []
    with open(tsv) as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            q = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if abs(float(q["mz"]) - mz) < tol: out.append(int(q["n_xic"]))
    return out


def digest(mzml):
    """Hash the spectrum data only -- the header carries a wall-clock stamp and parameters."""
    h, on = hashlib.sha256(), False
    for line in open(mzml, "rb"):
        if b"<spectrumList" in line: on = True
        if on: h.update(line)
        if b"</spectrumList>" in line: break
    return h.hexdigest()


def main():
    if len(sys.argv) < 2: sys.exit(__doc__)
    binary = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="spextract_test_")
    os.makedirs(work, exist_ok=True)
    inp = os.path.join(work, "synth.mzML")
    z2_mono, z1_mono = synth(inp)
    print(f"synthetic input: {inp}  (z=2 mono {z2_mono}, z=1 mono {z1_mono})")
    fails = []

    def check(name, fn):
        try:
            fn(); print(f"  PASS  {name}")
        except AssertionError as e:
            print(f"  FAIL  {name}: {e}"); fails.append(name)

    # 1 + 2: the default drops the z=1 envelope and keeps one precursor per envelope
    out_def = os.path.join(work, "default.mzML")
    run(binary, inp, out_def)
    pd = precursors(out_def)

    def c1():
        z2 = [p for p in pd if p[1] == 2]
        assert z2, f"no z=2 precursor emitted; got {pd}"
        assert all(abs(m - z2_mono) < 0.02 for m, _ in z2), \
            f"a z=2 precursor is not at the monoisotope {z2_mono}: {z2} -- an M+1/M+2 peak reseeded"
    check("isotope ownership: one z=2 precursor, at the mono", c1)

    def c2():
        assert not [p for p in pd if p[1] == 1], \
            f"charge:min_charge=2 is the default but a z=1 precursor was emitted: {pd}"
    check("min_charge default 2 drops the z=1 envelope", c2)

    out_z1 = os.path.join(work, "z1.mzML")
    run(binary, inp, out_z1, extra=["-charge:min_charge", "1"])
    p1 = precursors(out_z1)

    def c3():
        assert [p for p in p1 if p[1] == 1], f"charge:min_charge=1 must emit the z=1 envelope: {p1}"
        assert len([p for p in p1 if p[1] == 2]) == len([p for p in pd if p[1] == 2]), \
            "lowering min_charge changed the z=2 result"
    check("min_charge 1 restores it and leaves z=2 alone", c3)

    # 3: a negative count must be rejected, not wrapped through Size
    def c4():
        run(binary, inp, os.path.join(work, "bad.mzML"),
            extra=["-assembly:min_fragments", "-1"], expect_fail=True)
    check("a negative min_fragments is rejected", c4)

    def c5():
        run(binary, inp, os.path.join(work, "bad2.mzML"),
            extra=["-gate:delta_rt", "-5"], expect_fail=True)
    check("a negative gate:delta_rt is rejected", c5)

    # 4: determinism across thread counts for a fixed binary
    out_t4 = os.path.join(work, "t4.mzML")
    run(binary, inp, out_t4, threads=4)

    def c6():
        a, b = digest(out_def), digest(out_t4)
        assert a == b, f"spectrum data differs between 1 and 4 threads:\n  {a}\n  {b}"
    check("byte-identical at 1 vs 4 threads", c6)

    # 5: trace extension across gaps -- the check every broken integer detector would have failed
    # 20 cycles, so that after a seven-cycle gap BOTH halves are long enough to pass the MS1
    # minimum trace length (a 3-frame remnant is correctly dropped by OpenMS, which is not the
    # behaviour under test here).
    for label, gaps, want_traces in (("one missing cycle -> ONE trace", (10,), 1),
                                     ("seven missing cycles -> TWO traces", tuple(range(7, 14)), 2)):
        inp_g = os.path.join(work, f"gap{len(gaps)}.mzML")
        synth(inp_g, n_cycles=20, gap_cycles=gaps)
        out_g = os.path.join(work, f"gap{len(gaps)}.out.mzML")
        dump = os.path.join(work, f"gap{len(gaps)}")
        run(binary, inp_g, out_g, extra=["-diag:dump_ms1_tsv", dump])

        def cg(label=label, want=want_traces, dump=dump, out_g=out_g, gaps=gaps):
            tr = ms1_traces_near(dump + ".traces.tsv", z2_mono)
            assert len(tr) == want, f"expected {want} MS1 trace(s) at {z2_mono} with gaps {gaps}, got {len(tr)} (n_xic {tr})"
            if want == 1:
                assert tr[0] >= 20 - len(gaps) - 1, f"the surviving trace is too short: n_xic {tr[0]}"
                # the gapped FRAGMENT must still reach the z=2 spectrum: a rule that ends a trace at
                # its first gap drops it, which is how 92% of the peptides went missing
                pk = peaks_of(out_g, z2_mono)
                assert any(abs(m - 233.11) < 0.02 for m in pk), f"gapped fragment 233.11 missing from the z=2 spectrum: {pk}"
        check(f"trace extension: {label}", cg)

    # 9 + 10: the SHIPPED detector, and the calibration failing closed.
    # Until 2026-09-04 every check above ran on spx:detector=openms -- the fallback -- because the
    # synthetic nativeIDs carried no "frame=" key and there was no tdf, so the integer detector
    # (the default since that morning) was never once exercised by the suite.
    def detector_of(path):
        import re
        m = re.search(r'spx:detector"[^>]*value="([a-z]+)"', open(path, encoding="utf-8",
                                                                  errors="replace").read())
        return m.group(1) if m else None

    def c9():
        tdf = synth_tdf(os.path.join(work, "analysis.tdf"), n_frames=64)
        out_i = os.path.join(work, "integer.mzML")
        env = dict(os.environ, SPEXTRACT_MZPEAK_TDF=tdf)
        run(binary, inp, out_i, env=env)
        got = detector_of(out_i)
        assert got == "integer", (
            f"the SHIPPED default detector did not run: spx:detector={got}. With a parseable "
            f"frame= nativeID and a valid tdf there is nothing left to fall back for.")
    check("the shipped detector (integer) actually runs", c9)

    def c10():
        # No tdf: the integer detector cannot calibrate, and must SAY so rather than invent one.
        out_f = os.path.join(work, "fallback.mzML")
        r = run(binary, inp, out_f, extra=("-trace:detector", "integer"))
        assert detector_of(out_f) == "openms", "no calibration available, yet integer still ran"
        assert "alling back" in (r.stdout + r.stderr), \
            "fell back to the OpenMS detector without warning -- a silent detector switch"
    check("no calibration -> falls back to openms, loudly", c10)

    total = 6 + 2 + 2
    print(f"\n{total - len(fails)}/{total} checks passed" + (f"; FAILED: {', '.join(fails)}" if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
