#!/usr/bin/env python3
"""SpeXtractor benchmark harness v2 -- run identity, plan binding, artifact production.

This file produces ARTIFACTS and PROVENANCE. It computes no measurements; every number lives in
collate2.py. The separation is deliberate: the failure class this rewrite exists to kill is
"one published row assembled from three different attempts", and that is a bookkeeping problem,
not an arithmetic one.

WHY THIS FILE EXISTS
--------------------
Ad-hoc benchmarking produced eleven data-handling failures in one session. The v1 harness added
eleven guards. Two independent adversarial reviews found that six of them were THEATRE -- the
claimed invariant had no implementing code -- and that the largest hole was one no guard
addressed: there was no immutable run identity, so a failed converter left a stale mzML that was
accepted, a failed Sage left an old TSV that was rescored, and the manifest was overwritten with
fresh metadata.

The root cause of that whole class is LOCATION ADDRESSING. `outroot/{sample}__{arm}/pseudo.mzML`
is a name that means different things at different times, so a write that fails leaves the
previous meaning in place and the next reader cannot tell which write it is looking at.

The fix is not a check. It is: the tool writes into a directory that did not exist one second
earlier, whose name embeds a nonce no other process knows. There is no pre-existing pseudo.mzML
to survive a nonzero exit, because there is no pre-existing directory. See commit_stage().

GUARD DISCIPLINE
----------------
Every guard raises Abort(guard_id, message) and is exercised by `--selftest`, which deliberately
triggers it and verifies THAT SPECIFIC guard fired. A guard with no selftest case is listed in
README.md under NOT SELF-TESTED. A guard that is not implemented is listed under UNCOVERED. It is
not described here as if it existed -- that is precisely what got six v1 guards classified as
theatre.

Tiers, used honestly throughout:
  STRUCTURAL  the wrong thing cannot be expressed; there is no code path that accepts it
  VERIFIED    checked against a digest committed before the data existed
  CHECKED     a runtime abort on an observable condition
  RECORDED    measured and stored, not prevented (say so; do not imply more)

Commands
--------
    bench2.py plan-check --plan P            # validate only; no side effects, no dirs created
    bench2.py pin        --sample dataset A        # write pinned content digests into samples.yaml
    bench2.py run        --plan P            # execute; prints the runset id
    bench2.py verify     --runset R [--deep] # rehash sealed artifacts; writes nothing
    bench2.py show       --runset R
    bench2.py selftest                       # trigger every guard, assert it fires
"""
from __future__ import annotations

import argparse
import base64
import errno
import fcntl
import hashlib
import json
import os
import platform
import re
import shutil
import sqlite3
import struct
import subprocess
import sys
import tempfile
import time
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
HARNESS_VERSION = "2.0.0"
CONTRACT_VERSION = 2

# Only these tokens may appear in the converter argv. Everything that can change the output must
# travel through the INI, because the INI is what gets hashed into the recipe. A parameter passed
# on the command line would override the INI and be invisible to no-op detection. (FAILURE 5.)
ALLOWED_CONVERTER_ARGV_FLAGS = ("-ini", "-in", "-out")

ARM_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
SAMPLE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")

# Files below this size are hashed whole. Larger ones (Bruker .tdf_bin is tens of GB) are sampled
# head+tail. This is a real coverage limit and is labelled as such in every recipe that carries
# it -- see raw_identity(). It is NOT described anywhere as "content identity".
FULL_HASH_LIMIT = 8 << 20
SAMPLE_CHUNK = 1 << 20

PLAN_TOP_KEYS = {"schema", "title", "binary", "threads", "tcmalloc", "samples", "arms", "sweeps"}
PLAN_ARM_KEYS = {"name", "role", "kind", "baseline", "params", "source"}
PLAN_SWEEP_KEYS = {"param", "baseline", "points", "scale"}
PLAN_ROLES = {"baseline", "treatment", "external_baseline"}
PLAN_KINDS = {"run", "reference_mzml"}


# --------------------------------------------------------------------------------- abort plumbing
class Abort(Exception):
    """Guards abort; they never warn.

    Carrying the guard id on the exception is what makes --selftest able to assert that THE
    SPECIFIC guard fired rather than that something, somewhere, failed. The v1 harness had no
    way to distinguish "the guard worked" from "an unrelated crash happened to prevent the bad
    outcome", which is how six guards shipped without ever being demonstrated to fire.
    """

    def __init__(self, guard, message):
        super(Abort, self).__init__("[%s] %s" % (guard, message))
        self.guard = guard
        self.message = message


def die(guard, message):
    raise Abort(guard, message)


def canon(obj):
    """The single canonical serializer.

    allow_nan=False is load-bearing: json.dumps emits bare NaN/Infinity by default, which is not
    valid JSON and which silently poisons every downstream reader. Killing it at the serializer
    rather than at each call site means a NaN cannot enter an artifact at all. (Review: "json.dumps
    can emit nonstandard NaN".)
    """
    try:
        return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                          ensure_ascii=False, allow_nan=False).encode("utf-8")
    except ValueError as exc:
        # G-B20
        die("G-B20", "refusing to serialize non-finite float (%s). A NaN in an artifact becomes "
                     "a plausible wrong number downstream." % exc)


def jload(path):
    """Strict JSON read. Rejects NaN/Infinity on the way IN as well as on the way out."""
    def _reject(tok):
        die("G-B20", "%s contains non-JSON constant %r" % (path, tok))
    return json.loads(Path(path).read_text(), parse_constant=_reject)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def short(digest):
    return digest[:12]


def utc_now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def utc_stamp():
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())


def new_attempt_id():
    return "%s-%s" % (utc_stamp(), uuid.uuid4().hex[:12])


def harness_code_digest():
    """Digest over harness CODE ONLY.

    samples.yaml is deliberately EXCLUDED. If the registry were part of this digest, pinning an
    unrelated new sample would change every recipe id and retroactively declare completed runsets
    to have used "different harness code". Registry content still binds per-sample, through the
    pinned raw_content_id / acquisition_id that enter each recipe individually.
    """
    parts = []
    for name in sorted(("bench2.py", "collate2.py")):
        p = HERE / name
        if p.exists():
            parts.append({"f": name, "h": sha256_file(p)})
    return sha256_bytes(canon(parts))


# ------------------------------------------------------------------------------------ filesystem
def fsync_dir(path):
    fd = os.open(str(path), os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def assert_same_filesystem(paths):
    """G-B08 -- rename atomicity is an ASSUMPTION unless the device ids are equal.

    The whole commit protocol rests on os.rename being atomic. Put staging/ on node-local disk
    for I/O speed (the obvious cluster optimization) and rename raises EXDEV; the reflex fix is
    shutil.move, which is a 12 GB copy -- non-atomic, and the design collapses silently. Assert
    instead of assuming.
    """
    devs = {}
    for p in paths:
        devs[str(p)] = os.stat(str(p)).st_dev
    if len(set(devs.values())) != 1:
        die("G-B08", "staging/, stages/ and failed/ are on different filesystems %s.\n"
                     "       os.rename() would raise EXDEV and the atomic-commit protocol -- the "
                     "entire basis of run identity -- would not hold." % devs)


def check_free_space(path, need_gb):
    """G-B19 -- /scratch quota exhaustion silently kills jobs on this cluster.

    A 6 h convert that dies at 95%% on a full filesystem is indistinguishable, after the fact,
    from a tool bug. Refuse up front.
    """
    usage = shutil.disk_usage(str(path))
    free_gb = usage.free / (1024.0 ** 3)
    if free_gb < need_gb:
        die("G-B19", "only %.1f GB free at %s; this plan needs at least %.1f GB.\n"
                     "       Quota exhaustion mid-run leaves a truncated artifact and wastes the "
                     "whole runset." % (free_gb, path, need_gb))
    return free_gb


class HostLock(object):
    """G-B07 -- exclusive access for HARNESS INVOCATIONS on this host.

    Replaces v1's `ps` substring scan, which was race-prone (two invocations could both pass
    preflight before either launched), BSD-fragile, defeated by a renamed binary, and had a real
    substring bug: `str(os.getpid()) not in line` false-matches pid 123 against pid 5123.

    HONEST SCOPE: this excludes other bench2 invocations. It does NOT exclude an unrelated job,
    another user's Sage, or general node load. Those are RECORDED in the seal (loadavg, memory,
    big processes) so a suspicious result can be interrogated -- they are not prevented, and this
    docstring does not claim they are. See README UNCOVERED.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.fd = None

    def acquire(self, runset_id):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fd = os.open(str(self.path), os.O_CREAT | os.O_RDWR, 0o644)
        try:
            fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except (IOError, OSError) as exc:
            if exc.errno in (errno.EACCES, errno.EAGAIN, errno.EWOULDBLOCK):
                try:
                    holder = os.read(self.fd, 4096).decode("utf-8", "replace")
                except Exception:
                    holder = "<unreadable>"
                os.close(self.fd)
                self.fd = None
                die("G-B07", "another bench2 run holds the host lock %s\n       holder: %s"
                    % (self.path, holder.strip()))
            raise
        token = canon({"pid": os.getpid(), "runset_id": runset_id,
                       "host": platform.node(), "since": utc_now()})
        os.ftruncate(self.fd, 0)
        os.lseek(self.fd, 0, os.SEEK_SET)
        os.write(self.fd, token)
        os.fsync(self.fd)
        # Readback: detects a lock file the filesystem did not actually give us (read-only mount,
        # broken NFS export). It does NOT prove flock exclusion semantics over NFS/Lustre -- that
        # requires the mount to honour flock and is listed in README UNCOVERED.
        os.lseek(self.fd, 0, os.SEEK_SET)
        if os.read(self.fd, len(token)) != token:
            die("G-B07", "host lock %s did not read back what was written; the filesystem is not "
                         "honouring the lock and exclusivity cannot be assumed." % self.path)
        return self

    def release(self):
        if self.fd is not None:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_UN)
            finally:
                os.close(self.fd)
                self.fd = None


# ------------------------------------------------------------------------- streaming validators
MZML_TERMINATORS = (b"</indexedmzML>", b"</mzML>")
_SPECTRUM_LIST_RE = re.compile(br'<spectrumList[^>]*\bcount="(\d+)"')
_SPECTRUM_START_RE = re.compile(br"<spectrum\b")
_SOURCEFILE_RE = re.compile(br'<sourceFile\b[^>]*\blocation="([^"]*)"')


def hash_and_validate_mzml(path):
    """G-B11 -- one streaming pass computes the digest AND proves the file is complete.

    v1 accepted any file that existed and was nonempty, which is exactly how a truncated 12 GB
    mzML from a killed converter got searched and published. Three independent completeness
    signals, all cheap because they ride along with the hash we have to compute anyway:
      * the file ends with a closing mzML tag  (truncation)
      * <spectrumList count="N"> equals the observed number of <spectrum starts  (partial write)
      * the declared sourceFile location is recovered  (corroboration only, see note)

    NOTE on sourceFile: this is a STRING recorded by the converter. It corroborates; it is not
    identity. Identity of the input is raw_content_id (G-B13). v1's mistake was treating a
    substring in metadata as a fingerprint; that mistake is not repeated here.
    """
    h = hashlib.sha256()
    declared = None
    seen = 0
    source_locations = []
    carry = b""
    tail = b""
    size = 0
    with open(str(path), "rb") as fh:
        while True:
            blk = fh.read(1 << 20)
            if not blk:
                break
            size += len(blk)
            h.update(blk)
            window = carry + blk
            if declared is None:
                m = _SPECTRUM_LIST_RE.search(window)
                if m:
                    declared = int(m.group(1))
            seen += len(_SPECTRUM_START_RE.findall(window)) - len(_SPECTRUM_START_RE.findall(carry))
            for m in _SOURCEFILE_RE.finditer(window):
                loc = m.group(1).decode("utf-8", "replace")
                if loc not in source_locations:
                    source_locations.append(loc)
            carry = window[-4096:]
            tail = window[-64:]
    if size == 0:
        die("G-B11", "%s is empty" % path)
    if not any(tail.rstrip().endswith(t) for t in MZML_TERMINATORS):
        die("G-B11", "%s does not end with </mzML> or </indexedmzML>: the file is TRUNCATED.\n"
                     "       A truncated mzML that is merely 'nonempty' is what v1 accepted and "
                     "searched." % path)
    if declared is None:
        die("G-B11", "%s has no <spectrumList count=...>; completeness cannot be established"
            % path)
    if declared != seen:
        die("G-B11", "%s declares %d spectra but contains %d. The converter did not finish "
                     "writing." % (path, declared, seen))
    return {"sha256": h.hexdigest(), "bytes": size,
            "spectrum_count_declared": declared, "spectrum_count_seen": seen,
            "source_file_locations": source_locations}


def hash_and_validate_tsv(path, required_columns):
    """G-B11 (TSV form) -- header contract, field-count contract, trailing newline.

    A Sage run killed at 90%% leaves a TSV whose last row is half a line. v1 counted it. The
    trailing-newline requirement is the cheapest possible truncation detector for a line-oriented
    format, and the per-row field count catches an embedded tab or a partial write mid-file.
    """
    h = hashlib.sha256()
    size = 0
    header = None
    nrows = 0
    bad_row = None
    ends_nl = False
    with open(str(path), "rb") as fh:
        pending = b""
        while True:
            blk = fh.read(1 << 20)
            if not blk:
                break
            size += len(blk)
            h.update(blk)
            ends_nl = blk.endswith(b"\n")
            pending += blk
            *lines, pending = pending.split(b"\n")
            for line in lines:
                if header is None:
                    header = line.decode("utf-8", "replace").rstrip("\r").split("\t")
                    continue
                if not line.strip():
                    continue
                nrows += 1
                if bad_row is None:
                    n = line.decode("utf-8", "replace").rstrip("\r").count("\t") + 1
                    if n != len(header):
                        bad_row = (nrows, n)
        if pending.strip():
            bad_row = bad_row or (nrows + 1, -1)
    if size == 0:
        die("G-B11", "%s is empty" % path)
    if header is None:
        die("G-B11", "%s has no header row" % path)
    if not ends_nl:
        die("G-B11", "%s does not end with a newline: the file is TRUNCATED and its last row is "
                     "a partial record." % path)
    if bad_row is not None:
        die("G-B11", "%s row %d has %d fields, header has %d"
            % (path, bad_row[0], bad_row[1], len(header)))
    missing = [c for c in required_columns if c not in header]
    if missing:
        die("G-B11", "%s is missing required column(s) %s.\n       header seen: %s\n"
                     "       v1 did `row.get('peptide_q', 1)` here, which turned a renamed column "
                     "into a published zero." % (path, missing, header))
    if nrows == 0:
        die("G-B11", "%s has a header but no rows" % path)
    return {"sha256": h.hexdigest(), "bytes": size, "header": header, "rows": nrows,
            "header_sha256": sha256_bytes("\t".join(header).encode("utf-8"))}


# --------------------------------------------------------------------------------- raw identity
def raw_identity(raw_dir):
    """G-B13 -- identity of a Bruker .d from its CONTENT, plus inode-level stat pinning.

    v1 asserted that the string "dataset A-A-3" appeared in the path and in the first 400 KB of the
    produced mzML. A symlink retargeted between preflight and execution defeats both: the
    converter reads dataset D and writes dataset A's path into the header. That check is DELETED here rather
    than patched, because its existence produced false assurance.

    What this returns instead:
      content_id      digest over a deterministic walk of every member file
      acquisition_id  digest over instrument/method facts read out of analysis.tdf
      stat_pins       (dev, ino, mtime_ns, size) per member, re-verified after the run

    COVERAGE LIMIT, stated in the artifact itself: files above FULL_HASH_LIMIT (i.e. .tdf_bin,
    which is the actual spectra) are sampled head+tail, not hashed whole. A modified middle is
    invisible to content_id. This is a deliberate cost tradeoff -- hashing 40 GB per recipe is not
    viable -- and it is why the returned dict carries `coverage`. It is never called "content
    identity" without that qualifier.
    """
    raw = Path(raw_dir)
    if not raw.exists():
        die("G-B13", "raw not found: %s" % raw)
    real = raw.resolve()
    members = []
    pins = []
    for p in sorted(real.rglob("*")):
        if not p.is_file():
            continue
        st = p.stat()
        rel = str(p.relative_to(real))
        pins.append({"p": rel, "dev": st.st_dev, "ino": st.st_ino,
                     "mtime_ns": st.st_mtime_ns, "n": st.st_size})
        if st.st_size <= FULL_HASH_LIMIT:
            members.append({"p": rel, "n": st.st_size, "h": sha256_file(p), "m": "full"})
        else:
            with open(str(p), "rb") as fh:
                head = fh.read(SAMPLE_CHUNK)
                fh.seek(max(0, st.st_size - SAMPLE_CHUNK))
                tailb = fh.read(SAMPLE_CHUNK)
            members.append({"p": rel, "n": st.st_size, "m": "head_tail_1MiB",
                            "h0": sha256_bytes(head), "h1": sha256_bytes(tailb)})
    if not members:
        die("G-B13", "%s contains no files" % raw)
    acq = acquisition_identity(real)
    return {
        "realpath": str(real),
        "content_id": sha256_bytes(canon(members)),
        "coverage": "full<=%dMiB; head+tail %dMiB above" % (FULL_HASH_LIMIT >> 20,
                                                            SAMPLE_CHUNK >> 20),
        "acquisition_id": acq["acquisition_id"],
        "method_id": acq["method_id"],
        # The window table itself travels in the seal so collate2 never has to re-open the .d
        # (and therefore can never fall back to a hardcoded table when it cannot).
        "method": acq["method"],
        "stat_pins": pins,
    }


def _tdf_path(raw_real):
    p = Path(raw_real) / "analysis.tdf"
    if not p.exists():
        die("G-B18", "%s has no analysis.tdf; acquisition method identity cannot be established "
                     "and there is deliberately no fallback window table." % raw_real)
    return p


def acquisition_identity(raw_real):
    """G-B18 -- per-sample acquisition method, read from that sample's own analysis.tdf.

    v1 hardcoded WINDOWS_S23 and applied it to dataset A and dataset D. Nothing recorded that those samples
    shared an acquisition method; nothing read DiaFrameMsMsWindows from their raw data. Spectra
    outside dataset B's tiles were silently skipped, which shows up as depressed recall rather than as
    an error. WINDOWS_S23 is deleted. There is no default table.
    """
    tdf = _tdf_path(raw_real)
    con = sqlite3.connect("file:%s?mode=ro" % tdf, uri=True)
    try:
        meta = {}
        try:
            for k, v in con.execute("SELECT Key, Value FROM GlobalMetadata"):
                meta[str(k)] = str(v)
        except sqlite3.Error as exc:
            die("G-B18", "%s: cannot read GlobalMetadata (%s)" % (tdf, exc))
        try:
            rows = list(con.execute(
                "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth "
                "FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin"))
        except sqlite3.Error as exc:
            die("G-B18", "%s: no readable DiaFrameMsMsWindows table (%s). This is not diaPASEF "
                         "data, or the .d is incomplete." % (tdf, exc))
        if not rows:
            die("G-B18", "%s: DiaFrameMsMsWindows is empty" % tdf)
        try:
            nframes = list(con.execute("SELECT COUNT(*) FROM Frames"))[0][0]
            maxscans = list(con.execute("SELECT MAX(NumScans) FROM Frames"))[0][0]
        except sqlite3.Error as exc:
            die("G-B18", "%s: cannot read Frames (%s)" % (tdf, exc))
    finally:
        con.close()

    for key in ("OneOverK0AcqRangeLower", "OneOverK0AcqRangeUpper"):
        if key not in meta:
            die("G-B18", "%s: GlobalMetadata lacks %s; ion-mobility bounds for the acquisition "
                         "tiles cannot be derived and 2-D window assignment would silently "
                         "degrade to 1-D." % (tdf, key))
    k0_lo = float(meta["OneOverK0AcqRangeLower"])
    k0_hi = float(meta["OneOverK0AcqRangeUpper"])
    if not maxscans or int(maxscans) < 2:
        die("G-B18", "%s: Frames.NumScans is %r; cannot map scan number to 1/K0"
            % (tdf, maxscans))
    nscans = int(maxscans)

    def k0_of_scan(scan):
        # Scan 0 is the HIGH 1/K0 end of the ramp. This is a LINEAR INTERPOLATION, not the Bruker
        # SDK conversion; the model name travels in method_id and in estimand_id so a number
        # computed under it can never share a column with one computed under a different model.
        # Its approximation error is documented in README UNCOVERED.
        frac = float(scan) / float(nscans - 1)
        frac = min(max(frac, 0.0), 1.0)
        return k0_hi - frac * (k0_hi - k0_lo)

    windows = []
    for wg, sb, se, imz, iw in rows:
        imz = float(imz)
        iw = float(iw)
        sb = int(sb)
        se = int(se)
        windows.append({
            "window_group": int(wg),
            "mz_lo": round(imz - iw / 2.0, 6), "mz_hi": round(imz + iw / 2.0, 6),
            "scan_begin": sb, "scan_end": se,
            "im_lo": round(min(k0_of_scan(sb), k0_of_scan(se)), 6),
            "im_hi": round(max(k0_of_scan(sb), k0_of_scan(se)), 6),
        })
    method = {"windows": windows, "im_model": "linear_scan_interp/v1",
              "k0_range": [k0_lo, k0_hi], "num_scans": nscans}
    acq_fields = {
        "instrument_serial": meta.get("InstrumentSerialNumber", ""),
        "acquisition_datetime": meta.get("AcquisitionDateTime", ""),
        "sample_name": meta.get("SampleName", ""),
        "frame_count": int(nframes),
        "method": method,
    }
    return {"acquisition_id": sha256_bytes(canon(acq_fields)),
            "method_id": sha256_bytes(canon(method))[:16],
            "method": method,
            "global_metadata": meta}


def verify_raw_unchanged(raw_dir, pinned):
    """G-B13 (post-run half) -- the .d must be the same bytes AFTER the 6 h run as before.

    The pre-run digest alone is a TOCTOU hole: the recipe pins dataset A's digest at construction and
    the converter opens the .d minutes later and reads it for hours. Retarget the symlink in that
    window and the recipe says dataset A while the tool read dataset D.

    This does NOT prevent the retarget. It DETECTS it, at seal time, and the attempt aborts -- so
    nothing is sealed and nothing can be published. A retarget-and-retarget-back inside the window
    still evades, and that residual is in README UNCOVERED.
    """
    real = Path(raw_dir).resolve()
    if str(real) != pinned["realpath"]:
        die("G-B13", "raw path %s now resolves to %s (was %s). The symlink was RETARGETED during "
                     "the run." % (raw_dir, real, pinned["realpath"]))
    now = {}
    for p in sorted(real.rglob("*")):
        if p.is_file():
            st = p.stat()
            now[str(p.relative_to(real))] = (st.st_dev, st.st_ino, st.st_mtime_ns, st.st_size)
    before = dict((m["p"], (m["dev"], m["ino"], m["mtime_ns"], m["n"]))
                  for m in pinned["stat_pins"])
    if set(now) != set(before):
        die("G-B13", "the raw .d gained or lost files during the run: added=%s removed=%s"
            % (sorted(set(now) - set(before))[:5], sorted(set(before) - set(now))[:5]))
    for rel in sorted(before):
        if now[rel] != before[rel]:
            die("G-B13", "raw member %s changed during the run (dev/ino/mtime/size %s -> %s). The "
                         "input was not stable for the duration of the measurement."
                % (rel, before[rel], now[rel]))
    acq = acquisition_identity(real)
    if acq["acquisition_id"] != pinned["acquisition_id"]:
        die("G-B13", "analysis.tdf acquisition identity changed during the run (%s -> %s)"
            % (short(pinned["acquisition_id"]), short(acq["acquisition_id"])))


# ------------------------------------------------------------------------------------- registry
def load_registry():
    try:
        import yaml
    except ImportError:
        die("G-B00", "need pyyaml: python3 -m pip install --user pyyaml")
    # $SPEXTRACTOR_REGISTRY exists so --selftest can stand up a synthetic registry. It is NOT a
    # general escape hatch: the runset records which registry projection every cell used.
    path = Path(os.environ.get("SPEXTRACTOR_REGISTRY") or (HERE / "samples.yaml"))
    if not path.exists():
        die("G-B00", "registry not found: %s" % path)
    reg = yaml.safe_load(path.read_text())
    if not isinstance(reg, dict):
        die("G-B00", "%s is not a mapping" % path)
    for key in ("samples", "search", "vault"):
        if key not in reg:
            die("G-B00", "%s lacks required top-level key '%s'" % (path, key))
    return reg, path


def reg_get(mapping, key, where):
    """Typed registry access. There is no .get() with a falsy default anywhere on this path.

    v1 did `Path(entry["references"].get("diann_lib", ""))`, and Path("") is ".", which exists,
    so a missing library became "read the current directory as parquet" and, downstream, a
    legitimate-looking 0%% purity.
    """
    if not isinstance(mapping, dict) or key not in mapping or mapping[key] in (None, ""):
        die("G-B06", "registry entry %s lacks required key '%s'" % (where, key))
    return mapping[key]


def sample_entry(reg, sample):
    samples = reg["samples"]
    if sample not in samples:
        die("G-B05", "sample '%s' is not in samples.yaml. A sample is an ID resolved from the "
                     "registry; bare paths are not accepted anywhere." % sample)
    entry = samples[sample]
    if not isinstance(entry, dict):
        die("G-B05", "registry entry for %s is malformed" % sample)
    return entry


def require_pinned(reg, sample):
    """G-B05 -- a sample must carry pinned content digests before it can be benchmarked.

    Pinning converts every later use from an existence PRINT (v1: "[preflight] dataset A.diann_report:
    OK") into a digest COMPARISON. Note the honest limit: pinning makes a reference STABLE, not
    CORRECT. If dataset A's diann_report already points at an dataset D parquet, pinning freezes that error.
    The semantic check that actually catches it lives in collate2 (G-C09: the DIA-NN Run column
    must name this sample's raw file).
    """
    entry = sample_entry(reg, sample)
    pinned = entry.get("pinned")
    if not pinned:
        die("G-B05", "sample '%s' is not pinned. Run:\n         bench2.py pin --sample %s\n"
                     "       Until then its raw identity is a path, and a path is not an identity."
            % (sample, sample))
    for key in ("raw_content_id", "acquisition_id", "method_id", "pinned_utc"):
        if key not in pinned:
            die("G-B05", "sample '%s' pin lacks '%s'; re-pin it" % (sample, key))
    return entry, pinned


def verify_pinned_references(reg, sample):
    """G-B06 -- declared references must exist AND hash to their pinned digest.

    v1 printed MISSING and continued. The claim that "absence is a fact about the registry"
    was false: nothing aborted, and truth silently skipped the sample.
    """
    entry = sample_entry(reg, sample)
    refs = entry.get("references")
    if not refs:
        die("G-B06", "sample '%s' declares no references" % sample)
    out = {}
    for name in sorted(refs):
        rec = refs[name]
        if not isinstance(rec, dict):
            die("G-B06", "reference %s.%s must be a mapping with path+sha256, got %r.\n"
                         "       A bare path is an existence claim, not an identity."
                % (sample, name, rec))
        path = Path(reg_get(rec, "path", "%s.references.%s" % (sample, name)))
        want = reg_get(rec, "sha256", "%s.references.%s" % (sample, name))
        if not path.exists():
            die("G-B06", "reference %s.%s missing: %s" % (sample, name, path))
        if path.stat().st_size == 0:
            die("G-B06", "reference %s.%s is empty: %s" % (sample, name, path))
        got = sha256_file(path)
        if got != want:
            die("G-B06", "reference %s.%s changed since it was pinned\n       path: %s\n"
                         "       pinned: %s\n       now:    %s"
                % (sample, name, path, short(want), short(got)))
        out[name] = {"path": str(path), "sha256": got, "bytes": path.stat().st_size}
    return out


def verify_search_toolchain(reg):
    """G-B06 (search half) -- Sage binary, config and FASTA are hashed.

    v1 hashed none of them. The path merely contained the string "v0.14.7", and reg['search']
    ['fasta'] was declared and never used by anything.
    """
    s = reg["search"]
    out = {}
    for name in ("sage_binary", "sage_config", "fasta"):
        rec = s.get(name)
        if not isinstance(rec, dict):
            die("G-B06", "search.%s must be a mapping with path+sha256 (got %r)" % (name, rec))
        path = Path(reg_get(rec, "path", "search.%s" % name))
        want = reg_get(rec, "sha256", "search.%s" % name)
        if not path.exists():
            die("G-B06", "search.%s missing: %s" % (name, path))
        got = sha256_file(path)
        if got != want:
            die("G-B06", "search.%s changed since pinning: %s -> %s"
                % (name, short(want), short(got)))
        out[name] = {"path": str(path), "sha256": got}
    # The Sage config names a FASTA. Hashing the config detects drift but not that the config
    # points somewhere other than the registry FASTA -- so parse it and compare.
    try:
        cfg = jload(Path(out["sage_config"]["path"]))
    except Abort:
        raise
    except Exception as exc:
        die("G-B06", "search.sage_config is not readable JSON: %s" % exc)
    named = None
    db = cfg.get("database") if isinstance(cfg, dict) else None
    if isinstance(db, dict):
        named = db.get("fasta")
    if named is None:
        named = cfg.get("fasta") if isinstance(cfg, dict) else None
    if named is not None:
        if os.path.realpath(str(named)) != os.path.realpath(out["fasta"]["path"]):
            die("G-B06", "sage_config names FASTA %s but the registry pins %s. The search would "
                         "not be against the database this runset claims."
                % (named, out["fasta"]["path"]))
    out["sage_config_names_fasta"] = named
    return out


# ----------------------------------------------------------------------------------------- plan
def _reject_unknown(got, allowed, where):
    """G-B01 -- an unknown key in the plan ABORTS.

    This is the typo class: `basline: base` parses fine, is silently ignored by `.get('baseline')`,
    and the arm publishes as if it were unpaired-by-design. A whitelist turns that into an abort
    at plan-check time, before any compute.
    """
    unknown = sorted(set(got) - set(allowed))
    if unknown:
        die("G-B01", "unknown key(s) %s in %s. Allowed: %s\n"
                     "       A typo here (e.g. 'basline') silently unpairs an arm in v1."
            % (unknown, where, sorted(allowed)))


def load_plan(path, reg):
    try:
        import yaml
    except ImportError:
        die("G-B00", "need pyyaml")
    p = Path(path)
    if not p.exists():
        die("G-B01", "plan not found: %s" % p)
    raw_bytes = p.read_bytes()
    raw = yaml.safe_load(raw_bytes.decode("utf-8"))
    if not isinstance(raw, dict):
        die("G-B01", "plan is not a mapping")
    _reject_unknown(raw.keys(), PLAN_TOP_KEYS, "plan top level")
    for key in ("schema", "title", "binary", "samples", "arms"):
        if key not in raw:
            die("G-B01", "plan lacks required key '%s'" % key)
    if raw["schema"] != "spextractor.plan/2":
        die("G-B01", "plan schema is %r; this harness reads only 'spextractor.plan/2'. There is no "
                     "forward-compatibility guessing." % raw["schema"])

    samples = raw["samples"]
    if not isinstance(samples, list) or not samples:
        die("G-B01", "plan.samples must be a non-empty list")
    for s in samples:
        if not SAMPLE_NAME_RE.match(str(s)):
            die("G-B02", "illegal sample name %r" % s)
    if len(set(samples)) != len(samples):
        die("G-B02", "duplicate sample in plan.samples: %s" % samples)

    arms = []
    seen_names = set()
    for a in raw["arms"]:
        if not isinstance(a, dict):
            die("G-B01", "each arm must be a mapping")
        _reject_unknown(a.keys(), PLAN_ARM_KEYS, "arm %r" % a.get("name"))
        name = a.get("name")
        # G-B02: arm names are validated because in v1 they were PATH COMPONENTS. A name
        # containing '/' or '..' escaped the output tree; a duplicate name silently overwrote
        # the previous arm's directory. Here names are never path components at all (stage
        # directories are named by digest), but a duplicate would still collide in the cell map.
        if not isinstance(name, str) or not ARM_NAME_RE.match(name):
            die("G-B02", "illegal arm name %r; must match %s" % (name, ARM_NAME_RE.pattern))
        if name in seen_names:
            die("G-B02", "duplicate arm name %r. In v1 this silently overwrote a directory and "
                         "one arm's numbers were published as the other's." % name)
        seen_names.add(name)
        role = a.get("role")
        if role not in PLAN_ROLES:
            die("G-B01", "arm %r has role %r; must be one of %s" % (name, role, sorted(PLAN_ROLES)))
        kind = a.get("kind", "run")
        if kind not in PLAN_KINDS:
            die("G-B01", "arm %r has kind %r" % (name, kind))
        params = a.get("params") or {}
        if not isinstance(params, dict):
            die("G-B01", "arm %r params must be a mapping" % name)
        if kind == "reference_mzml":
            if not a.get("source"):
                die("G-B01", "arm %r is kind reference_mzml and must declare `source` "
                             "(e.g. references.reference_mzml)" % name)
            if params:
                die("G-B01", "arm %r is an external mzML; it has no converter and therefore no "
                             "parameters" % name)
        arms.append({"name": name, "role": role, "kind": kind,
                     "baseline": a.get("baseline"), "params": dict(params),
                     "source": a.get("source")})

    arms.extend(_expand_sweeps(raw.get("sweeps") or [], seen_names))
    _validate_pairing(arms)
    _validate_scaled_params(reg, arms, raw.get("sweeps") or [])

    plan = {
        "schema": raw["schema"],
        "title": raw["title"],
        "binary": str(raw["binary"]),
        "threads": int(raw.get("threads", 120)),
        "tcmalloc": bool(raw.get("tcmalloc", True)),
        "samples": list(samples),
        "arms": arms,
    }
    normalized = {
        "schema": plan["schema"], "binary": plan["binary"], "threads": plan["threads"],
        "tcmalloc": plan["tcmalloc"], "samples": sorted(plan["samples"]),
        "arms": sorted([{k: v for k, v in arm.items()} for arm in arms],
                       key=lambda x: x["name"]),
    }
    plan["file_sha256"] = sha256_bytes(raw_bytes)
    plan["normalized_sha256"] = sha256_bytes(canon(normalized))
    plan["source_bytes"] = raw_bytes
    plan["path"] = str(p.resolve())
    return plan


def _expand_sweeps(sweeps, seen_names):
    """G-B04 (generator half) -- a sweep is the ONLY way a scaled parameter can enter a plan."""
    out = []
    for sw in sweeps:
        if not isinstance(sw, dict):
            die("G-B01", "each sweep must be a mapping")
        _reject_unknown(sw.keys(), PLAN_SWEEP_KEYS, "sweep %r" % sw.get("param"))
        for key in ("param", "baseline", "points"):
            if key not in sw:
                die("G-B01", "sweep lacks required key '%s'" % key)
        points = sw["points"]
        if not isinstance(points, list):
            die("G-B01", "sweep %r points must be a list" % sw["param"])
        # FAILURE 6: ambiguity_margin 0.1 was compared at a SINGLE point against an integer
        # partner count. One point cannot distinguish "this parameter does nothing" from "this
        # parameter is on the wrong scale". Three is the minimum that shows a trend.
        if len(points) < 3:
            die("G-B04", "sweep of %r has %d point(s); at least 3 are required.\n"
                         "       This is FAILURE 6: a single-point comparison of a rescaled "
                         "parameter cannot distinguish a null effect from a scale error."
                % (sw["param"], len(points)))
        vals = [float(x) for x in points]
        if sorted(vals) != vals and sorted(vals, reverse=True) != vals:
            die("G-B04", "sweep of %r points are not monotone: %s" % (sw["param"], points))
        if len(set(vals)) != len(vals):
            die("G-B04", "sweep of %r has duplicate points" % sw["param"])
        slug = re.sub(r"[^A-Za-z0-9]+", "", sw["param"].split(":")[-1])[:20] or "sweep"
        for pt in points:
            name = ("%s_%s" % (slug, str(pt).replace(".", "p").replace("-", "m")))[:64]
            if name in seen_names:
                die("G-B02", "sweep generated duplicate arm name %r" % name)
            seen_names.add(name)
            out.append({"name": name, "role": "treatment", "kind": "run",
                        "baseline": sw["baseline"], "params": {sw["param"]: pt}, "source": None})
    return out


def _validate_pairing(arms):
    """G-B03 -- every treatment resolves to a baseline IN THIS PLAN.

    v1's preflight did `if not arm.baseline: continue`, which skipped ALL validation for an arm
    that simply omitted the field -- and then ran it and published it anyway. That is FAILURE 2:
    "+36.7%/+37.4% confirmed" used three different baseline configs.
    """
    by_name = dict((a["name"], a) for a in arms)
    for a in arms:
        if a["role"] == "baseline":
            if a["baseline"]:
                die("G-B03", "arm %r has role baseline and must not declare a baseline"
                    % a["name"])
            continue
        if a["role"] == "external_baseline":
            continue
        if not a["baseline"]:
            pass  # MUTANT: G-B03 disabled
        if a["baseline"] not in by_name:
            die("G-B03", "arm %r names baseline %r which is not in this plan.\n"
                         "       A pair is a pair of cells in ONE runset; cross-plan pairing is "
                         "not representable." % (a["name"], a["baseline"]))
        parent = by_name[a["baseline"]]
        if parent["role"] != "baseline":
            die("G-B03", "arm %r names %r as baseline but that arm has role %r; baseline chains "
                         "deeper than 1 are refused" % (a["name"], parent["name"], parent["role"]))
    if not any(a["role"] == "baseline" for a in arms):
        die("G-B03", "plan declares no baseline arm")


def _validate_scaled_params(reg, arms, sweeps):
    """G-B04 (enforcement half) -- a scaled parameter may ONLY appear via a sweep.

    The registry, not the plan, declares which parameters have a scale that has bitten before.
    A plan cannot opt out by hand-writing the arm, which is exactly the bypass a plan-level
    declaration would leave open.
    """
    scaled = set((reg.get("tool") or {}).get("scaled_params") or [])
    if not scaled:
        return
    swept = set(sw["param"] for sw in sweeps)
    for a in arms:
        for k in a["params"]:
            if k in scaled and k not in swept:
                die("G-B04", "arm %r sets scaled parameter %r directly.\n"
                             "       Registry `tool.scaled_params` lists it because its units are "
                             "not self-evident (FAILURE 6: a 0.1 FRACTION was compared against an "
                             "integer COUNT). It may only enter through a >=3-point sweep."
                    % (a["name"], k))


# ------------------------------------------------------------------------------------ tool + ini
def canonical_ini(path):
    """G-B14 (config half) -- typed extraction from the tool's OWN resolved INI.

    v1 parsed `--helphelp` with a regex and diffed strings. That misses equivalent spellings
    (0.1 vs 0.10), boolean aliases (true vs 1), aliased flags, and any parameter the tool
    normalizes -- and a help-format change silently degrades it.

    Here the tool writes its own resolved INI, we extract <ITEM name/type/value> as a typed map,
    and hash that. Never a regex over the XML: a text canonicalizer that under-strips makes two
    identical configs hash differently, which would make a real no-op look like a real difference
    -- failing in the dangerous direction.
    """
    try:
        tree = ET.parse(str(path))
    except ET.ParseError as exc:
        die("G-B14", "resolved INI %s is not parseable XML: %s" % (path, exc))
    items = {}

    def walk(node, prefix):
        for child in node:
            tag = child.tag.split("}")[-1]
            name = child.get("name", "")
            if tag == "ITEM":
                key = "%s:%s" % (prefix, name) if prefix else name
                items[key] = [child.get("type", ""), child.get("value", "")]
            elif tag == "ITEMLIST":
                key = "%s:%s" % (prefix, name) if prefix else name
                items[key] = ["itemlist", [g.get("value", "") for g in child]]
            elif tag in ("NODE", "PARAMETERS", "ROOT"):
                walk(child, "%s:%s" % (prefix, name) if (prefix and name) else (name or prefix))

    walk(tree.getroot(), "")
    if not items:
        die("G-B14", "resolved INI %s contained no <ITEM> elements; no-op detection would be "
                     "disabled and this harness does not run with a disabled guard." % path)
    return items, sha256_bytes(canon(items))


def write_resolved_ini(binary, params, threads, dest):
    """Ask the tool for its own resolved configuration, then apply the arm's parameters."""
    proc = subprocess.run([str(binary), "-write_ini", str(dest)],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
    if proc.returncode != 0 or not dest.exists():
        die("G-B14", "%s -write_ini failed (exit %s):\n%s"
            % (binary, proc.returncode, proc.stdout.decode("utf-8", "replace")[-2000:]))
    tree = ET.parse(str(dest))
    root = tree.getroot()
    index = {}

    def walk(node, prefix):
        for child in node:
            tag = child.tag.split("}")[-1]
            name = child.get("name", "")
            if tag == "ITEM":
                key = "%s:%s" % (prefix, name) if prefix else name
                index[key] = child
                index[name] = index.get(name, child)
            elif tag in ("NODE", "PARAMETERS", "ROOT", "ITEMLIST"):
                walk(child, "%s:%s" % (prefix, name) if (prefix and name) else (name or prefix))

    walk(root, "")
    applied = {}
    for key, value in sorted(params.items()):
        node = index.get(key) or index.get(key.split(":")[-1])
        if node is None:
            # v1 probed parameters by looking for the literal phrase "Unknown option" in output
            # and ignoring the return code. Absence from the tool's own resolved INI is the
            # authoritative test.
            die("G-B16", "parameter %r does not exist in this binary's resolved configuration.\n"
                         "       Known keys include: %s"
                % (key, ", ".join(sorted(k for k in index if ":" in k)[:12])))
        node.set("value", str(value))
        applied[key] = str(value)
    tnode = index.get("threads") or index.get("1:threads")
    if tnode is not None:
        tnode.set("value", str(threads))
    tree.write(str(dest), encoding="UTF-8", xml_declaration=True)
    return applied


def validate_converter_argv(argv):
    """G-B16 -- only -ini/-in/-out may appear.

    If a parameter could ride on the command line it would override the INI and be invisible to
    ini_sha256, so no-op detection (FAILURE 5) would be aspirational rather than exact.
    """
    flags = [tok for tok in argv[1:] if tok.startswith("-")]
    bad = [f for f in flags if f not in ALLOWED_CONVERTER_ARGV_FLAGS]
    if bad:
        die("G-B16", "converter argv carries disallowed flag(s) %s. Every output-affecting "
                     "parameter must travel through the hashed INI, or no-op detection is blind "
                     "to it." % bad)
    return argv


def tool_identity(binary):
    p = Path(binary)
    if not p.exists() or not p.is_file():
        die("G-B09", "binary not found: %s" % binary)
    real = p.resolve()
    st = real.stat()
    if not os.access(str(real), os.X_OK):
        die("G-B09", "binary is not executable: %s" % real)
    return {"realpath": str(real), "sha256": sha256_file(real), "bytes": st.st_size,
            "mtime_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)),
            "dev": st.st_dev, "ino": st.st_ino}


def assert_binary_unchanged(binary, before):
    """G-B09 -- the binary is hashed immediately before AND after EVERY exec.

    v1 hashed once at preflight and stamped that cached value into every manifest. Replace the
    binary between arm 1 and arm 2 and both manifests carry binary A's hash while two different
    binaries executed -- and collate's "one binary across the set" check passes trivially. That
    is FAILURE 7, and the v1 comment claimed identity was "asserted equal" when no assertion
    existed anywhere in the file.
    """
    after = tool_identity(binary)
    for key in ("realpath", "sha256", "ino"):
        if after[key] != before[key]:
            die("G-B09", "the tool changed across this exec: %s %s -> %s\n"
                         "       tar silently skips overwriting a RUNNING executable, which is how "
                         "a 6 h benchmark ran a 2-day-old build (FAILURE 7)."
                % (key, before[key], after[key]))
    return after


# ------------------------------------------------------------------------------- commit protocol
class Bench(object):
    """The one global content-addressed store.

    Not per-plan: two plans with the same filename stem shared a scratch root and a vault filename
    in v1, so one silently clobbered the other. Nothing here is addressed by a plan's filename.
    """

    def __init__(self, root):
        self.root = Path(root)
        for sub in ("staging", "failed", "stages", "runsets", "locks"):
            (self.root / sub).mkdir(parents=True, exist_ok=True)
        for kind in ("convert", "adopt", "search", "metrics"):
            (self.root / "stages" / kind).mkdir(parents=True, exist_ok=True)
        assert_same_filesystem([self.root / "staging", self.root / "stages",
                                self.root / "failed"])

    def stage_dir(self, kind, stage_id):
        return self.root / "stages" / kind / stage_id

    def load_stage(self, kind, stage_id):
        """The ONLY way to reach a sealed stage. There is no glob anywhere in this harness.

        v1's cmd_truth and cmd_collate did root.glob("*__*") under a directory derived from the
        PLAN FILENAME, so they collated whatever happened to be present -- including a copied
        directory carrying a foreign manifest, and including nothing at all (an empty root
        produced an apparently valid empty vault report).
        """
        d = self.stage_dir(kind, stage_id)
        if not d.is_dir():
            die("G-B15", "stage %s/%s does not exist. Stages are named by digest and are reached "
                         "only through the runset contract." % (kind, stage_id))
        seal = jload(d / "SEAL.json")
        verify_seal(seal, d)
        return d, seal


def seal_digest(seal):
    body = dict((k, v) for k, v in seal.items() if k != "seal")
    return sha256_bytes(canon(body))


def verify_seal(seal, where):
    """G-B12 -- the seal detects accidental edits of the sealed metadata.

    HONEST: this is ACCIDENT-EVIDENCE, not tamper-evidence. The digest is unkeyed and the
    function that computes it ships in this repository, so anyone who can write the directory can
    recompute it. It detects drift, a half-written file, and a hand-edit; it does not resist a
    determined operator. Saying "bypassable only by forging SHA-256" would be exactly the kind of
    overclaim that got v1's guards called theatre.
    """
    if "seal" not in seal:
        die("G-B12", "%s has no seal" % where)
    if seal_digest(seal) != seal["seal"]:
        die("G-B12", "%s: SEAL.json does not match its own seal -- it was edited after commit."
            % where)
    if seal.get("status") != "ok":
        die("G-B12", "%s: seal status is %r; only sealed-ok stages are readable"
            % (where, seal.get("status")))
    if seal.get("seal_version") != CONTRACT_VERSION:
        die("G-B12", "%s: seal_version %r, expected %d" % (where, seal.get("seal_version"),
                                                           CONTRACT_VERSION))


class Attempt(object):
    """One attempt to perform one computation.

    The attempt directory is created fresh with a uuid4 nonce in its name. No other process knows
    that name; nothing was ever written there before. THAT is why a nonzero converter exit cannot
    leave a stale artifact behind to be accepted: there is no previous artifact in a directory
    that did not exist a second ago.
    """

    def __init__(self, bench, kind, recipe):
        self.bench = bench
        self.kind = kind
        self.attempt_id = new_attempt_id()
        self.dir = bench.root / "staging" / self.attempt_id
        if self.dir.exists():
            die("G-B12", "attempt directory %s already exists (impossible)" % self.dir)
        self.dir.mkdir(parents=True)
        self.work = self.dir / "work"
        self.work.mkdir()
        (self.dir / "recipe.json").write_bytes(canon(recipe))
        self.recipe = recipe
        self.recipe_id = sha256_bytes(canon(recipe))
        self.stage_id = "%s-%s-%s" % (kind, self.recipe_id[:16], self.attempt_id)
        self.started = utc_now()
        self.t0 = time.time()
        self.exit_codes = {}
        self.artifacts = {}
        self.checks = {}

    def fail(self, reason, detail):
        """Failed attempts keep their bytes for forensics and are NEVER on any input path."""
        dest = self.bench.root / "failed" / self.attempt_id
        info = {"attempt_id": self.attempt_id, "kind": self.kind, "recipe_id": self.recipe_id,
                "reason": reason, "detail": str(detail), "exit_codes": self.exit_codes,
                "started_utc": self.started, "failed_utc": utc_now()}
        try:
            (self.dir / "FAILURE.json").write_bytes(canon(info))
            os.rename(str(self.dir), str(dest))
        except OSError:
            dest = self.dir
        return dest

    def run_process(self, argv, env, logname, timeout=None):
        """G-B10 -- a nonzero return code aborts. There is no 'but the file looks fine' branch.

        v1: `if not mzml.exists() or size == 0: die(...)`. A converter that exited nonzero but
        left the PREVIOUS run's nonempty mzML in place sailed through, the manifest was
        overwritten with fresh metadata, and the stale artifact was searched and published.
        Sage's return code was discarded outright.
        """
        logp = self.dir / logname
        with open(str(logp), "wb") as log:
            proc = subprocess.run(argv, cwd=str(self.work), stdout=log,
                                  stderr=subprocess.STDOUT, env=env, timeout=timeout)
        self.exit_codes[Path(argv[0]).name] = proc.returncode
        if proc.returncode != 0:
            tail = logp.read_bytes()[-4000:].decode("utf-8", "replace")
            dest = self.fail("nonzero_exit", proc.returncode)
            die("G-B10", "%s exited %d. NOTHING is sealed; the attempt was moved to %s.\n"
                         "       v1 accepted a nonzero exit whenever a file from a PREVIOUS "
                         "attempt happened to be present.\n       ---- log tail ----\n%s"
                % (Path(argv[0]).name, proc.returncode, dest, tail))
        return proc

    def declare(self, name, info, role="primary"):
        rec = {"role": role}
        rec.update(info)
        self.artifacts[name] = rec

    def hash_incidentals(self):
        """Every file in the sealed directory is hashed, not only the declared outputs.

        Sage writes results.sage.pin, results.json and friends. If only declared outputs were
        hashed those could be altered undetectably and later read by a human or a future stage.
        """
        for p in sorted(self.dir.rglob("*")):
            if not p.is_file():
                continue
            rel = str(p.relative_to(self.dir))
            if rel in ("recipe.json", "SEAL.json") or rel in self.artifacts:
                continue
            key = rel
            if key in self.artifacts:
                continue
            self.artifacts[key] = {"role": "incidental", "sha256": sha256_file(p),
                                   "bytes": p.stat().st_size}

    def commit(self, extra=None):
        self.hash_incidentals()
        seal = {
            "seal_version": CONTRACT_VERSION,
            "kind": self.kind,
            "stage_id": self.stage_id,
            "attempt_id": self.attempt_id,
            "recipe_id": self.recipe_id,
            "recipe_sha256": sha256_file(self.dir / "recipe.json"),
            "harness_version": HARNESS_VERSION,
            "harness_code_digest": harness_code_digest(),
            "status": "ok",
            "exit_codes": self.exit_codes,
            "artifacts": self.artifacts,
            "checks": self.checks,
            "started_utc": self.started,
            "finished_utc": utc_now(),
            "wall_s": round(time.time() - self.t0, 3),
            "host": host_facts(),
        }
        if extra:
            seal.update(extra)
        seal["seal"] = seal_digest(seal)
        (self.dir / "SEAL.json").write_bytes(canon(seal))
        commit_stage(self.bench, self.kind, self.stage_id, self.dir)
        return self.stage_id, seal


def commit_stage(bench, kind, stage_id, attempt_dir):
    """The atomic commit. This is the whole design in twelve lines.

    Death BEFORE the rename  -> the directory stays in staging/, which no reader enumerates.
    Death DURING the rename  -> rename is atomic; either it happened or it did not.
    Death AFTER the rename   -> contents were fsynced first, so a node crash cannot expose a
                                sealed directory with unflushed bytes. (v1 had no fsync at all.)
    """
    final = bench.stage_dir(kind, stage_id)
    if final.exists():
        die("G-B12", "stage_id collision at %s -- refusing. stage_ids embed a uuid4 nonce, so "
                     "this should be unreachable." % final)
    for p in sorted(attempt_dir.rglob("*")):
        if p.is_file():
            fd = os.open(str(p), os.O_RDONLY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
            os.chmod(str(p), 0o444)
    fsync_dir(attempt_dir)
    final.parent.mkdir(parents=True, exist_ok=True)
    # Rename BEFORE chmod: renaming a directory rewrites its ".." entry and therefore needs write
    # permission on the directory being moved. chmod 0555 first would make the commit impossible.
    os.rename(str(attempt_dir), str(final))
    for p in sorted(final.rglob("*")):
        if p.is_dir():
            os.chmod(str(p), 0o555)
    os.chmod(str(final), 0o555)
    fsync_dir(final.parent)
    return final


def host_facts():
    """RECORDED, not prevented.

    Third-party CPU and memory contention on a shared node cannot be excluded from userspace.
    What is possible is to record enough that a suspicious result can be interrogated after the
    fact. The docstring says RECORDED because that is what it is.
    """
    facts = {"node": platform.node(), "system": platform.system(),
             "release": platform.release(), "python": platform.python_version()}
    try:
        facts["loadavg"] = [round(x, 2) for x in os.getloadavg()]
    except (OSError, AttributeError):
        facts["loadavg"] = None
    try:
        facts["ncpu"] = os.cpu_count()
    except Exception:
        facts["ncpu"] = None
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    facts["mem_available_kb"] = int(line.split()[1])
                    break
    except (IOError, OSError):
        facts["mem_available_kb"] = None
    return facts


# ------------------------------------------------------------------------------------- runsets
def runset_path(bench, runset_id, revision):
    return bench.root / "runsets" / ("%s.r%d.json" % (runset_id, revision))


def write_runset(bench, obj):
    """Runsets are APPEND-ONLY revisions; nothing is ever mutated in place.

    v1 overwrote manifest.json on every attempt, so a rerun replaced the provenance of the run
    whose artifacts were still on disk. Here .r1 is written by `run` and .r2 by `truth`; the
    manifest-overwrite failure has no code path.
    """
    obj["seal"] = seal_digest(obj)
    p = runset_path(bench, obj["runset_id"], obj["revision"])
    if p.exists():
        die("G-B15", "%s already exists; runset revisions are append-only" % p)
    tmp = p.with_suffix(".tmp%d" % os.getpid())
    tmp.write_bytes(canon(obj))
    os.rename(str(tmp), str(p))
    os.chmod(str(p), 0o444)
    fsync_dir(p.parent)
    return p


def load_runset(bench, ref):
    """Resolve `<runset_id>` (latest revision) or `<runset_id>.rN`."""
    if ".r" in ref:
        rid, rev = ref.rsplit(".r", 1)
        p = runset_path(bench, rid, int(rev))
        if not p.exists():
            die("G-B15", "no such runset revision: %s" % p)
    else:
        cands = sorted((bench.root / "runsets").glob("%s.r*.json" % ref))
        if not cands:
            die("G-B15", "no runset %r under %s" % (ref, bench.root / "runsets"))
        p = max(cands, key=lambda x: int(x.name.rsplit(".r", 1)[1].split(".")[0]))
    obj = jload(p)
    if obj.get("runset_version") != CONTRACT_VERSION:
        die("G-B15", "%s: runset_version %r" % (p, obj.get("runset_version")))
    if seal_digest(obj) != obj.get("seal"):
        die("G-B15", "%s was modified after it was written. The runset contract is the only "
                     "thing that says which runs belong to this plan; a mutated contract is not "
                     "readable." % p)
    return obj, p


def ledger_append(bench, record):
    """Pre-registration: a runset is recorded BEFORE its results exist.

    Without this, nothing stops running a plan three times and publishing the best one. The
    ledger makes the count of attempts for a plan digest visible in the published report.
    """
    p = bench.root / "runsets" / "LEDGER.jsonl"
    with open(str(p), "ab") as fh:
        fh.write(canon(record) + b"\n")
        fh.flush()
        os.fsync(fh.fileno())


def ledger_count(bench, plan_normalized_sha256):
    p = bench.root / "runsets" / "LEDGER.jsonl"
    if not p.exists():
        return 0
    n = 0
    for line in p.read_bytes().splitlines():
        if not line.strip():
            continue
        try:
            rec = json.loads(line.decode("utf-8"))
        except ValueError:
            continue
        if rec.get("plan_normalized_sha256") == plan_normalized_sha256 \
                and rec.get("event") == "run_started":
            n += 1
    return n


# ------------------------------------------------------------------------------------- stages
def bench_root_from(reg, override=None):
    root = override or os.environ.get("SPEXTRACTOR_BENCH_ROOT") or reg.get("bench_root")
    if not root:
        die("G-B00", "no bench root: set `bench_root:` in samples.yaml or $SPEXTRACTOR_BENCH_ROOT")
    return Path(root)


def do_convert(bench, plan, reg, sample, arm, tool_pre, refs, entry, pinned):
    raw_dir = reg_get(entry, "raw", "samples.%s" % sample)
    ident = raw_identity(raw_dir)
    # G-B05: content, not a path substring. The pinned digest is what makes a retargeted symlink
    # a hard failure at recipe-construction time instead of a wrong number six hours later.
    if ident["content_id"] != pinned["raw_content_id"]:
        die("G-B05", "sample %s raw content digest %s does not match the pin %s.\n"
                     "       Either the .d changed or the path now resolves elsewhere. v1 checked "
                     "only that the string %r appeared in the path (FAILURE 1)."
            % (sample, short(ident["content_id"]), short(pinned["raw_content_id"]), sample))
    if ident["acquisition_id"] != pinned["acquisition_id"]:
        die("G-B05", "sample %s acquisition identity %s does not match the pin %s"
            % (sample, short(ident["acquisition_id"]), short(pinned["acquisition_id"])))

    recipe = {
        "recipe_version": CONTRACT_VERSION,
        "kind": "convert",
        "harness_version": HARNESS_VERSION,
        "harness_code_digest": harness_code_digest(),
        "tool": tool_pre,
        "input": {"sample": sample, "raw_realpath": ident["realpath"],
                  "raw_content_id": ident["content_id"],
                  "raw_coverage": ident["coverage"],
                  "acquisition_id": ident["acquisition_id"],
                  "method_id": ident["method_id"]},
        "config": {"explicit_params": dict((k, str(v)) for k, v in sorted(arm["params"].items())),
                   "threads_requested": plan["threads"]},
        "arm": arm["name"],
    }
    att = Attempt(bench, "convert", recipe)
    ini = att.work / "resolved.ini"
    applied = write_resolved_ini(tool_pre["realpath"], arm["params"], plan["threads"], ini)
    ini_items, ini_sha = canonical_ini(ini)
    # The recipe's config digest is the hash of the config that will literally execute, because
    # the very next line runs the tool with -ini pointing at this file.
    recipe["config"]["ini_sha256"] = ini_sha
    recipe["config"]["applied"] = applied
    (att.dir / "recipe.json").write_bytes(canon(recipe))
    att.recipe = recipe
    att.recipe_id = sha256_bytes(canon(recipe))
    att.stage_id = "convert-%s-%s" % (att.recipe_id[:16], att.attempt_id)

    out = att.work / "pseudo.mzML"
    argv = validate_converter_argv([tool_pre["realpath"], "-ini", str(ini),
                                    "-in", ident["realpath"], "-out", str(out)])
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(plan["threads"])
    if plan["tcmalloc"]:
        tc = "/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"
        if Path(tc).exists():
            env["LD_PRELOAD"] = tc
            recipe.setdefault("env", {})["LD_PRELOAD"] = tc
            recipe["env"]["ld_preload_sha256"] = sha256_file(tc)

    att.run_process(argv, env, "run.log")
    assert_binary_unchanged(tool_pre["realpath"], tool_pre)
    verify_raw_unchanged(raw_dir, ident)

    if not out.exists():
        dest = att.fail("missing_output", "pseudo.mzML")
        die("G-B11", "converter exited 0 but produced no pseudo.mzML (attempt kept at %s)" % dest)
    info = hash_and_validate_mzml(out)
    att.declare("work/pseudo.mzML", {"sha256": info["sha256"], "bytes": info["bytes"]})
    # The ARTIFACT digest is the file's own sha256. `ini_sha256` is the digest of the CANONICAL
    # typed parameter map extracted from it -- a different thing, used for no-op detection. Using
    # one where the other belongs makes verify --deep fail on an untouched file.
    att.declare("work/resolved.ini",
                {"sha256": sha256_file(ini), "bytes": ini.stat().st_size,
                 "canonical_sha256": ini_sha}, "config")
    att.checks = {"spectrum_count_declared": info["spectrum_count_declared"],
                  "spectrum_count_seen": info["spectrum_count_seen"],
                  "source_file_locations": info["source_file_locations"],
                  "ini_sha256": ini_sha,
                  "method_id": ident["method_id"]}
    return att.commit({"ini_sha256": ini_sha, "method_id": ident["method_id"],
                       "method": ident["method"], "raw_realpath": ident["realpath"],
                       "sample": sample, "arm": arm["name"]})


def do_adopt(bench, plan, reg, sample, arm, refs, entry, pinned):
    """An external mzML (the reference implementation) becomes a first-class arm.

    v1 registered reference_mzml, printed its existence in preflight, and never read it -- so
    every published SpeXtractor-vs-the reference implementation ratio entered through an unguarded external procedure.
    An adopt stage gives the external file a sealed identity so the SAME search, scoring and
    truth stages run against it.
    """
    key = arm["source"]
    if not key.startswith("references."):
        die("G-B01", "arm %r source %r must be references.<name>" % (arm["name"], key))
    name = key.split(".", 1)[1]
    if name not in refs:
        die("G-B06", "arm %r wants reference %r for sample %s, which the registry does not "
                     "declare" % (arm["name"], name, sample))
    src = Path(refs[name]["path"])
    # An external mzML has no converter, but it is still SCORED against this sample's own
    # acquisition geometry -- which is read from this sample's .d, exactly as for a convert arm.
    raw_dir = reg_get(entry, "raw", "samples.%s" % sample)
    acq = acquisition_identity(Path(raw_dir).resolve())
    if acq["acquisition_id"] != pinned["acquisition_id"]:
        die("G-B05", "sample %s acquisition identity drifted from its pin" % sample)
    recipe = {
        "recipe_version": CONTRACT_VERSION,
        "kind": "adopt",
        "harness_version": HARNESS_VERSION,
        "harness_code_digest": harness_code_digest(),
        "input": {"sample": sample, "source_key": key, "source_path": str(src),
                  "source_sha256": refs[name]["sha256"],
                  "raw_realpath": str(Path(raw_dir).resolve()),
                  "acquisition_id": pinned["acquisition_id"],
                  "method_id": pinned["method_id"]},
        "arm": arm["name"],
    }
    att = Attempt(bench, "adopt", recipe)
    dest = att.work / "pseudo.mzML"
    shutil.copyfile(str(src), str(dest))
    info = hash_and_validate_mzml(dest)
    if info["sha256"] != refs[name]["sha256"]:
        att.fail("copy_mismatch", info["sha256"])
        die("G-B06", "adopted mzML digest %s != pinned %s"
            % (short(info["sha256"]), short(refs[name]["sha256"])))
    att.declare("work/pseudo.mzML", {"sha256": info["sha256"], "bytes": info["bytes"]})
    att.checks = {"spectrum_count_declared": info["spectrum_count_declared"],
                  "spectrum_count_seen": info["spectrum_count_seen"],
                  "adopted_from": str(src), "method_id": pinned["method_id"]}
    return att.commit({"ini_sha256": None, "method_id": pinned["method_id"],
                       "method": acq["method"],
                       "raw_realpath": str(Path(raw_dir).resolve()),
                       "sample": sample, "arm": arm["name"]})


SAGE_REQUIRED_COLUMNS = ("psm_id", "filename", "scannr", "rank", "label", "peptide",
                         "charge", "spectrum_q", "peptide_q", "protein_q", "hyperscore")


def do_search(bench, reg, sample, arm, produce_kind, produce_stage_id, search_tools):
    pdir, pseal = bench.load_stage(produce_kind, produce_stage_id)
    mzml_src = pdir / "work" / "pseudo.mzML"
    mzml_digest = pseal["artifacts"]["work/pseudo.mzML"]["sha256"]
    recipe = {
        "recipe_version": CONTRACT_VERSION,
        "kind": "search",
        "harness_version": HARNESS_VERSION,
        "harness_code_digest": harness_code_digest(),
        # Named by DIGEST, not by path. This is the link that makes cross-attempt assembly
        # impossible: score's input digest must equal search's output digest, which must equal
        # convert's output digest, all inside one cell.
        "input": {"produce_kind": produce_kind, "produce_stage_id": produce_stage_id,
                  "pseudo_mzml_sha256": mzml_digest, "sample": sample},
        "sage": {"binary_sha256": search_tools["sage_binary"]["sha256"],
                 "config_sha256": search_tools["sage_config"]["sha256"],
                 "fasta_sha256": search_tools["fasta"]["sha256"]},
        "arm": arm["name"],
    }
    att = Attempt(bench, "search", recipe)
    sage = search_tools["sage_binary"]["path"]
    argv = [sage, search_tools["sage_config"]["path"], "-o", str(att.work),
            "--disable-telemetry-i-dont-want-to-improve-sage", str(mzml_src)]
    att.run_process(argv, dict(os.environ), "sage.log")
    tsv = att.work / "results.sage.tsv"
    if not tsv.exists():
        dest = att.fail("missing_output", "results.sage.tsv")
        die("G-B11", "sage exited 0 but wrote no results.sage.tsv (attempt kept at %s)" % dest)
    info = hash_and_validate_tsv(tsv, SAGE_REQUIRED_COLUMNS)
    att.declare("work/results.sage.tsv", {"sha256": info["sha256"], "bytes": info["bytes"]})
    att.checks = {"rows": info["rows"], "header_sha256": info["header_sha256"],
                  "header": info["header"]}
    return att.commit({"sample": sample, "arm": arm["name"],
                       "header_sha256": info["header_sha256"]})


# ------------------------------------------------------------------------------------- commands
def cmd_plan_check(args):
    reg, _ = load_registry()
    plan = load_plan(args.plan, reg)
    print("plan OK: %s" % plan["title"])
    print("  normalized digest %s" % short(plan["normalized_sha256"]))
    print("  samples: %s" % ", ".join(plan["samples"]))
    for a in plan["arms"]:
        print("  arm %-24s role=%-18s baseline=%-10s params=%s"
              % (a["name"], a["role"], a["baseline"] or "-", a["params"] or "{}"))
    for s in plan["samples"]:
        entry, pinned = require_pinned(reg, s)
        verify_pinned_references(reg, s)
        print("  sample %-5s raw=%s acq=%s method=%s"
              % (s, short(pinned["raw_content_id"]), short(pinned["acquisition_id"]),
                 pinned["method_id"]))
    verify_search_toolchain(reg)
    print("preflight would pass. No directories were created.")
    return 0


def cmd_pin(args):
    try:
        import yaml
    except ImportError:
        die("G-B00", "need pyyaml")
    reg, regpath = load_registry()
    sample = args.sample
    entry = sample_entry(reg, sample)
    raw = reg_get(entry, "raw", "samples.%s" % sample)
    ident = raw_identity(raw)
    refs = entry.get("references") or {}
    new_refs = {}
    for name in sorted(refs):
        rec = refs[name]
        path = Path(rec["path"] if isinstance(rec, dict) else rec)
        if not path.exists():
            die("G-B06", "cannot pin %s.%s: %s does not exist" % (sample, name, path))
        new_refs[name] = {"path": str(path), "sha256": sha256_file(path)}
    entry["pinned"] = {"raw_content_id": ident["content_id"],
                       "acquisition_id": ident["acquisition_id"],
                       "method_id": ident["method_id"],
                       "raw_coverage": ident["coverage"],
                       "pinned_utc": utc_now(),
                       "pinned_by": "%s@%s" % (os.environ.get("USER", "?"), platform.node())}
    entry["references"] = new_refs
    entry.pop("fingerprint", None)  # superseded by acquisition_id; a substring is not an identity
    regpath.write_text(yaml.safe_dump(reg, default_flow_style=False, sort_keys=True))
    print("pinned %s: raw=%s acq=%s method=%s"
          % (sample, short(ident["content_id"]), short(ident["acquisition_id"]),
             ident["method_id"]))
    for name in sorted(new_refs):
        print("  ref %-16s %s" % (name, short(new_refs[name]["sha256"])))
    return 0


def cmd_run(args):
    reg, _ = load_registry()
    plan = load_plan(args.plan, reg)
    bench = Bench(bench_root_from(reg, args.bench_root))
    check_free_space(bench.root, float(args.need_gb))

    runset_id = "%s-%s" % (utc_stamp(), uuid.uuid4().hex[:12])
    lock = HostLock(bench.root / "locks" / "host.lock")
    lock.acquire(runset_id)
    try:
        ledger_append(bench, {"event": "run_started", "utc": utc_now(),
                              "runset_id": runset_id, "plan_path": plan["path"],
                              "plan_normalized_sha256": plan["normalized_sha256"],
                              "host": platform.node()})
        tool_pre = tool_identity(plan["binary"])
        search_tools = verify_search_toolchain(reg)
        prepared = {}
        for s in plan["samples"]:
            entry, pinned = require_pinned(reg, s)
            prepared[s] = (entry, pinned, verify_pinned_references(reg, s))

        # Run order is randomized within each sample and RECORDED. This mitigates the fixed-order
        # confound (cache warming, thermal state, node pressure masquerading as an arm effect).
        # It does NOT estimate variance -- n is still 1. See README UNCOVERED.
        import random
        seed = int(uuid.uuid4().int % (2 ** 31))
        rng = random.Random(seed)
        arm_order = {}
        cells = []
        for s in plan["samples"]:
            order = list(plan["arms"])
            rng.shuffle(order)
            arm_order[s] = [a["name"] for a in order]
            entry, pinned, refs = prepared[s]
            for arm in order:
                print("[run] %s/%s ..." % (s, arm["name"]))
                if arm["kind"] == "reference_mzml":
                    kind = "adopt"
                    stage_id, pseal = do_adopt(bench, plan, reg, s, arm, refs, entry, pinned)
                else:
                    kind = "convert"
                    stage_id, pseal = do_convert(bench, plan, reg, s, arm, tool_pre, refs,
                                                 entry, pinned)
                search_id, sseal = do_search(bench, reg, s, arm, kind, stage_id, search_tools)
                cells.append({"sample": s, "arm": arm["name"], "role": arm["role"],
                              "kind": arm["kind"], "baseline_arm": arm["baseline"],
                              "produce_kind": kind, "produce_stage_id": stage_id,
                              "produce_sha256": pseal["artifacts"]["work/pseudo.mzML"]["sha256"],
                              "ini_sha256": pseal.get("ini_sha256"),
                              "method_id": pseal.get("method_id"),
                              "search_stage_id": search_id,
                              "metrics_stage_id": None})
                print("      produce %s  search %s" % (short(stage_id.split("-")[-1]),
                                                       short(search_id.split("-")[-1])))
        detect_noop_arms(cells)

        obj = {
            "runset_version": CONTRACT_VERSION, "runset_id": runset_id, "revision": 1,
            "parent_revision_sha256": None,
            "created_utc": utc_now(), "status": "produced",
            "plan_path": plan["path"], "plan_title": plan["title"],
            "plan_file_sha256": plan["file_sha256"],
            "plan_normalized_sha256": plan["normalized_sha256"],
            # The plan bytes are EMBEDDED. If the plan file is later moved or edited, the runset
            # is still collatable; the on-disk file becomes optional corroboration.
            "plan_source_b64": base64.b64encode(plan["source_bytes"]).decode("ascii"),
            "plan_arms": plan["arms"], "plan_samples": plan["samples"],
            "harness_version": HARNESS_VERSION,
            "harness_code_digest": harness_code_digest(),
            "tool": tool_pre, "search_tools": search_tools,
            "host": host_facts(), "host_node": platform.node(),
            "arm_order": arm_order, "order_seed": seed,
            "cells": cells,
        }
        p = write_runset(bench, obj)
        ledger_append(bench, {"event": "run_completed", "utc": utc_now(),
                              "runset_id": runset_id,
                              "plan_normalized_sha256": plan["normalized_sha256"]})
        print("\n[runset] %s" % runset_id)
        print("         %s" % p)
        print("         next: collate2.py truth --runset %s" % runset_id)
        return 0
    finally:
        lock.release()


def detect_noop_arms(cells):
    """G-B14 -- an arm whose OUTPUT BYTES equal its baseline's is a no-op. Exact, zero false
    positives, no help-text parsing.

    FAILURE 5: three arms passed values that were already defaults and burned ~2.5 h producing
    byte-identical output. v1 tried to catch this by regexing --helphelp and diffing strings,
    which misses 0.1-vs-0.10, boolean aliases, aliased flags, tool-side clamping, and any
    parameter ignored under another mode.

    Two tests, both applied:
      (a) resolved INI equality  -- catches it BEFORE the run would be worth anything
      (b) output digest equality -- catches everything (a) misses, after the fact
    """
    by_sample = {}
    for c in cells:
        by_sample.setdefault(c["sample"], {})[c["arm"]] = c
    for sample, arms in sorted(by_sample.items()):
        for arm, cell in sorted(arms.items()):
            base = cell.get("baseline_arm")
            if not base or base not in arms:
                continue
            b = arms[base]
            if cell["produce_sha256"] == b["produce_sha256"]:
                die("G-B14", "arm %r on %s produced output BYTE-IDENTICAL to its baseline %r "
                             "(%s).\n       This is FAILURE 5: the arm tested nothing. Output-"
                             "digest equality is exact -- there is no interpretation in which "
                             "these two runs differ."
                    % (arm, sample, base, short(cell["produce_sha256"])))
            if cell["ini_sha256"] and cell["ini_sha256"] == b["ini_sha256"]:
                die("G-B14", "arm %r on %s resolved to a configuration identical to baseline %r "
                             "(ini %s). FAILURE 5."
                    % (arm, sample, base, short(cell["ini_sha256"])))


def assert_cells_match_plan(runset):
    """G-B17 -- exact set equality between the plan's cells and the runset's cells.

    No missing, no extra, no duplicates. v1's truth/collate globbed every *__* directory under a
    root derived from the PLAN FILENAME and did not restrict to the plan's samples, arms or
    binary; a copied directory was scored, and an empty root produced a valid-looking empty
    report.
    """
    required = set()
    for s in runset["plan_samples"]:
        for a in runset["plan_arms"]:
            required.add((s, a["name"]))
    seen = {}
    for c in runset["cells"]:
        key = (c["sample"], c["arm"])
        if key in seen:
            die("G-B17", "runset contains TWO cells for (%s, %s). A rerun after a failure must "
                         "not leave two completed runs for one plan cell; exactly one is "
                         "required." % key)
        seen[key] = c
    if set(seen) != required:
        missing = sorted(required - set(seen))
        extra = sorted(set(seen) - required)
        die("G-B17", "runset does not match its plan.\n       missing cells: %s\n"
                     "       extra cells:   %s" % (missing, extra))
    return seen


def cmd_verify(args):
    reg, _ = load_registry()
    bench = Bench(bench_root_from(reg, args.bench_root))
    runset, path = load_runset(bench, args.runset)
    cells = assert_cells_match_plan(runset)
    nfiles = 0
    nbytes = 0
    for (sample, arm), c in sorted(cells.items()):
        stages = [(c["produce_kind"], c["produce_stage_id"]), ("search", c["search_stage_id"])]
        if c.get("metrics_stage_id"):
            stages.append(("metrics", c["metrics_stage_id"]))
        for kind, sid in stages:
            d, seal = bench.load_stage(kind, sid)
            if sha256_file(d / "recipe.json") != seal["recipe_sha256"]:
                die("G-B12", "%s/%s: recipe.json does not match its sealed digest" % (kind, sid))
            for name, rec in sorted(seal["artifacts"].items()):
                p = d / name
                if not p.exists():
                    die("G-B12", "%s/%s: sealed artifact %s is gone" % (kind, sid, name))
                if p.stat().st_size != rec["bytes"]:
                    die("G-B12", "%s/%s: %s changed size after sealing" % (kind, sid, name))
                if args.deep:
                    got = sha256_file(p)
                    if got != rec["sha256"]:
                        die("G-B12", "%s/%s: %s content changed after sealing (%s -> %s)"
                            % (kind, sid, name, short(rec["sha256"]), short(got)))
                    nbytes += rec["bytes"]
                nfiles += 1
        # Chain binding: each derived stage names its parent by digest.
        _, pseal = bench.load_stage(c["produce_kind"], c["produce_stage_id"])
        sdir, sseal = bench.load_stage("search", c["search_stage_id"])
        srecipe = jload(sdir / "recipe.json")
        if srecipe["input"]["pseudo_mzml_sha256"] != \
                pseal["artifacts"]["work/pseudo.mzML"]["sha256"]:
            die("G-B12", "cell (%s,%s): the search stage was run against a DIFFERENT mzML than "
                         "the convert stage in this cell produced. This is the three-attempts "
                         "failure." % (sample, arm))
    mode = "deep" if args.deep else "fast"
    print("verify(%s) OK: %d cells, %d sealed artifacts%s"
          % (mode, len(cells), nfiles,
             (", %.2f GB rehashed" % (nbytes / (1024.0 ** 3))) if args.deep else ""))
    print("runset %s r%d  plan %s  host %s"
          % (runset["runset_id"], runset["revision"],
             short(runset["plan_normalized_sha256"]), runset["host_node"]))
    return 0


def cmd_show(args):
    reg, _ = load_registry()
    bench = Bench(bench_root_from(reg, args.bench_root))
    runset, path = load_runset(bench, args.runset)
    print("runset %s revision %d  status=%s" % (runset["runset_id"], runset["revision"],
                                                runset["status"]))
    print("plan    %s (%s)" % (runset["plan_title"], short(runset["plan_normalized_sha256"])))
    print("tool    %s  host %s" % (short(runset["tool"]["sha256"]), runset["host_node"]))
    print("order   %s (seed %d)" % (runset["arm_order"], runset["order_seed"]))
    for c in runset["cells"]:
        print("  %-5s %-24s produce=%s search=%s metrics=%s"
              % (c["sample"], c["arm"], short(c["produce_stage_id"].split("-")[-1]),
                 short(c["search_stage_id"].split("-")[-1]),
                 short(c["metrics_stage_id"].split("-")[-1]) if c.get("metrics_stage_id")
                 else "-"))
    return 0


# ============================================================================== TEST FIXTURES
# Prefixed st_ and used only by --selftest in this file and in collate2.py.

def st_build_mzml(path, specs, source_location, count_override=None, truncate=False,
                  terminator=True):
    """Emit a small but genuinely valid mzML. Verified against pyteomics."""
    def b64(vals):
        return base64.b64encode(struct.pack("<%dd" % len(vals), *vals)).decode("ascii")

    parts = ['<?xml version="1.0" encoding="utf-8"?>',
             '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">',
             ' <cvList count="1"><cv id="MS" fullName="PSI-MS" URI="http://psidev.info"/></cvList>',
             ' <fileDescription><fileContent/><sourceFileList count="1">',
             '  <sourceFile id="RAW1" name="%s" location="%s"/>'
             % (os.path.basename(source_location), source_location),
             ' </sourceFileList></fileDescription>',
             ' <run id="r1">',
             '  <spectrumList count="%d" defaultDataProcessingRef="dp">'
             % (count_override if count_override is not None else len(specs))]
    for i, sp in enumerate(specs):
        mzs = [p[0] for p in sp["peaks"]]
        ins = [p[1] for p in sp["peaks"]]
        emz, ein = b64(mzs), b64(ins)
        rt_unit = sp.get("rt_unit", "second")
        if rt_unit is None:
            rt_attr = ''
        elif rt_unit == "second":
            rt_attr = ' unitAccession="UO:0000010" unitName="second"'
        elif rt_unit == "minute":
            rt_attr = ' unitAccession="UO:0000031" unitName="minute"'
        else:
            rt_attr = ' unitAccession="XX:9999999" unitName="%s"' % rt_unit
        parts.append('   <spectrum index="%d" id="scan=%d" defaultArrayLength="%d">'
                     % (i, i + 1, len(mzs)))
        parts.append('    <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="%d"/>'
                     % sp.get("ms_level", 2))
        parts.append('    <scanList count="1"><scan>')
        parts.append('     <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" '
                     'value="%s"%s/>' % (sp["rt"], rt_attr))
        parts.append('    </scan></scanList>')
        parts.append('    <precursorList count="1"><precursor><isolationWindow>')
        if sp.get("iso_target") is not None:
            parts.append('     <cvParam cvRef="MS" accession="MS:1000827" name="isolation window '
                         'target m/z" value="%s"/>' % sp["iso_target"])
        parts.append('    </isolationWindow><selectedIonList count="1"><selectedIon>')
        if sp.get("mz") is not None:
            parts.append('     <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" '
                         'value="%s"/>' % sp["mz"])
        if sp.get("z") is not None:
            parts.append('     <cvParam cvRef="MS" accession="MS:1000041" name="charge state" '
                         'value="%d"/>' % sp["z"])
        if sp.get("im") is not None:
            parts.append('     <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced '
                         'ion mobility" value="%s"/>' % sp["im"])
        parts.append('    </selectedIon></selectedIonList></precursor></precursorList>')
        parts.append('    <binaryDataArrayList count="2">')
        for acc, nm, enc, unit in (("MS:1000514", "m/z array", emz,
                                    ' unitAccession="MS:1000040" unitName="m/z"'),
                                   ("MS:1000515", "intensity array", ein,
                                    ' unitAccession="MS:1000131" unitName="number of counts"')):
            parts.append('     <binaryDataArray encodedLength="%d">' % len(enc))
            parts.append('      <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/>')
            parts.append('      <cvParam cvRef="MS" accession="MS:1000576" '
                         'name="no compression"/>')
            parts.append('      <cvParam cvRef="MS" accession="%s" name="%s"%s/>'
                         % (acc, nm, unit))
            parts.append('      <binary>%s</binary>' % enc)
            parts.append('     </binaryDataArray>')
        parts.append('    </binaryDataArrayList>')
        parts.append('   </spectrum>')
    parts.append('  </spectrumList>')
    parts.append(' </run>')
    if terminator:
        parts.append('</mzML>')
    text = "\n".join(parts) + "\n"
    if truncate:
        text = text[:int(len(text) * 0.6)]
    Path(path).write_text(text)
    return path


def st_make_fake_d(path, mz_centers=(500.0, 600.0), width=25.0, serial="SER-1",
                   dt="2026-07-01T10:00:00", sample_name="X", frames=10, num_scans=709,
                   k0=(0.6, 1.6), full_mobility_span=False):
    """A .d whose analysis.tdf is a real sqlite database with the tables we read."""
    d = Path(path)
    d.mkdir(parents=True, exist_ok=True)
    (d / "analysis.tdf_bin").write_bytes(b"\x00" * 4096)
    tdf = d / "analysis.tdf"
    if tdf.exists():
        tdf.unlink()
    con = sqlite3.connect(str(tdf))
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)", [
        ("InstrumentSerialNumber", serial), ("AcquisitionDateTime", dt),
        ("SampleName", sample_name),
        ("OneOverK0AcqRangeLower", str(k0[0])), ("OneOverK0AcqRangeUpper", str(k0[1]))])
    con.execute("CREATE TABLE Frames (Id INTEGER, NumScans INTEGER)")
    con.executemany("INSERT INTO Frames VALUES (?,?)",
                    [(i + 1, num_scans) for i in range(frames)])
    con.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, "
                "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL)")
    rows = []
    for i, c in enumerate(mz_centers):
        if full_mobility_span:
            sb, se = 0, num_scans - 1
        else:
            sb, se = 0 + i * 100, 300 + i * 100
        rows.append((i + 1, sb, se, float(c), float(width)))
    con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return d


FAKETOOL_SRC = r'''#!/usr/bin/env python3
"""Stand-in for spextractor, driven by $FAKETOOL_MODE. Used only by --selftest."""
import os, sys, xml.etree.ElementTree as ET
sys.path.insert(0, %(harness)r)
import bench2

def write_ini(dest, extra):
    root = ET.Element("PARAMETERS", {"version": "1.7.0"})
    node = ET.SubElement(root, "NODE", {"name": "trace"})
    ET.SubElement(node, "ITEM", {"name": "ms1_split_valleys", "value": "0.0", "type": "double"})
    ET.SubElement(node, "ITEM", {"name": "ambiguity_margin", "value": "0.0", "type": "double"})
    ET.SubElement(root, "ITEM", {"name": "threads", "value": "1", "type": "int"})
    ET.ElementTree(root).write(dest, encoding="UTF-8", xml_declaration=True)

def read_ini(path):
    vals = {}
    for it in ET.parse(path).getroot().iter():
        if it.tag.endswith("ITEM"):
            vals[it.get("name")] = it.get("value")
    return vals

def main():
    a = sys.argv[1:]
    mode = os.environ.get("FAKETOOL_MODE", "ok")
    if "-write_ini" in a:
        write_ini(a[a.index("-write_ini") + 1], None)
        return 0
    ini = a[a.index("-ini") + 1] if "-ini" in a else None
    raw = a[a.index("-in") + 1] if "-in" in a else ""
    out = a[a.index("-out") + 1] if "-out" in a else ""
    vals = read_ini(ini) if ini else {}
    if mode == "fail":
        sys.stderr.write("simulated converter failure\n")
        return 3
    if mode == "mutate_self":
        with open(__file__, "a") as fh:
            fh.write("# mutated mid-run\n")
    if mode == "retarget":
        import sqlite3
        con = sqlite3.connect(os.path.join(raw, "analysis.tdf"))
        con.execute("UPDATE GlobalMetadata SET Value='OTHER' "
                    "WHERE Key='InstrumentSerialNumber'")
        con.commit(); con.close()
    n = int(os.environ.get("FAKETOOL_SPECTRA", "3"))
    if mode == "noop":
        n = max(n - 1, 1)
    elif float(vals.get("ms1_split_valleys", "0") or 0) > 0:
        n = n + 1
    mz0 = float(os.environ.get("FAKETOOL_MZ0", "500.0"))
    step = float(os.environ.get("FAKETOOL_MZSTEP", "1.0"))
    im = float(os.environ.get("FAKETOOL_IM", "1.0"))
    indexed = os.environ.get("FAKETOOL_PEAKMODE") == "indexed"
    specs = []
    for i in range(n):
        k = i %% 300
        if indexed:
            # Fragment layout chosen so that no accidental coincidence exists between the
            # target reference sets and their +11.003 Da shifts: bases are 7.0 Da apart and
            # fragments 0.37 Da apart, so 11.003 never lands within 20 ppm of a real fragment.
            b = 1000.0 + k * 7.0
            peaks = ([(b + j * 0.37, 100.0) for j in range(6)] +
                     [(b + 21.0 + j * 0.37, 50.0) for j in range(6)] +
                     [(b + 11.003 + j * 0.37, 5.0) for j in range(4)] +
                     [(b + 32.003 + j * 0.37, 4.0) for j in range(4)])
            spec_im = 0.7 + k * 0.002
        else:
            peaks = [(300.0 + j, 10.0 * (j + 1)) for j in range(25)]
            spec_im = im
        specs.append({"mz": mz0 + k * step, "z": 2, "im": spec_im, "rt": 600.0,
                      "iso_target": mz0, "peaks": peaks})
    bench2.st_build_mzml(out, specs, "file://" + os.path.abspath(raw),
                         count_override=(n + 1) if mode == "partial" else None,
                         truncate=(mode == "truncate"))
    return 0

sys.exit(main())
'''


def st_write_faketool(path):
    Path(path).write_text(FAKETOOL_SRC % {"harness": str(HERE)})
    os.chmod(str(path), 0o755)
    return path


FAKESAGE_SRC = r'''#!/usr/bin/env python3
"""Stand-in for sage, driven by $FAKESAGE_MODE."""
import os, sys
mode = os.environ.get("FAKESAGE_MODE", "ok")
a = sys.argv[1:]
outdir = a[a.index("-o") + 1]
if mode == "fail":
    sys.stderr.write("simulated sage failure\n"); sys.exit(4)
cols = ["psm_id","filename","scannr","rank","label","peptide","charge",
        "spectrum_q","peptide_q","protein_q","hyperscore"]
if mode == "drop_peptide_q":
    cols.remove("peptide_q")
rows = []
for i in range(12):
    rec = {"psm_id": str(i), "filename": "pseudo.mzML", "scannr": "scan=%d" % i, "rank": "1",
           "label": "1", "peptide": "PEPTIDE%d" % (i % 5), "charge": "2",
           "spectrum_q": "0.001", "peptide_q": "0.002", "protein_q": "0.001",
           "hyperscore": "30.0"}
    if mode == "na_qvalue" and i == 3:
        rec["peptide_q"] = "NA"
    rows.append("\t".join(rec[c] for c in cols))
body = "\t".join(cols) + "\n" + "\n".join(rows) + "\n"
if mode == "no_trailing_newline":
    body = body.rstrip("\n")
if mode == "ragged":
    body = body.replace("\nPEPTIDE", "\nEXTRA\tPEPTIDE", 1)
    body = "\t".join(cols) + "\n" + "\n".join(rows[:-1]) + "\n" + rows[-1] + "\tSURPLUS\n"
open(os.path.join(outdir, "results.sage.tsv"), "w").write(body)
sys.exit(0)
'''


def st_write_fakesage(path):
    Path(path).write_text(FAKESAGE_SRC)
    os.chmod(str(path), 0o755)
    return path


# ================================================================================== SELF TEST
class SelfTest(object):
    """--selftest deliberately triggers each guard and verifies THAT guard fired.

    The v1 harness was committed without a single guard ever having been demonstrated to fire.
    Six turned out to be theatre. A guard that has never been observed to abort is a comment.
    """

    def __init__(self, label):
        self.label = label
        self.results = []

    def expect_abort(self, guard, name, fn, match=None):
        """Assert that `fn` aborts with `guard` -- and, when `match` is given, that the
        SPECIFIC check fired.

        Matching on the guard id alone is a theatre mechanism in its own right. Mutation
        testing showed that disabling `raw not found` at bench2.py:426 still left this
        assertion green, because a sibling die("G-B13") further down fired instead and
        reported an unrelated condition. The suite could not tell the two apart, so a
        distinct check could be deleted with the tests still passing. `match` pins the
        assertion to the message, making each site independently killable.
        """
        try:
            fn()
        except Abort as exc:
            if exc.guard != guard:
                self.results.append(("FAIL", guard, name,
                                     "wrong guard fired: %s" % exc.guard))
            elif match and match.lower() not in exc.message.lower():
                self.results.append(("FAIL", guard, name,
                                     "right guard, WRONG CHECK: expected %r, got %r"
                                     % (match, exc.message.splitlines()[0][:60])))
            else:
                self.results.append(("PASS", guard, name, exc.message.splitlines()[0][:78]))
            return
        except Exception as exc:  # a crash is NOT a guard firing
            self.results.append(("FAIL", guard, name,
                                 "crashed instead of aborting: %s: %s"
                                 % (type(exc).__name__, exc)))
            return
        self.results.append(("FAIL", guard, name, "guard did NOT fire"))

    def expect_ok(self, guard, name, fn):
        """Positive control. A harness that refuses everything gets bypassed by hand, and then
        it guarantees nothing -- so every abort path is paired with a case that must NOT abort."""
        try:
            fn()
        except Abort as exc:
            self.results.append(("FAIL", guard, name, "false positive: %s" % exc.message[:70]))
            return
        except Exception as exc:
            self.results.append(("FAIL", guard, name, "crashed: %s: %s"
                                 % (type(exc).__name__, exc)))
            return
        self.results.append(("PASS", guard, name, "no false positive"))

    def report(self):
        width = max(len(r[2]) for r in self.results) if self.results else 10
        print("\n%s" % self.label)
        print("=" * (width + 34))
        for status, guard, name, detail in self.results:
            print("%-4s %-7s %-*s  %s" % (status, guard, width, name, detail))
        npass = sum(1 for r in self.results if r[0] == "PASS")
        nfail = len(self.results) - npass
        print("=" * (width + 34))
        print("%d passed, %d FAILED, %d guard cases total" % (npass, nfail, len(self.results)))
        return 0 if nfail == 0 else 1


def _write_yaml(path, obj):
    import yaml
    Path(path).write_text(yaml.safe_dump(obj, default_flow_style=False, sort_keys=True))


def _selftest_env(tmp):
    """Build a complete synthetic world: registry, raws, tool, sage, bench root."""
    tmp = Path(tmp)
    raws = tmp / "raw"
    raws.mkdir()
    d08 = st_make_fake_d(raws / "dataset A.d", mz_centers=(500.0, 600.0), serial="SER-08",
                         sample_name="dataset A")
    d30 = st_make_fake_d(raws / "dataset D.d", mz_centers=(500.0, 600.0), serial="SER-30",
                         sample_name="dataset D")
    tool = st_write_faketool(tmp / "faketool.py")
    sage = st_write_fakesage(tmp / "fakesage.py")
    fasta = tmp / "human.fasta"
    fasta.write_text(">sp|X|X\nPEPTIDE\n")
    sage_cfg = tmp / "sage.json"
    sage_cfg.write_text(json.dumps({"database": {"fasta": str(fasta)}}))
    refs = {}
    for name, d in (("dataset A", d08), ("dataset D", d30)):
        rep = tmp / ("%s_report.parquet" % name)
        rep.write_bytes(b"PAR1placeholder")
        lib = tmp / ("%s_lib.parquet" % name)
        lib.write_bytes(b"PAR1placeholderlib")
        refs[name] = {"diann_report": {"path": str(rep), "sha256": sha256_file(rep)},
                      "diann_lib": {"path": str(lib), "sha256": sha256_file(lib)}}
    reg = {
        "bench_root": str(tmp / "bench"),
        "vault": str(tmp / "vault"),
        "tool": {"scaled_params": ["trace:ambiguity_margin"]},
        "search": {
            "sage_binary": {"path": str(sage), "sha256": sha256_file(sage)},
            "sage_config": {"path": str(sage_cfg), "sha256": sha256_file(sage_cfg)},
            "fasta": {"path": str(fasta), "sha256": sha256_file(fasta)},
        },
        "samples": {
            "dataset A": {"raw": str(d08), "references": refs["dataset A"]},
            "dataset D": {"raw": str(d30), "references": refs["dataset D"]},
        },
    }
    for name, d in (("dataset A", d08), ("dataset D", d30)):
        ident = raw_identity(d)
        reg["samples"][name]["pinned"] = {
            "raw_content_id": ident["content_id"],
            "acquisition_id": ident["acquisition_id"],
            "method_id": ident["method_id"],
            "raw_coverage": ident["coverage"],
            "pinned_utc": utc_now(), "pinned_by": "selftest"}
    return {"tmp": tmp, "reg": reg, "tool": tool, "sage": sage, "d08": d08, "d30": d30}


def _plan(env, arms, samples=("dataset A",), sweeps=None, extra=None):
    p = {"schema": "spextractor.plan/2", "title": "selftest", "binary": str(env["tool"]),
         "threads": 2, "tcmalloc": False, "samples": list(samples), "arms": arms}
    if sweeps:
        p["sweeps"] = sweeps
    if extra:
        p.update(extra)
    path = env["tmp"] / ("plan_%s.yaml" % uuid.uuid4().hex[:6])
    _write_yaml(path, p)
    return path


BASE_ARM = {"name": "base", "role": "baseline", "params": {}}
TREAT_ARM = {"name": "split", "role": "treatment", "baseline": "base",
             "params": {"trace:ms1_split_valleys": 7.0}}


def selftest_bench(env):
    st = SelfTest("bench2.py -- run identity, plan binding, artifact production")
    reg = env["reg"]
    tmp = env["tmp"]

    # ---- G-B20 canonical serializer refuses NaN -------------------------------------------
    st.expect_abort("G-B20", "NaN cannot enter an artifact",
                    lambda: canon({"x": float("nan")}))
    st.expect_ok("G-B20", "finite floats serialize", lambda: canon({"x": 1.5}))

    # ---- G-B13 raw identity refuses a missing / retargeted .d -----------------------------
    # Was PURE THEATRE: mutation sweep killed 0 of 6 G-B13 sites. This calls the PRODUCTION
    # function rather than reimplementing its comparison, so deleting the guard fails the suite.
    st.expect_abort("G-B13", "raw_identity on a nonexistent .d aborts",
                    lambda: raw_identity(str(tmp / "no-such-run.d")), match="raw not found")
    # and a path that exists but is a plain file, not a .d directory
    _f = tmp / "not_a_dir.d"
    _f.write_text("x")
    # a plain file passes exists() and falls through to the empty-members check at :447 --
    # asserted against what the code ACTUALLY does, which message-pinning is how I found out
    st.expect_abort("G-B13", "raw_identity on a file with no .d members aborts",
                    lambda: raw_identity(str(_f)), match="contains no files")

    # ---- G-B13 post-run stability: the TOCTOU half ---------------------------------------
    # Was UNVERIFIED (mutation: sites 567/581 survived). These are the sites the critical gate
    # names, because they cover the retarget-during-the-run hole: the recipe pins dataset A at
    # construction and the converter reads the .d for HOURS afterwards. Both call the PRODUCTION
    # verify_raw_unchanged() and are message-pinned, so neither can be satisfied by a sibling.
    def _mk_d(path, iso_mz):
        """Minimal but REAL .d: acquisition_identity reads GlobalMetadata and
        DiaFrameMsMsWindows out of a genuine sqlite analysis.tdf, so a byte fixture is not
        enough. Two tiles are sufficient to exercise the window table."""
        path.mkdir(exist_ok=True)
        tdf = path / "analysis.tdf"
        if tdf.exists():
            tdf.unlink()
        con = sqlite3.connect(str(tdf))
        con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
        con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)",
                        [("InstrumentName", "timsTOF"), ("AcquisitionSoftware", "selftest"),
                         ("OneOverK0AcqRangeLower", "0.6"), ("OneOverK0AcqRangeUpper", "1.6")])
        con.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INT, ScanNumBegin INT, "
                    "ScanNumEnd INT, IsolationMz REAL, IsolationWidth REAL)")
        con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?)",
                        [(1, 34, 602, iso_mz, 25.0), (1, 602, 944, iso_mz + 300.0, 25.0)])
        con.execute("CREATE TABLE Frames (Id INT, Time REAL, MsMsType INT, NumScans INT)")
        con.executemany("INSERT INTO Frames VALUES (?,?,?,?)",
                        [(1, 0.5, 0, 944), (2, 1.2, 9, 944)])
        con.commit(); con.close()
        (path / "analysis.tdf_bin").write_bytes(b"spectra")
        return path

    _d = _mk_d(tmp / "stable.d", 700.0)
    _pin = raw_identity(str(_d))

    # unchanged input must NOT abort -- a guard that refuses everything gets bypassed by hand
    st.expect_ok("G-B13", "an unmodified .d passes the post-run check",
                 lambda: verify_raw_unchanged(str(_d), _pin))

    # :567 -- the path now resolves somewhere else (symlink retargeted mid-run)
    _other = _mk_d(tmp / "other.d", 900.0)
    _link = tmp / "link.d"
    if _link.is_symlink() or _link.exists():
        _link.unlink()
    _link.symlink_to(_other)
    st.expect_abort("G-B13", "raw path retargeted during the run is detected",
                    lambda: verify_raw_unchanged(str(_link), _pin),
                    match="RETARGETED")

    # :581 -- same path, but a member's bytes changed under us
    (_d / "analysis.tdf_bin").write_bytes(b"spectra-MODIFIED-mid-run")
    st.expect_abort("G-B13", "a raw member modified during the run is detected",
                    lambda: verify_raw_unchanged(str(_d), _pin),
                    match="changed during the run")

    # :577 -- a member APPEARED or vanished (a partial re-copy of the .d mid-run)
    (_d / "analysis.tdf_bin").write_bytes(b"spectra")          # undo the :581 modification
    (_d / "stray.bin").write_bytes(b"appeared mid-run")
    st.expect_abort("G-B13", "a raw member added during the run is detected",
                    lambda: verify_raw_unchanged(str(_d), _pin),
                    match="gained or lost files")
    (_d / "stray.bin").unlink()

    # :586 -- acquisition identity drift. Rewriting the window table also changes the file's
    # bytes, so a physical simulation trips :581 first (message-pinning caught exactly that).
    # The check under test is the COMPARISON, so the pin is constructed with a wrong
    # acquisition_id and correct stat pins: every earlier check passes and only :586 can fire.
    # FRESH fixture: earlier cases mutated _d, and a stat pin includes mtime_ns, so even
    # rewriting identical bytes leaves _d permanently dirty against _pin.
    _acq = _mk_d(tmp / "acq.d", 700.0)
    _pin_bad = json.loads(json.dumps(raw_identity(str(_acq))))
    _pin_bad["acquisition_id"] = "f" * 64
    st.expect_abort("G-B13", "acquisition identity changing mid-run is detected",
                    lambda: verify_raw_unchanged(str(_acq), _pin_bad),
                    match="acquisition identity changed")

    # ---- G-B05 cross-sample: the PRODUCTION comparison, not a reimplementation -------------
    # FAILURE 1 -- a published ratio divided an dataset D result by dataset A's the reference implementation reference; the
    # class of error occurred 3x. The existing case for this guard copies the comparison into
    # its own body and calls die() itself, so disabling the production check at :1384 left the
    # suite fully green (mutation-verified). These call do_convert / do_adopt directly. The
    # digest check is the first statement after raw_identity(), so the later arguments are
    # never reached and may be None.
    _d0 = _mk_d(tmp / "b05.d", 700.0)          # no walrus: this interpreter is 3.7
    _id0 = raw_identity(str(_d0))
    # the PIN shape is not raw_identity's return shape -- do_convert reads raw_content_id
    _pin_ok = {"raw_content_id": _id0["content_id"],
               "acquisition_id": _id0["acquisition_id"],
               "realpath": str(Path(str(_d0)).resolve()),
               "stat_pins": _id0["stat_pins"]}
    _entry = {"raw": str(_d0)}

    def _bad_content():
        bad = json.loads(json.dumps(_pin_ok))
        bad["raw_content_id"] = "0" * 64
        do_convert(None, None, None, "dataset A", None, None, None, _entry, bad)
    st.expect_abort("G-B05", "raw content digest mismatch reaches the PRODUCTION check",
                    _bad_content, match="raw content digest")

    def _bad_acq():
        bad = json.loads(json.dumps(_pin_ok))
        bad["acquisition_id"] = "0" * 64
        do_convert(None, None, None, "dataset A", None, None, None, _entry, bad)
    st.expect_abort("G-B05", "acquisition identity mismatch reaches the PRODUCTION check",
                    _bad_acq, match="acquisition identity")

    # :1478 -- the adopt path (an external mzML such as the reference implementation's) is scored against THIS
    # sample's acquisition geometry, so it needs the same pin check.
    def _bad_adopt():
        bad = json.loads(json.dumps(_pin_ok))
        bad["acquisition_id"] = "0" * 64
        # arm is a dict whose "source" names the reference; the acquisition-pin check runs
        # after that resolution, so both must be well-formed to reach it
        do_adopt(None, None, None, "dataset A",
                 {"name": "dt", "source": "references.reference_mzml"},
                 {"reference_mzml": {"path": "/nonexistent"}}, _entry, bad)
    st.expect_abort("G-B05", "adopt path checks acquisition pin before using an external mzML",
                    _bad_adopt, match="acquisition identity drifted")

    # ---- G-B01 plan schema whitelist -------------------------------------------------------
    st.expect_abort("G-B01", "unknown top-level plan key", lambda: load_plan(
        _plan(env, [BASE_ARM, TREAT_ARM], extra={"outroot": "/tmp/x"}), reg))
    st.expect_abort("G-B01", "typo 'basline' is not silently ignored", lambda: load_plan(
        _plan(env, [BASE_ARM, {"name": "t", "role": "treatment", "basline": "base",
                               "params": {"trace:ms1_split_valleys": 7.0}}]), reg))
    st.expect_abort("G-B01", "wrong plan schema version", lambda: load_plan(
        _write_plan_raw(env, {"schema": "spextractor.plan/1", "title": "x",
                              "binary": str(env["tool"]), "samples": ["dataset A"],
                              "arms": [BASE_ARM]}), reg))

    # ---- G-B02 names -----------------------------------------------------------------------
    st.expect_abort("G-B02", "arm name containing '/'", lambda: load_plan(
        _plan(env, [BASE_ARM, {"name": "../escape", "role": "treatment", "baseline": "base",
                               "params": {"trace:ms1_split_valleys": 7.0}}]), reg))
    st.expect_abort("G-B02", "duplicate arm names", lambda: load_plan(
        _plan(env, [BASE_ARM, dict(TREAT_ARM), dict(TREAT_ARM)]), reg))

    # ---- G-B03 pairing ---------------------------------------------------------------------
    st.expect_abort("G-B03", "treatment with no baseline", lambda: load_plan(
        _plan(env, [BASE_ARM, {"name": "t", "role": "treatment",
                               "params": {"trace:ms1_split_valleys": 7.0}}]), reg))
    st.expect_abort("G-B03", "baseline names a missing arm", lambda: load_plan(
        _plan(env, [BASE_ARM, {"name": "t", "role": "treatment", "baseline": "nope",
                               "params": {"trace:ms1_split_valleys": 7.0}}]), reg))
    st.expect_ok("G-B03", "well-formed pair accepted",
                 lambda: load_plan(_plan(env, [BASE_ARM, TREAT_ARM]), reg))

    # ---- G-B04 scale / sweeps --------------------------------------------------------------
    st.expect_abort("G-B04", "scaled param set outside a sweep", lambda: load_plan(
        _plan(env, [BASE_ARM, {"name": "amb", "role": "treatment", "baseline": "base",
                               "params": {"trace:ambiguity_margin": 0.1}}]), reg))
    st.expect_abort("G-B04", "sweep with only 2 points", lambda: load_plan(
        _plan(env, [BASE_ARM], sweeps=[{"param": "trace:ambiguity_margin",
                                        "baseline": "base", "points": [0.1, 0.2]}]), reg))
    st.expect_ok("G-B04", "3-point sweep accepted", lambda: load_plan(
        _plan(env, [BASE_ARM], sweeps=[{"param": "trace:ambiguity_margin", "baseline": "base",
                                        "points": [0.1, 0.2, 0.4]}]), reg))

    # ---- G-B05 pinning / raw identity ------------------------------------------------------
    def unpinned():
        r = json.loads(json.dumps(reg))
        r["samples"]["dataset A"].pop("pinned")
        require_pinned(r, "dataset A")
    st.expect_abort("G-B05", "unpinned sample refused", unpinned)

    def wrong_pin():
        r = json.loads(json.dumps(reg))
        r["samples"]["dataset A"]["pinned"]["raw_content_id"] = "0" * 64
        entry, pinned = require_pinned(r, "dataset A")
        ident = raw_identity(entry["raw"])
        if ident["content_id"] != pinned["raw_content_id"]:
            die("G-B05", "raw content digest mismatch (selftest path)")
    st.expect_abort("G-B05", "raw content digest mismatch (symlink retarget)", wrong_pin)

    st.expect_abort("G-B05", "sample not in registry", lambda: require_pinned(reg, "S99"))

    # ---- G-B06 references / toolchain ------------------------------------------------------
    def bad_ref():
        r = json.loads(json.dumps(reg))
        r["samples"]["dataset A"]["references"]["diann_report"]["sha256"] = "1" * 64
        verify_pinned_references(r, "dataset A")
    st.expect_abort("G-B06", "reference digest changed since pinning", bad_ref)

    def missing_ref():
        r = json.loads(json.dumps(reg))
        r["samples"]["dataset A"]["references"]["diann_report"]["path"] = str(tmp / "nope.parquet")
        verify_pinned_references(r, "dataset A")
    st.expect_abort("G-B06", "reference file absent aborts (not 'MISSING' print)", missing_ref)

    def bare_path_ref():
        r = json.loads(json.dumps(reg))
        r["samples"]["dataset A"]["references"]["diann_lib"] = "/some/path.parquet"
        verify_pinned_references(r, "dataset A")
    st.expect_abort("G-B06", "bare-path reference (no digest) refused", bare_path_ref)

    def wrong_fasta():
        r = json.loads(json.dumps(reg))
        other = tmp / "other.fasta"
        other.write_text(">y\nY\n")
        cfg = tmp / "sage_bad.json"
        cfg.write_text(json.dumps({"database": {"fasta": str(other)}}))
        r["search"]["sage_config"] = {"path": str(cfg), "sha256": sha256_file(cfg)}
        verify_search_toolchain(r)
    st.expect_abort("G-B06", "sage_config names a non-registry FASTA", wrong_fasta)
    st.expect_ok("G-B06", "consistent toolchain accepted",
                 lambda: verify_search_toolchain(reg))

    # ---- G-B07 host lock -------------------------------------------------------------------
    def double_lock():
        lock_a = HostLock(tmp / "lockdir" / "host.lock").acquire("A")
        try:
            HostLock(tmp / "lockdir" / "host.lock").acquire("B")
        finally:
            lock_a.release()
    st.expect_abort("G-B07", "second concurrent run refused by flock", double_lock)

    # ---- G-B08 filesystem assumption -------------------------------------------------------
    class _FakeStat(object):
        pass

    def cross_device():
        real_stat = os.stat
        paths = [tmp / "a", tmp / "b"]
        for p in paths:
            p.mkdir(exist_ok=True)

        def fake(path, *a, **k):
            s = real_stat(path, *a, **k)

            class S(object):
                st_dev = 1 if str(path).endswith("a") else 2
            return S()
        os.stat = fake
        try:
            assert_same_filesystem(paths)
        finally:
            os.stat = real_stat
    st.expect_abort("G-B08", "staging and stages on different devices", cross_device)

    # ---- G-B19 disk space ------------------------------------------------------------------
    st.expect_abort("G-B19", "insufficient free space refuses to start",
                    lambda: check_free_space(tmp, 10 ** 9))
    st.expect_ok("G-B19", "adequate free space passes",
                 lambda: check_free_space(tmp, 0.001))

    # ---- G-B16 argv whitelist --------------------------------------------------------------
    st.expect_abort("G-B16", "converter argv carrying an extra flag", lambda:
                    validate_converter_argv(["tool", "-ini", "x", "-in", "y", "-out", "z",
                                             "-trace:ms1_split_valleys", "7.0"]))
    st.expect_ok("G-B16", "clean converter argv", lambda:
                 validate_converter_argv(["tool", "-ini", "x", "-in", "y", "-out", "z"]))
    st.expect_abort("G-B16", "unknown parameter absent from resolved INI", lambda:
                    write_resolved_ini(env["tool"], {"trace:does_not_exist": 1}, 2,
                                       tmp / "probe.ini"))

    # ---- G-B11 output completeness ---------------------------------------------------------
    good = tmp / "good.mzML"
    st_build_mzml(good, [{"mz": 500.0, "z": 2, "im": 1.0, "rt": 60.0, "iso_target": 500.0,
                          "peaks": [(100.0, 5.0)]}], "file:///x.d")
    st.expect_ok("G-B11", "complete mzML validates",
                 lambda: hash_and_validate_mzml(good))
    trunc = tmp / "trunc.mzML"
    st_build_mzml(trunc, [{"mz": 500.0, "z": 2, "im": 1.0, "rt": 60.0, "iso_target": 500.0,
                           "peaks": [(100.0, 5.0)]}], "file:///x.d", truncate=True)
    st.expect_abort("G-B11", "truncated mzML (no closing tag)",
                    lambda: hash_and_validate_mzml(trunc))
    miscount = tmp / "miscount.mzML"
    st_build_mzml(miscount, [{"mz": 500.0, "z": 2, "im": 1.0, "rt": 60.0, "iso_target": 500.0,
                              "peaks": [(100.0, 5.0)]}], "file:///x.d", count_override=7)
    st.expect_abort("G-B11", "spectrumList count != spectra written",
                    lambda: hash_and_validate_mzml(miscount))
    empty = tmp / "empty.mzML"
    empty.write_bytes(b"")
    st.expect_abort("G-B11", "empty output", lambda: hash_and_validate_mzml(empty))

    tsv_ok = tmp / "ok.tsv"
    tsv_ok.write_text("\t".join(SAGE_REQUIRED_COLUMNS) + "\n" +
                      "\t".join(["1", "a", "b", "1", "1", "P", "2", "0.001", "0.001",
                                 "0.001", "9"]) + "\n")
    st.expect_ok("G-B11", "complete TSV validates",
                 lambda: hash_and_validate_tsv(tsv_ok, SAGE_REQUIRED_COLUMNS))
    tsv_nonl = tmp / "nonl.tsv"
    tsv_nonl.write_text(tsv_ok.read_text().rstrip("\n"))
    st.expect_abort("G-B11", "TSV without trailing newline is truncated",
                    lambda: hash_and_validate_tsv(tsv_nonl, SAGE_REQUIRED_COLUMNS))
    tsv_missing = tmp / "missing.tsv"
    cols = [c for c in SAGE_REQUIRED_COLUMNS if c != "peptide_q"]
    tsv_missing.write_text("\t".join(cols) + "\n" + "\t".join(["x"] * len(cols)) + "\n")
    st.expect_abort("G-B11", "TSV missing peptide_q column",
                    lambda: hash_and_validate_tsv(tsv_missing, SAGE_REQUIRED_COLUMNS))
    tsv_ragged = tmp / "ragged.tsv"
    tsv_ragged.write_text("\t".join(SAGE_REQUIRED_COLUMNS) + "\n" +
                          "\t".join(["1"] * (len(SAGE_REQUIRED_COLUMNS) - 2)) + "\n")
    st.expect_abort("G-B11", "TSV row with wrong field count",
                    lambda: hash_and_validate_tsv(tsv_ragged, SAGE_REQUIRED_COLUMNS))

    # ---- G-B14 INI canonicalization --------------------------------------------------------
    ini_a = tmp / "a.ini"
    ini_b = tmp / "b.ini"
    write_resolved_ini(env["tool"], {"trace:ms1_split_valleys": 0.10}, 2, ini_a)
    write_resolved_ini(env["tool"], {"trace:ms1_split_valleys": 0.1}, 2, ini_b)
    _, sha_a = canonical_ini(ini_a)
    _, sha_b = canonical_ini(ini_b)
    st.results.append(("PASS" if sha_a == sha_b else "FAIL", "G-B14",
                       "0.10 and 0.1 hash identically",
                       "%s vs %s" % (short(sha_a), short(sha_b))))
    st.expect_abort("G-B14", "unparseable INI refuses to hash",
                    lambda: canonical_ini(tsv_ok))

    # ---- G-B14 no-op detection (output-digest equality) -------------------------------------
    st.expect_abort("G-B14", "arm byte-identical to its baseline", lambda: detect_noop_arms([
        {"sample": "dataset A", "arm": "base", "baseline_arm": None, "produce_sha256": "aa",
         "ini_sha256": "i1"},
        {"sample": "dataset A", "arm": "split", "baseline_arm": "base", "produce_sha256": "aa",
         "ini_sha256": "i2"}]))
    st.expect_abort("G-B14", "arm resolving to its baseline's config", lambda: detect_noop_arms([
        {"sample": "dataset A", "arm": "base", "baseline_arm": None, "produce_sha256": "aa",
         "ini_sha256": "i1"},
        {"sample": "dataset A", "arm": "split", "baseline_arm": "base", "produce_sha256": "bb",
         "ini_sha256": "i1"}]))
    st.expect_ok("G-B14", "genuinely different arm accepted", lambda: detect_noop_arms([
        {"sample": "dataset A", "arm": "base", "baseline_arm": None, "produce_sha256": "aa",
         "ini_sha256": "i1"},
        {"sample": "dataset A", "arm": "split", "baseline_arm": "base", "produce_sha256": "bb",
         "ini_sha256": "i2"}]))

    # ---- G-B17 plan/runset set equality ----------------------------------------------------
    def rs(cells, samples=("dataset A",), arms=("base", "split")):
        return {"plan_samples": list(samples),
                "plan_arms": [{"name": a} for a in arms], "cells": cells}
    st.expect_abort("G-B17", "missing cell", lambda: assert_cells_match_plan(
        rs([{"sample": "dataset A", "arm": "base"}])))
    st.expect_abort("G-B17", "extra undeclared cell", lambda: assert_cells_match_plan(
        rs([{"sample": "dataset A", "arm": "base"}, {"sample": "dataset A", "arm": "split"},
            {"sample": "dataset D", "arm": "base"}])))
    st.expect_abort("G-B17", "two completed runs for one plan cell",
                    lambda: assert_cells_match_plan(
                        rs([{"sample": "dataset A", "arm": "base"},
                            {"sample": "dataset A", "arm": "base"},
                            {"sample": "dataset A", "arm": "split"}])))
    st.expect_ok("G-B17", "exact match accepted", lambda: assert_cells_match_plan(
        rs([{"sample": "dataset A", "arm": "base"}, {"sample": "dataset A", "arm": "split"}])))

    # ---- G-B12 seal integrity ---------------------------------------------------------------
    def tampered_seal():
        obj = {"seal_version": CONTRACT_VERSION, "status": "ok", "x": 1}
        obj["seal"] = seal_digest(obj)
        obj["x"] = 2
        verify_seal(obj, "<selftest>")
    st.expect_abort("G-B12", "hand-edited SEAL.json", tampered_seal)
    st.expect_ok("G-B12", "intact seal verifies", lambda: verify_seal(
        _sealed({"seal_version": CONTRACT_VERSION, "status": "ok", "x": 1}), "<selftest>"))

    # ---- G-B18 acquisition method ------------------------------------------------------------
    def no_windows():
        d = tmp / "nowin.d"
        d.mkdir(exist_ok=True)
        con = sqlite3.connect(str(d / "analysis.tdf"))
        con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
        con.commit()
        con.close()
        acquisition_identity(d)
    st.expect_abort("G-B18", "no DiaFrameMsMsWindows -> no fallback table", no_windows)
    st.expect_abort("G-B18", "no analysis.tdf at all",
                    lambda: acquisition_identity(tmp / "raw"))
    st.expect_ok("G-B18", "windows extracted per sample",
                 lambda: acquisition_identity(env["d08"]))

    # A different acquisition window table MUST produce a different method_id, or
    # WINDOWS_S23-applied-to-every-sample could recur without anything noticing.
    other_d = st_make_fake_d(tmp / "othermethod.d", mz_centers=(400.0, 700.0),
                             serial="SER-08", sample_name="dataset A")
    mid_a = acquisition_identity(env["d08"])["method_id"]
    mid_b = acquisition_identity(other_d)["method_id"]
    st.results.append(("PASS" if mid_a != mid_b else "FAIL", "G-B18",
                       "different tiles -> different method_id",
                       "%s vs %s" % (mid_a[:8], mid_b[:8])))

    # ---- END-TO-END: the meta-failure ------------------------------------------------------
    _selftest_end_to_end(st, env)
    return st


def _sealed(obj):
    obj = dict(obj)
    obj["seal"] = seal_digest(obj)
    return obj


def _write_plan_raw(env, obj):
    path = env["tmp"] / ("planraw_%s.yaml" % uuid.uuid4().hex[:6])
    _write_yaml(path, obj)
    return path


def _selftest_end_to_end(st, env):
    """The failure no v1 guard addressed, driven end to end.

    Sequence: a good run commits a sealed convert stage; then the SAME recipe is attempted again
    with the converter failing. We assert (1) the second attempt aborts on the return code,
    (2) no second stage was sealed, and (3) attempt A's bytes are byte-for-byte untouched.
    """
    reg = env["reg"]
    tmp = env["tmp"]
    bench = Bench(tmp / "bench")
    plan_path = _plan(env, [BASE_ARM, TREAT_ARM])
    plan = load_plan(plan_path, reg)
    tool_pre = tool_identity(plan["binary"])
    entry, pinned = require_pinned(reg, "dataset A")
    refs = verify_pinned_references(reg, "dataset A")

    os.environ["FAKETOOL_MODE"] = "ok"
    stage_id, seal = do_convert(bench, plan, reg, "dataset A", BASE_ARM, tool_pre, refs, entry, pinned)
    good_dir = bench.stage_dir("convert", stage_id)
    good_digest = sha256_file(good_dir / "work" / "pseudo.mzML")
    st.results.append(("PASS" if good_dir.is_dir() else "FAIL", "G-B12",
                       "good run seals a stage", short(stage_id.split("-")[-1])))

    n_before = len(list((bench.root / "stages" / "convert").iterdir()))

    def failing_rerun():
        os.environ["FAKETOOL_MODE"] = "fail"
        try:
            do_convert(bench, plan, reg, "dataset A", BASE_ARM, tool_pre, refs, entry, pinned)
        finally:
            os.environ["FAKETOOL_MODE"] = "ok"
    st.expect_abort("G-B10", "nonzero converter exit aborts the attempt", failing_rerun)

    n_after = len(list((bench.root / "stages" / "convert").iterdir()))
    st.results.append(("PASS" if n_after == n_before else "FAIL", "G-B12",
                       "failed attempt sealed NOTHING",
                       "%d stages before, %d after" % (n_before, n_after)))
    st.results.append(("PASS" if sha256_file(good_dir / "work" / "pseudo.mzML") == good_digest
                       else "FAIL", "G-B12", "prior stage bytes untouched by the failure",
                       short(good_digest)))
    st.results.append(("PASS" if any((bench.root / "failed").iterdir()) else "FAIL", "G-B10",
                       "failed attempt preserved for forensics",
                       "under %s" % (bench.root / "failed")))

    def truncated_output():
        os.environ["FAKETOOL_MODE"] = "truncate"
        try:
            do_convert(bench, plan, reg, "dataset A", TREAT_ARM, tool_pre, refs, entry, pinned)
        finally:
            os.environ["FAKETOOL_MODE"] = "ok"
    st.expect_abort("G-B11", "exit 0 with a truncated mzML still aborts", truncated_output)

    def miscounted_output():
        os.environ["FAKETOOL_MODE"] = "partial"
        try:
            do_convert(bench, plan, reg, "dataset A", TREAT_ARM, tool_pre, refs, entry, pinned)
        finally:
            os.environ["FAKETOOL_MODE"] = "ok"
    st.expect_abort("G-B11", "exit 0 with a short spectrumList still aborts", miscounted_output)

    def binary_mutated():
        os.environ["FAKETOOL_MODE"] = "mutate_self"
        try:
            do_convert(bench, plan, reg, "dataset A", TREAT_ARM, tool_pre, refs, entry, pinned)
        finally:
            os.environ["FAKETOOL_MODE"] = "ok"
            st_write_faketool(env["tool"])
    st.expect_abort("G-B09", "binary changed across the exec (FAILURE 7)", binary_mutated)

    def raw_retargeted():
        os.environ["FAKETOOL_MODE"] = "retarget"
        try:
            do_convert(bench, plan, reg, "dataset A", TREAT_ARM, tool_pre, refs, entry, pinned)
        finally:
            os.environ["FAKETOOL_MODE"] = "ok"
            st_make_fake_d(env["d08"], mz_centers=(500.0, 600.0), serial="SER-08",
                           sample_name="dataset A")
    st.expect_abort("G-B13", "raw .d altered during the run (symlink retarget class)",
                    raw_retargeted)

    # Re-pin dataset A after the retarget test rewrote its tdf.
    ident = raw_identity(env["d08"])
    reg["samples"]["dataset A"]["pinned"].update({
        "raw_content_id": ident["content_id"], "acquisition_id": ident["acquisition_id"],
        "method_id": ident["method_id"]})
    entry, pinned = require_pinned(reg, "dataset A")

    def sage_fails():
        os.environ["FAKESAGE_MODE"] = "fail"
        try:
            do_search(bench, reg, "dataset A", BASE_ARM, "convert", stage_id,
                      verify_search_toolchain(reg))
        finally:
            os.environ["FAKESAGE_MODE"] = "ok"
    st.expect_abort("G-B10", "sage return code is NOT discarded", sage_fails)

    def sage_drops_column():
        os.environ["FAKESAGE_MODE"] = "drop_peptide_q"
        try:
            do_search(bench, reg, "dataset A", BASE_ARM, "convert", stage_id,
                      verify_search_toolchain(reg))
        finally:
            os.environ["FAKESAGE_MODE"] = "ok"
    st.expect_abort("G-B11", "sage TSV missing peptide_q aborts (not a published zero)",
                    sage_drops_column)

    st.expect_ok("G-B10", "healthy search seals a stage", lambda: do_search(
        bench, reg, "dataset A", BASE_ARM, "convert", stage_id, verify_search_toolchain(reg)))

    st.expect_abort("G-B15", "unknown stage id is not reachable",
                    lambda: bench.load_stage("convert", "convert-deadbeef-nope"))


def cmd_selftest(args):
    tmp = Path(tempfile.mkdtemp(prefix="spextractor-selftest-"))
    keep = bool(args.keep)
    try:
        env = _selftest_env(tmp)
        st = selftest_bench(env)
        rc = st.report()
        if keep:
            print("\nfixtures kept at %s" % tmp)
        return rc
    finally:
        if not keep:
            shutil.rmtree(str(tmp), ignore_errors=True)


# ---------------------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("plan-check", help="validate a plan; no side effects")
    p.add_argument("--plan", required=True)

    p = sub.add_parser("pin", help="write pinned content digests into samples.yaml")
    p.add_argument("--sample", required=True)

    p = sub.add_parser("run", help="execute a plan")
    p.add_argument("--plan", required=True)
    p.add_argument("--bench-root")
    p.add_argument("--need-gb", default="50")

    p = sub.add_parser("verify", help="rehash sealed artifacts; writes nothing")
    p.add_argument("--runset", required=True)
    p.add_argument("--bench-root")
    p.add_argument("--deep", action="store_true")

    p = sub.add_parser("show", help="print a runset")
    p.add_argument("--runset", required=True)
    p.add_argument("--bench-root")

    p = sub.add_parser("selftest", help="deliberately trigger every guard and verify it fires")
    p.add_argument("--keep", action="store_true")

    args = ap.parse_args(argv)
    if not args.cmd:
        ap.print_help()
        return 2
    table = {"plan-check": cmd_plan_check, "pin": cmd_pin, "run": cmd_run,
             "verify": cmd_verify, "show": cmd_show, "selftest": cmd_selftest}
    try:
        return table[args.cmd](args)
    except Abort as exc:
        sys.stderr.write("\nABORT %s\n%s\n\n" % (exc.guard, exc.message))
        return 2


if __name__ == "__main__":
    sys.exit(main())
