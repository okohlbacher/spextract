#!/usr/bin/env python3
"""Semantic digest of an mzML: sha256 of the DATA (from <spectrumList onwards), ignoring the header.

Two runs of a perfectly deterministic tool can never be byte-identical, for two reasons that have
now each cost this project a false "determinism FAIL":
  1. OpenMS stamps MS:1000747 "completion time" from the wall clock;
  2. the header records the invocation's own parameters -- so a threads-8 run and a threads-100 run
     differ in the header by construction, which is precisely the comparison a thread-scaling
     determinism test wants to make.
Hashing from <spectrumList onwards compares what the science depends on: the spectra themselves.
(This also retroactively explains the July r1-vs-imw005 md5 mismatch recorded as unexplained.)
"""
import hashlib, sys
MARK = b"<spectrumList"
END = b"</spectrumList>"
def digest(path):
    h = hashlib.sha256()
    started = False
    ended = False
    tail = b""
    pending = b""
    with open(path, "rb") as f:
        while True:
            c = f.read(1 << 22)
            if not c:
                break
            if not started:
                buf = tail + c
                i = buf.find(MARK)
                if i < 0:
                    tail = buf[-len(MARK):]
                    continue
                started = True
                c = buf[i:]
            # stop at </spectrumList>: the trailing <indexList> holds byte OFFSETS, which shift with any
            # header-length change (a 4-byte longer <software> version stamp moved every offset and
            # falsely flagged two spectrum-identical files as different, 2026-09-02).
            # The tag can straddle a chunk boundary: search the carry + chunk, hash everything but a
            # tag-length carry, and keep that carry for the next chunk (review finding, 2026-09-06).
            pending += c
            j = pending.find(END)
            if j >= 0:
                h.update(pending[:j + len(END)])
                ended = True
                break
            keep = len(END) - 1
            h.update(pending[:-keep])
            pending = pending[-keep:]
    if not started:
        raise SystemExit(f"{path}: no {MARK!r} found")
    if not ended:
        raise SystemExit(f"{path}: {MARK!r} without {END!r}: truncated or not an mzML spectrum list")
    return h.hexdigest()
if __name__ == "__main__":
    ds = [(p, digest(p)) for p in sys.argv[1:]]
    for p, d in ds:
        print(f"{d}  {p}")
    if len(ds) > 1:
        same = len({d for _, d in ds}) == 1
        print("SPECTRUM DATA " + ("IDENTICAL -> deterministic" if same else "DIFFERS -> NOT deterministic"))
        sys.exit(0 if same else 1)
