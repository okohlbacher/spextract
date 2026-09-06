#!/usr/bin/env python3
"""Diff everything OUTSIDE the spectrum list of two mzML files, with the completion-time stamp masked.

bench/semantic_digest.py deliberately hashes only <spectrumList>..</spectrumList>; this is the
complementary check for a change to the WRITER (run-level userParams, dataProcessing, instrument and
software blocks) that the digest cannot see. The trailing <indexList> is skipped: it holds byte
offsets that move with any header-length change. Exit 1 on any difference.
"""
import difflib, re, sys
STAMP = re.compile(rb'(name="completion time" value=")[^"]*(")')
def header(path):
    with open(path, "rb") as f:
        head = f.read(1 << 24)                       # the header is far below 16 MiB
    i = head.find(b"<spectrumList")
    if i < 0: raise SystemExit(f"{path}: no <spectrumList found in the first 16 MiB")
    return STAMP.sub(rb'\1<masked>\2', head[:i]).decode("utf-8", "replace").splitlines()
if __name__ == "__main__":
    if len(sys.argv) != 3: raise SystemExit("usage: mzml_header_diff.py A.mzML B.mzML")
    a, b = header(sys.argv[1]), header(sys.argv[2])
    d = list(difflib.unified_diff(a, b, sys.argv[1], sys.argv[2], lineterm=""))
    print("\n".join(d) if d else "HEADERS IDENTICAL (completion time masked)")
    sys.exit(1 if d else 0)
