#!/usr/bin/env python3
"""Add value="" to every valueless <cvParam> so MSFragger's batmass-io can read an OpenMS mzML.

OpenMS writes boolean-type cvParams with NO value attribute (legal mzML 1.1: e.g.
<cvParam ... name="centroid spectrum"/>). batmass-io (MSFragger's reader) calls String.trim() on
the missing value, NPEs on the spectrum, and silently reports "Scans = 0" -- the file looks empty
and MSFragger searches nothing. This is the reader's bug, not the file's, but the fix lives here:
give every cvParam a value.

Streaming, line-oriented: the 14.8 GB file never lands in memory. Only cvParam TAGS are touched;
the base64 binary arrays (the bulk) pass through untouched. A cvParam that already has value= is
left exactly as-is, so RT/mz/intensity values are never altered.
"""
import re, sys

SRC, DST = sys.argv[1], sys.argv[2]
# match a whole <cvParam ...> or <cvParam .../> tag
TAG = re.compile(rb"<cvParam\b[^>]*?/?>")


def fix(tag):
    if b"value=" in tag:
        return tag                                  # already has a value -> untouched
    # insert value="" just before the closing (/> or >), preserving self-closing form
    if tag.endswith(b"/>"):
        return tag[:-2].rstrip() + b' value=""/>'
    return tag[:-1].rstrip() + b' value=""'  + b'>'


n_lines = n_fixed = 0
with open(SRC, "rb") as fi, open(DST, "wb") as fo:
    for line in fi:
        n_lines += 1
        if b"<cvParam" in line and b"value=" not in line:
            new, k = TAG.subn(lambda m: fix(m.group(0)), line)
            n_fixed += k
            fo.write(new)
        else:
            fo.write(line)
        if n_lines % 20_000_000 == 0:
            sys.stderr.write("  ...%d lines, %d cvParams fixed\n" % (n_lines, n_fixed))
            sys.stderr.flush()
sys.stderr.write("done: %d lines, %d valueless cvParams patched -> %s\n" % (n_lines, n_fixed, DST))
