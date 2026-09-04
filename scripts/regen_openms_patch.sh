#!/bin/bash
# Regenerate patches/openms-brukertims-mz-calibration.patch from the cluster's OpenMS tree (NOT a git checkout).
# Every edited file keeps a pristine copy next to it (<file>.orig); this diffs each against its .orig (p0 paths)
# and appends the committed header hunk for files that have no .orig (the BrukerTimsFile.h accessor).
# Usage (on the cluster): scripts/regen_openms_patch.sh /path/to/scratch/OpenMS > new.patch; then grep-verify:
#   par-load, ms1-par, perf-load, TdfMzCalibration x11, lastMzCalibration, load_err, ms1_err, pick-mz
set -eu
TREE=${1:-/path/to/scratch/OpenMS}
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd "$TREE"
for o in $(find src -name "*.orig" | sort); do
  f=${o%.orig}
  diff -u --label "$f" --label "$f" "$o" "$f" || true
done
# files edited without a .orig: keep their hunks from the committed patch
for f in $(grep -oE "^--- [^ ]+" "$HERE/patches/openms-brukertims-mz-calibration.patch" | cut -c5- | sort -u); do
  [ -f "$f.orig" ] && continue
  awk -v f="$f" 'BEGIN{p=0} /^--- /{p=($2==f)} p' "$HERE/patches/openms-brukertims-mz-calibration.patch"
done
