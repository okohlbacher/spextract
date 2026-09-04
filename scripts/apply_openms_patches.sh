#!/bin/bash
# Apply SpeXtract's OpenMS patches to an OpenMS source tree and install the shared calibration header.
#
# WITHOUT these patches the tool builds and runs, but the Bruker .d reader falls back to a two-point
# linear-in-sqrt TOF->m/z chord that is -5..-11 ppm biased on the files we measured -- which costs
# ~6-11% of closed-search peptide identifications (docs/BENCHMARK-MATRIX-2026-09-01.md). The build is
# silent about it; only the run log says which calibration was used, and every emitted mzML carries
# the answer in the spx:mz_calibration userParam.
#
# Usage: scripts/apply_openms_patches.sh /path/to/OpenMS
set -euo pipefail
cd "$(dirname "$0")/.."
TREE=${1:?usage: apply_openms_patches.sh /path/to/OpenMS-source-tree}
[ -f "$TREE/src/openms/source/FORMAT/BrukerTimsFile.cpp" ] || {
  echo "not an OpenMS source tree (no src/openms/source/FORMAT/BrukerTimsFile.cpp): $TREE" >&2; exit 1; }

# the calibration model is header-only and lives in THIS repo; the patch includes it
install -m 0644 src/TdfMzCalibration.h "$TREE/src/openms/include/OpenMS/FORMAT/TdfMzCalibration.h"
echo "installed src/TdfMzCalibration.h -> $TREE/src/openms/include/OpenMS/FORMAT/"

for p in patches/*.patch; do
  if patch -p0 -N --dry-run -d "$TREE" < "$p" >/dev/null 2>&1; then
    patch -p0 -N -d "$TREE" < "$p"; echo "applied  $p"
  elif patch -p0 -R --dry-run -d "$TREE" < "$p" >/dev/null 2>&1; then
    echo "already applied, skipping  $p"
  else
    echo "FAILED to apply $p -- the tree may be a different OpenMS version" >&2; exit 1
  fi
done
echo "done. Rebuild the spextract target, then verify the run log says:"
echo "  TOF m/z calibration: TDF MzCalibration table model (exact, license-free)"
