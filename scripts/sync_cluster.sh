#!/bin/bash
# Deploy src/spextractor.cpp to the cluster OpenMS tree, rebuild, record provenance.
# The ONLY sanctioned deploy path (2026-09-01 review, hardened per kimi F7/F8/F11 + codex #17-22):
#   - refuses staged OR unstaged dirt (git status --porcelain), --allow-dirty to override
#   - remote lock so deploys cannot interleave with each other
#   - old source kept as .bak; restored if the build fails (no silent src/binary drift)
#   - provenance records OpenMS version + EPD-patch marker count, not just our source
#   - every deploy appended to deploys.log on the node
set -euo pipefail
cd "$(dirname "$0")/.."
NODE=${1:-ibminode05}
SRC=src/spextractor.cpp
REMOTE_SRC=/scratch/kohlbach/OpenMS/src/topp/spextractor.cpp
if [ -n "$(git status --porcelain -- "$SRC" src/TdfMzCalibration.h src/TdfLoad.h src/MzPeakStreamLoad.h)" ]; then
  if [ "${2:-}" != "--allow-dirty" ]; then
    echo "REFUSED: $SRC has staged or unstaged changes. Commit first, or pass --allow-dirty." >&2; exit 1
  fi
  DIRTY="+dirty"
else
  DIRTY=""
fi
SHA=$(shasum -a 256 "$SRC" | cut -d' ' -f1)
HDR_SHA=$(shasum -a 256 src/TdfMzCalibration.h | cut -d' ' -f1)
HDRS_SHA=$(cat src/TdfLoad.h src/MzPeakStreamLoad.h | shasum -a 256 | cut -d' ' -f1)
GITREV="$(git rev-parse --short HEAD)${DIRTY}"
scp -q "$SRC" "$NODE":"${REMOTE_SRC}.new"
# the calibration model is header-only and lives in THIS repo; the OpenMS patch includes it, so it
# must travel with every deploy or the cluster build silently uses a stale copy
scp -q src/TdfMzCalibration.h "$NODE":/scratch/kohlbach/OpenMS/src/openms/include/OpenMS/FORMAT/
# The tool's other headers live next to its source in the OpenMS topp directory. They were copied by
# hand once and then drifted: a rename left the cluster building against a stale TdfLoad.h. Sync them.
scp -q src/TdfLoad.h src/MzPeakStreamLoad.h "$NODE":/scratch/kohlbach/OpenMS/src/topp/
# The end-to-end suite lives at the deployment point and is what every verification run executes.
# It drifted once already: a newly added check silently did not run because the copy was stale.
scp -q test/test_spextractor.py "$NODE":/ceph/ibmi/abi/oliver/spextract/tools/
# Refuse to deploy into a tree that has lost the calibration patch: it would build fine, silently
# revert masses to the -5..-11 ppm chord, and record a clean-looking provenance [kimi review].
ssh "$NODE" "grep -q TdfMzCalibration /scratch/kohlbach/OpenMS/src/openms/source/FORMAT/BrukerTimsFile.cpp" || {
  echo "REFUSED: the cluster OpenMS tree has no calibration patch applied." >&2
  echo "  run: scripts/apply_openms_patches.sh /scratch/kohlbach/OpenMS   (on the node)" >&2
  exit 1; }

ssh "$NODE" "set -eu
  mkdir /scratch/kohlbach/.deploy.lock 2>/dev/null || { echo 'REFUSED: another deploy holds the lock' >&2; exit 1; }
  trap 'rmdir /scratch/kohlbach/.deploy.lock' EXIT
  cp -p '$REMOTE_SRC' '${REMOTE_SRC}.bak' 2>/dev/null || true
  mv '${REMOTE_SRC}.new' '$REMOTE_SRC'
  cd /scratch/kohlbach/OpenMS/build
  if ! make spextractor -j32 > /scratch/kohlbach/build_sync_\$(date +%Y%m%d-%H%M%S).log 2>&1; then
    mv '${REMOTE_SRC}.bak' '$REMOTE_SRC' 2>/dev/null || true
    echo 'BUILD FAILED - source restored, binary and provenance unchanged' >&2; exit 2
  fi
  OMSVER=\$(grep -m1 CF_OPENMS_PACKAGE_VERSION_FULLSTRING CMakeCache.txt | cut -d= -f2)
  EPD=/scratch/kohlbach/OpenMS/src/openms/source/FEATUREFINDER/ElutionPeakDetection.cpp
  { printf 'src_sha256=%s\ncalibration_header_sha256=$HDR_SHA\ntool_headers_sha256=$HDRS_SHA\ngit=%s\nbuilt=%s\nbinary_sha256=%s\nopenms=%s\nepd_sha256=%s\nepd_patch_markers=%s\nmasstrace_move_markers=%s\nlibopenms_sha256=%s\n' \
      '$SHA' '$GITREV' \"\$(date -Is)\" \"\$(sha256sum bin/spextractor | cut -d' ' -f1)\" \
      \"\$OMSVER\" \"\$(sha256sum \$EPD | cut -d' ' -f1)\" \"\$(grep -c 'SpeXtractor\|lock' \$EPD || true)\" \
      \"\$(grep -c 'MassTrace(MassTrace &&)' /scratch/kohlbach/OpenMS/src/openms/include/OpenMS/KERNEL/MassTrace.h || true)\" \"\$(sha256sum lib/libOpenMS.so | cut -d' ' -f1)\"
  } > bin/spextractor.provenance.tmp && mv bin/spextractor.provenance.tmp bin/spextractor.provenance
  cat bin/spextractor.provenance | tee -a /scratch/kohlbach/deploys.log"
