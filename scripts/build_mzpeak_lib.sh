#!/bin/bash
# Build the mzPeak C++ library (github.com/OpenMS/mzpeak; fork okohlbacher/mzpeak-openms trunk) for
# SpeXtract's streaming .mzpeak input, against the SAME Arrow/Parquet/libzip OpenMS was built with.
#   MZPEAK_SRC   checkout of the library                     (default /path/to/scratch/mzpeak-cpp)
#   PKGCONF      pkg-config dir with arrow/parquet/libzip .pc (default OpenMS contrib)
#   BOOST_ROOT   Boost >= 1.89 (json, compat)                 (default: OpenDIAlyzer env on ceph)
#   TOOLS        dir with meson + ninja on PATH               (default: same env)
# GCC 13 lacks std::ranges::to (GCC 14 has it): the four `| std::ranges::to<C>()` sites are rewritten to
# `| spx_to<C>()` and patches/spx_ranges_to.h is force-included. The library pins arrow>=24; OpenMS
# contrib ships 23.0.0 and the code compiles against it, so the pin is relaxed here.
set -eu
MZPEAK_SRC=${MZPEAK_SRC:-/path/to/scratch/mzpeak-cpp}
PKGCONF=${PKGCONF:-/path/to/scratch/contrib/build/lib/pkgconfig}
BOOST_ROOT=${BOOST_ROOT:-/path/to/shared/AI/OpenDIAlyzer/opt/env}
TOOLS=${TOOLS:-$BOOST_ROOT/bin}
HERE=$(cd "$(dirname "$0")/.." && pwd)
export PATH=$TOOLS:$PATH PKG_CONFIG_PATH=$PKGCONF BOOST_ROOT
cd "$MZPEAK_SRC"
cp "$HERE/patches/spx_ranges_to.h" .
sed -i "s/std::ranges::to</spx_to</g" src/wavelength_spectra.cpp src/data/array_index.cpp src/chromatograms.cpp src/spectra.cpp
sed -i "s/arrow_ver_str = '>=24.0.0'/arrow_ver_str = '>=23.0.0'/" meson.build
python3 "$HERE/patches/mzpeak-cpp-arrow23.py" "$MZPEAK_SRC/src/util/metadata_model.cpp"   # one Result-API site
rm -rf build
meson setup build --buildtype=release -Dcpp_args="-I$BOOST_ROOT/include -include $MZPEAK_SRC/spx_ranges_to.h" -Dcpp_link_args="-L$BOOST_ROOT/lib -Wl,-rpath,$BOOST_ROOT/lib"
ninja -C build -j${JOBS:-16}
ls -la build/libmzpeak.a build/subprojects/msnumpress/*.a
