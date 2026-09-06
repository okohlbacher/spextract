# Append the mzPeak link block for spextractor to the OpenMS tree's TOPP CMakeLists (cluster).
# Idempotent. Paths: fork lib at /path/to/scratch/mzpeak-cpp/build, Boost.JSON 1.89 from the ceph env,
# Arrow/Parquet/libzip from contrib (what libOpenMS itself links).
p = "/path/to/scratch/OpenMS/src/topp/CMakeLists.txt"; s = open(p).read()
if "SPEXTRACTOR_WITH_MZPEAK" in s: print("already applied"); raise SystemExit
s += '''
# [SpeXtractor mzpeak] streaming .mzpeak input through the mzPeak C++ library (fork trunk).
if(TARGET spextractor AND EXISTS /path/to/scratch/mzpeak-cpp/build/libmzpeak.a)
  set(_spx_env /path/to/shared/AI/OpenDIAlyzer/opt/env)
  target_compile_definitions(spextractor PRIVATE SPEXTRACTOR_WITH_MZPEAK)
  target_compile_features(spextractor PRIVATE cxx_std_23)
  target_include_directories(spextractor PRIVATE /path/to/scratch/mzpeak-cpp/include /path/to/scratch/mzpeak-cpp/boost189inc ${CMAKE_SOURCE_DIR}/src/openms/extern/SQLiteCpp/include ${CMAKE_SOURCE_DIR}/src/openms/extern/SQLiteCpp/sqlite3)
  file(GLOB _spx_numpress /path/to/scratch/mzpeak-cpp/build/subprojects/msnumpress/*.a)
  target_link_libraries(spextractor
    /path/to/scratch/mzpeak-cpp/build/libmzpeak.a ${_spx_numpress}
    ${_spx_env}/lib/libboost_json.so
    /path/to/scratch/contrib/build/lib/libparquet.a /path/to/scratch/contrib/build/lib/libarrow.a
    /path/to/scratch/contrib/build/lib/libarrow_bundled_dependencies.a
    /path/to/scratch/contrib/build/lib/libzip.so /path/to/scratch/contrib/build/lib/libz.so /path/to/scratch/contrib/build/lib/libbz2.a -Wl,--exclude-libs,libarrow_bundled_dependencies.a)
  set_target_properties(spextractor PROPERTIES BUILD_RPATH "${_spx_env}/lib;/path/to/scratch/contrib/build/lib")
endif()
'''
open(p, "w").write(s); print("topp CMake block appended")
