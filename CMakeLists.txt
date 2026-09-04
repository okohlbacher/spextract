cmake_minimum_required(VERSION 3.16)
project(SpeXtract VERSION 1.1.0 LANGUAGES CXX)

# SpeXtract is a STANDALONE application that links against an existing OpenMS
# installation. It is deliberately not part of the OpenMS source tree: OpenMS is a
# prerequisite here, not a host.
#
# Point CMake at an OpenMS build or install with:
#   cmake -DOpenMS_DIR=/path/to/OpenMS/build ..

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# The tool itself cannot be built without a PATCHED OpenMS: it includes
# <OpenMS/FORMAT/TdfMzCalibration.h>, which scripts/apply_openms_patches.sh installs into the
# OpenMS tree, and calls BrukerTimsFile::lastMzCalibration, which the same patch adds. No released
# OpenMS package contains either. The calibration model, however, is a self-contained header, so
# it can be compiled and tested anywhere -- which is what CI does on every platform.
option(SPEXTRACT_TESTS_ONLY "Build only the OpenMS-independent tests (portability CI)" OFF)
if(SPEXTRACT_TESTS_ONLY)
  add_executable(test_calibration_cpp tests/test_calibration_cpp.cpp)
  target_include_directories(test_calibration_cpp PRIVATE src)
  enable_testing()
  add_test(NAME calibration_golden
           COMMAND test_calibration_cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests/calibration_golden.json)
  return()
endif()

find_package(OpenMS REQUIRED)
find_package(OpenMP)

add_executable(spextract src/spextract.cpp)
target_link_libraries(spextract PRIVATE OpenMS)

if(OpenMP_CXX_FOUND)
  # The window/band loop is OpenMP-parallel; without it the tool still builds and runs,
  # but serially -- and the memory-admission logic degrades to one window at a time.
  target_link_libraries(spextract PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found: spextract will run single-threaded.")
endif()

install(TARGETS spextract RUNTIME DESTINATION bin)

# Optional streaming .mzpeak input through the mzPeak C++ library (github.com/OpenMS/mzpeak).
# -DMZPEAK_ROOT=<checkout with include/ and a built libmzpeak in build*/>. Without it the tool
# still builds and reads .mzpeak only through OpenMS' non-streaming MzPeakFile.
if(MZPEAK_ROOT)
  find_library(MZPEAK_LIB NAMES mzpeak PATHS ${MZPEAK_ROOT} PATH_SUFFIXES build build-release buildfix REQUIRED)
  find_package(Arrow REQUIRED)
  find_package(Parquet REQUIRED)
  target_compile_definitions(spextract PRIVATE SPEXTRACT_WITH_MZPEAK)
  target_include_directories(spextract PRIVATE ${MZPEAK_ROOT}/include)
  target_link_libraries(spextract PRIVATE ${MZPEAK_LIB} Parquet::parquet_static Arrow::arrow_static zip)
endif()
