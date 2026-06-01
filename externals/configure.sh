#!/bin/bash
set -euo pipefail

PREFIX=$1
FLAGS="$2"
PYTHON_EXE="${3:-python}"       # optional 3rd arg: path to Python executable
BUILD_TYPE="${4:-Release}"      # optional 4th arg: CMake build type (Debug/Release/…)
AD_MODE="${5:-forward_ad}"      # optional 5th arg: one of forward, forward_ad, backward
MODULE_NAME="${6:-_core}"       # optional 6th arg: nanobind module target name

if [[ -z "$PREFIX" ]]; then
  echo "Usage: $0 <prefix> [flags] [python_exe] [build_type] [ad_mode] [module_name]"
  exit 1
fi

case "$AD_MODE" in
  forward)
    G4HEPEM_AD_FLAGS=""
    ;;
  forward_ad)
    G4HEPEM_AD_FLAGS="-DCODI_FORWARD=yes"
    ;;
  backward)
    G4HEPEM_AD_FLAGS="-DCODI_REVERSE=yes"
    ;;
  *)
    echo "Unknown ad_mode '$AD_MODE'. Expected: forward, forward_ad, backward"
    exit 1
    ;;
esac

# Always run relative to the externals/ directory
cd "$(dirname "$0")"

# macOS installs to lib/, Linux typically to lib64/
if [[ "$(uname -s)" == "Darwin" ]]; then
  LIB_DIR="lib"
  # On macOS, build with the install RPATH from the start so make install is
  # idempotent (avoids install_name_tool -delete_rpath failures on re-install)
  RPATH_FLAGS="-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_INSTALL_RPATH=@loader_path"
else
  LIB_DIR="lib64"
  # On Linux, $ORIGIN means "same directory as the loading file"
  RPATH_FLAGS="-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON -DCMAKE_INSTALL_RPATH=\$ORIGIN"
fi

cd g4hepem
mkdir -p build_$PREFIX
cd       build_$PREFIX
cmake .. \
 -DG4HepEm_GEANT4_BUILD=OFF \
 -DBUILD_SHARED_LIBS=OFF \
 -DBUILD_STATIC_LIBS=ON \
 -DCMAKE_INSTALL_PREFIX=$PWD/../install_$PREFIX \
 -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
 $G4HEPEM_AD_FLAGS \
 -DCoDiPack_DIR=$PWD/../../CoDiPack/cmake \
 -DCMAKE_CXX_FLAGS="-include algorithm $FLAGS" \
 $RPATH_FLAGS

make -j install
cd ../..

cd ../hepemshow
mkdir -p build_$PREFIX
cd       build_$PREFIX
# Remove stale cache so Python_EXECUTABLE and other hints take effect cleanly
rm -f CMakeCache.txt
cmake .. \
  -DG4HepEm_DIR=$PWD/../../externals/g4hepem/install_$PREFIX/$LIB_DIR/cmake/G4HepEm/ \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DPYHEPEMSHOW_MODULE_NAME=$MODULE_NAME \
  -DCoDiPack_DIR=$PWD/../../externals/CoDiPack/cmake \
  -DCMAKE_CXX_FLAGS="-include algorithm $FLAGS" \
  -DPython_EXECUTABLE="$PYTHON_EXE" \
  $RPATH_FLAGS

make -j
cd ../..
