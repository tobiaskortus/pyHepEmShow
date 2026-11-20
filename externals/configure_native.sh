#!/bin/bash
set -euo pipefail

AD_MODE="${1:-}"
DSBO_INCLUDE_DIR="${2:-}"
HDF5_ROOT="${3:-}"
BUILD_TYPE="${4:-Release}"

if [[ -z "$AD_MODE" || -z "$DSBO_INCLUDE_DIR" ]]; then
  echo "Usage: $0 <forward|backward> <dsbo-include-dir> [hdf5-root] [build-type]"
  exit 1
fi

DSBO_INCLUDE_DIR="$(cd "$DSBO_INCLUDE_DIR" && pwd)"
if [[ -n "$HDF5_ROOT" ]]; then
  HDF5_ROOT="$(cd "$HDF5_ROOT" && pwd)"
fi

case "$AD_MODE" in
  forward)
    PREFIX="forward"
    G4HEPEM_AD_FLAGS=""
    ;;
  backward)
    PREFIX="backward"
    G4HEPEM_AD_FLAGS="-DCODI_REVERSE=yes"
    ;;
  *)
    echo "Unknown mode '$AD_MODE'. Expected: forward, backward"
    exit 1
    ;;
esac

# Always run relative to the externals/ directory
cd "$(dirname "$0")"

if [[ ! -d "g4hepem/install_$PREFIX" ]]; then
  cd g4hepem
  mkdir -p "build_$PREFIX"
  cd "build_$PREFIX"
  cmake .. \
   -DG4HepEm_GEANT4_BUILD=OFF \
   -DBUILD_SHARED_LIBS=OFF \
   -DBUILD_STATIC_LIBS=ON \
   -DCMAKE_INSTALL_PREFIX="$PWD/../install_$PREFIX" \
   -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
   $G4HEPEM_AD_FLAGS \
   -DCoDiPack_DIR="$PWD/../../CoDiPack/cmake" \
   -DCMAKE_CXX_FLAGS="-include algorithm"
  cmake --build . --parallel --target install
  cd ../..
fi

# macOS installs to lib/, Linux typically to lib64/, but some platforms
# (e.g. Ubuntu arm64) install to lib/. Resolve from the actual install tree.
if [[ -d "g4hepem/install_$PREFIX/lib/cmake/G4HepEm" ]]; then
  LIB_DIR="lib"
else
  LIB_DIR="lib64"
fi

cd ../hepemshow
mkdir -p "build_native_$PREFIX"
cd "build_native_$PREFIX"
rm -f CMakeCache.txt

G4HEPEM_LIB_DIR="$(cd ../.. && pwd)/externals/g4hepem/install_$PREFIX/$LIB_DIR"
RPATH_PATHS="$G4HEPEM_LIB_DIR"
if [[ -n "$HDF5_ROOT" ]]; then
  RPATH_PATHS="$RPATH_PATHS;$HDF5_ROOT/lib"
fi

CMAKE_EXTRA_ARGS=()
if [[ -n "$HDF5_ROOT" ]]; then
  CMAKE_EXTRA_ARGS+=("-DHDF5_ROOT=$HDF5_ROOT")
fi
CMAKE_EXTRA_ARGS+=("-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON" "-DCMAKE_INSTALL_RPATH=$RPATH_PATHS")

cmake .. \
  -DG4HepEm_DIR="$PWD/../../externals/g4hepem/install_$PREFIX/$LIB_DIR/cmake/G4HepEm/" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DPYHEPEMSHOW_BUILD_PYTHON=OFF \
  -DPYHEPEMSHOW_BUILD_NATIVE=ON \
  -DPYHEPEMSHOW_DSBO_INCLUDE_DIR="$DSBO_INCLUDE_DIR" \
  -DCoDiPack_DIR="$PWD/../../externals/CoDiPack/cmake" \
  -DCMAKE_CXX_FLAGS="-include algorithm" \
  "${CMAKE_EXTRA_ARGS[@]}"

cmake --build . --parallel --target HepEmShowNative
cd ../..
