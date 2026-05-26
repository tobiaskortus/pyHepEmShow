#!/bin/bash
PREFIX=$1
FLAGS="$2"
PYTHON_EXE="${3:-python}"   # optional 3rd arg: path to Python executable

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
 -DCMAKE_INSTALL_PREFIX=$PWD/../install_$PREFIX \
 -DCMAKE_BUILD_TYPE=Debug -DCODI_FORWARD=yes \
 -DCoDiPack_DIR=$PWD/../../CoDiPack/cmake \
 -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3 -ggdb3 -fno-omit-frame-pointer -fno-inline" \
 -DCMAKE_CXX_FLAGS="$FLAGS" \
 $RPATH_FLAGS

make -j install
cd ../..

cd ../hepemshow
mkdir -p build_$PREFIX
cd       build_$PREFIX
cmake .. \
  -DG4HepEm_DIR=$PWD/../../externals/g4hepem/install_$PREFIX/$LIB_DIR/cmake/G4HepEm/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCODI_FORWARD=yes \
  -DCoDiPack_DIR=$PWD/../../externals/CoDiPack/cmake \
  -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3 -ggdb3 -fno-omit-frame-pointer -fno-inline" \
  -DCMAKE_CXX_FLAGS="$FLAGS" \
  -DPython_EXECUTABLE="$PYTHON_EXE" \
  $RPATH_FLAGS

make -j
cd ../..
