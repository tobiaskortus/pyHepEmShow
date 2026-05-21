#!/bin/bash
PREFIX=$1
FLAGS="$2"

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
  RPATH_FLAGS=""
fi

cd g4hepem
mkdir build_$PREFIX
cd    build_$PREFIX
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
mkdir build_$PREFIX
cd    build_$PREFIX
cmake .. \
  -DG4HepEm_DIR=$PWD/../../externals/g4hepem/install_$PREFIX/$LIB_DIR/cmake/G4HepEm/ \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCODI_FORWARD=yes \
  -DCoDiPack_DIR=$PWD/../../externals/CoDiPack/cmake \
  -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3 -ggdb3 -fno-omit-frame-pointer -fno-inline" \
  -DCMAKE_CXX_FLAGS="$FLAGS"

make -j
cd ../..
