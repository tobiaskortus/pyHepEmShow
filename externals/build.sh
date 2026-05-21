#!/bin/bash
PREFIX=$1

# Always run relative to the externals/ directory
cd "$(dirname "$0")"

cd g4hepem
cd    build_$PREFIX
make -j install
cd ../..

cd ../hepemshow
cd    build_$PREFIX
make -j
cd ../..
