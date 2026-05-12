#!/usr/bin/env bash
set -euo pipefail

# Quick helper: configure, build with small demo parameters, then run common tests
# Usage: ./run-small.sh

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="$ROOT_DIR/build-small"

# Small demo params (fast): n=8, q=257, b=2, lambda=8, N_id=3, d=1
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
    -DUNIFIED_PARAMS=ON

cmake --build "$BUILD_DIR" -j

pushd "$BUILD_DIR" > /dev/null

# Run all executables in this build directory (skip common CMake files)
shopt -s nullglob
EXCLUDE_NAMES=("CMakeFiles" "CMakeCache.txt" "cmake_install.cmake" "Makefile" "build.ninja" "rules.ninja" "CTestTestfile.cmake" "cmake" "cpack" "ctest")

files=(./*)
# Sort by name for deterministic order
IFS=$'\n' files_sorted=($(sort <<<"${files[*]}"))
for path in "${files_sorted[@]}"; do
    exe=$(basename "$path")
    # only regular files and executable
    if [ ! -f "$path" ] || [ ! -x "$path" ]; then
        continue
    fi
    skip=false
    for ex in "${EXCLUDE_NAMES[@]}"; do
        if [ "$exe" = "$ex" ]; then skip=true; break; fi
    done
    if $skip; then continue; fi

    echo
    echo "╔══════════════════════════════════════════╗"
    printf "║  %-40s║\n" "$exe"
    echo "╚══════════════════════════════════════════╝"

    ./$exe || { echo "$exe failed"; exit 1; }
done
shopt -u nullglob

popd > /dev/null

echo "All small-parameter tests finished."