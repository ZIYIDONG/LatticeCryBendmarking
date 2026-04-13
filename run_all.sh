#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/build"

TARGETS=(
    LatticeCryBenchmarking
    test_powersof
    test_powersof_modswitch
    test_frd
    debug_frd
    bench_matops
)

for exe in "${TARGETS[@]}"; do
    echo
    echo "╔══════════════════════════════════════════╗"
    printf  "║  %-40s║\n" "$exe"
    echo "╚══════════════════════════════════════════╝"
    ./"$exe"
done

echo
echo "All tests finished."
