#!/usr/bin/env bash
# Profile measure_batch across three coefficient-count regimes (few ~1, mid ~100,
# large ~1000 coefficients) with samply, which needs NO sudo on macOS. Saves a
# Firefox-Profiler profile per regime under target/profiles/ and prints a
# top-functions-by-self-time summary via scripts/samply_top.py.
#
# Open a saved profile interactively with:  samply load target/profiles/<name>.json.gz
#
# Run from the repo root:
#     ./scripts/profile_measure_batch.sh
set -euo pipefail

export PATH="$HOME/.cargo/bin:$PATH"
OUT_DIR="${OUT_DIR:-target/profiles}"
SECS="${FLAME_SECS:-6}"
mkdir -p "$OUT_DIR"

cargo build --release -p ppvm-tableau --example profile_measure_batch
BIN="target/release/examples/profile_measure_batch"

for w in few mid large; do
  echo "=== samply: measure_batch / $w ==="
  FLAME_SECS="$SECS" samply record --save-only --unstable-presymbolicate \
    -o "$OUT_DIR/measure_batch_${w}.json.gz" -- "$BIN" "$w" flame
done

echo
for w in few mid large; do
  python3 scripts/samply_top.py "$OUT_DIR/measure_batch_${w}.json.gz" 20
  echo
done
