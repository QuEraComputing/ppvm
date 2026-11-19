#!/usr/bin/env bash
# Fetch clifft-bench workload circuits into target/clifft-bench/ (not committed).
set -euo pipefail
dest="$(git rev-parse --show-toplevel)/target/clifft-bench"
mkdir -p "$dest"
base=https://raw.githubusercontent.com/unitaryfoundation/clifft-bench/main/workloads/circuits
for name in "${@:-pure_surface_d7_r7_p1e-3}"; do
  curl -fsSL "$base/$name.stim" -o "$dest/$name.stim"
  echo "$dest/$name.stim"
done
