#!/usr/bin/env bash
# Build ppvm-cupauliprop with `--features cuda` in Docker and run its tests.
#
#   crates/ppvm-cupauliprop/docker/check.sh          # no GPU: build, link, host-only tests
#   GPU=1 crates/ppvm-cupauliprop/docker/check.sh    # --gpus all: also the GPU tests
#
# The repo is mounted at /work; the Linux target dir and cargo registry live
# in Docker volumes so they don't clash with the host build.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
IMAGE="${IMAGE:-ppvm-cupauliprop}"

docker build -t "$IMAGE" "$HERE"

run_args=(--rm -v "$REPO:/work" -v ppvm-cupp-target:/work/target -v ppvm-cupp-cargo:/opt/cargo/registry)
if [[ "${GPU:-0}" == 1 ]]; then
  run_args+=(--gpus all)
  test_args=(-- --include-ignored)
  extra_ld=""
else
  test_args=()
  extra_ld=":/opt/cuda-stub"
fi

docker run "${run_args[@]}" -e LD_LIBRARY_PATH="/opt/cupauliprop/lib:/usr/local/cuda/lib64$extra_ld" "$IMAGE" bash -euxc "
  cargo build --release -p ppvm-cupauliprop --features cuda --examples
  cargo test --release -p ppvm-cupauliprop -p ppvm-cupauliprop-sys --features cuda ${test_args[*]:-}
  ldd target/release/examples/trotter_gpu | grep -E 'cupauliprop|cudart|libcuda'
"
