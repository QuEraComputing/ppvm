# ppvm-cupauliprop

GPU backend for ppvm's Heisenberg-picture Pauli propagation on NVIDIA
[cuPauliProp](https://docs.nvidia.com/cuda/cuquantum/latest/cupauliprop/)
(cuQuantum), plus CPU and GPU Trotter benchmarks.

- `CudaPauliSum` (feature `cuda`): device-resident Pauli sum with ppvm gate
  names (`rx`, `rzz`, `rxx`, `ryy`, `pauli_error`, `overlap_with_zero`, `len`).
- `examples/trotter_cpu.rs`: ppvm `PauliSum`, configured like Python's
  `ppvm.PauliSum`.
- `examples/trotter_gpu.rs`: the same circuits on `CudaPauliSum`.

The circuits are copies of the TFIM and Heisenberg benchmarks in
`ppvm-benchmarks/trotter-benchmarks`, and the defaults reproduce that repo's
scan. Results JSON uses its `results*.json` schema.

## Setup (Linux + NVIDIA GPU)

Needs the CUDA 12 toolkit (`libcudart`) and cuPauliProp. The library ships in
the `cupauliprop-cu12` wheel on PyPI (a zip with `cuquantum/lib` and
`cuquantum/include`):

```bash
curl -fsSL https://pypi.org/pypi/cupauliprop-cu12/0.5.0/json \
  | jq -r --arg a "$(uname -m)" '.urls[] | select(.filename | endswith($a + ".whl")) | .url' \
  | xargs curl -fsSLO
unzip cupauliprop_cu12-*.whl 'cuquantum/*' -d ~/cupauliprop
export CUPAULIPROP_LIB_DIR=~/cupauliprop/cuquantum/lib
export CUDA_HOME=/usr/local/cuda            # or CUDART_LIB_DIR=<dir of libcudart.so>
export LD_LIBRARY_PATH=$CUPAULIPROP_LIB_DIR:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

On aarch64 hosts (e.g. GH200) gxhash also needs AES enabled explicitly
(`.cargo/config.toml` only does this for x86_64):

```bash
export RUSTFLAGS="-C target-feature=+aes,+neon"
```

## Running

Check the GPU is not shared (`nvidia-smi --query-compute-apps=pid,name
--format=csv` lists nothing, utilization ~0%), then:

```bash
cargo test --release -p ppvm-cupauliprop --features cuda -- --include-ignored

export RAYON_NUM_THREADS=1
for model in tfim heisenberg; do
  cargo run --release -p ppvm-cupauliprop --example trotter_cpu -- \
      --model $model --out results_cpu_$model.json
  cargo run --release -p ppvm-cupauliprop --features cuda --example trotter_gpu -- \
      --model $model --out results_gpu_$model.json
done
```

Scale up with the CLI, e.g. `--sizes 128,256,512 --cutoff 1e-7`. See
`--help` for all options (`--steps`, `--dt`, `--noise`, `--rounds`, …).
The GPU binary splits 90% of free device memory between the workspace and
its two expansion buffers by default; `--capacity <terms>` and
`--workspace-mib` override that. Running out of either fails with a message
naming the required size.

## Docker

A Docker image with all dependencies (CUDA toolkit, cuPauliProp, Rust) and a
script that runs the tests and both benchmarks against a ppvm checkout live
in ppvm-benchmarks under `trotter-benchmarks/cupauliprop/`.
