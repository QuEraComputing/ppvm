# TFIM Trotter scaling benchmark

Runtime-per-Trotter-run vs qubit count for the ppvm Pauli-propagation backend
under different hashers, alongside [PauliPropagation.jl][pp] as a reference.
This is the harness behind the "storage cliff" investigation: with `fxhash`,
the cached-hash low bits cluster `hashbrown`'s buckets at high fill, so runtime
balloons toward the top of a storage tier and then drops when the next (wider)
tier kicks in. The high-bit fold (and `gxhash`) remove that bump.

The collected CSVs and the rendered plot are **not** checked in — they are
specific to one machine/run. Only the scripts live here.

## Reproduce

All three series use the same circuit (TFIM, h=1, dt=0.1, truncation 1e-6,
depolarizing 1e-4) and the same qubit-count sweep. The bump is a high-fill
effect, so drive the state large with `J=1.0 STEPS=20`.

```bash
mkdir -p /tmp/tfim_sweep

# 1. ppvm: fxhash (no fold = pre-PR), fxhash (folded = this PR), gxhash.
#    gxhash needs AES at compile time.
RUSTFLAGS="-C target-feature=+aes" J=1.0 STEPS=20 \
  QUBITS="8,16,24,32,40,44,48,52,56,60,64,72,80,88,96,104,112,120,122" ITERS=2 \
  cargo run --release -p ppvm-runtime --example trotter_qubit_sweep \
  > /tmp/tfim_sweep/ppvm.csv

# 2. PauliPropagation.jl reference (single-threaded to match ppvm).
cd julia-benchmarks
J=1.0 STEPS=20 \
  QUBITS="8,16,24,32,40,44,48,52,56,60,64,72,80,88,96,104,112,120,122" ITERS=2 \
  julia --project=@. -t1 benches/trotter_sweep.jl > /tmp/tfim_sweep/pp.csv
cd ..

# 3. Plot (log-y).
uv run --with matplotlib python benchmarks/plot_tfim_sweep.py \
  --ppvm /tmp/tfim_sweep/ppvm.csv \
  --pp   /tmp/tfim_sweep/pp.csv \
  --out  /tmp/tfim_sweep/tfim_trotter_scaling.png
```

## Files

- `../crates/ppvm-runtime/examples/trotter_qubit_sweep.rs` — ppvm sweep,
  replicating the Python storage-tier dispatch for `[u8; N]` storage.
- `../julia-benchmarks/benches/trotter_sweep.jl` — PauliPropagation.jl sweep.
- `plot_tfim_sweep.py` — renders the log-y comparison from the two CSVs.

[pp]: https://github.com/MSRudolph/PauliPropagation.jl
