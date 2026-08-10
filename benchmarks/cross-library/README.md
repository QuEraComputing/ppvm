# Cross-library Pauli-propagation benchmark

`ppvm-*-2` against its own predecessor and the three other single-threaded
Pauli-propagation engines we know of, on two workloads:

| library | language | entry point used |
|---|---|---|
| **`ppvm-2`** (this repo) | Rust | `Sum<HashMapStore<PauliWord<[u8; N]>, f64>, CoefficientThreshold>` + the `RotationOne`/`RotationTwo` gates |
| **`ppvm-1`** (this repo, pre-refactor) | Rust | `PauliSum<config::fxhash::Byte<N, f64, CoefficientThreshold>>` — the legacy `ppvm-pauli-sum` |
| [PauliPropagation.jl][pp] | Julia | `propagate(PauliRotation(...), psum; min_abs_coeff)` |
| [PauliStrings.jl][ps] | Julia | `trotter_step!(O, gates; truncation, truncate_every)` |
| [pauli-prop][qk] (Qiskit) | Rust-accelerated Python | `propagate_through_circuit(op, qc, max_terms, atol, frame="h")` |

```bash
# Everything, with the term-for-term agreement check first.
uv run --no-project python3 benchmarks/cross-library/run_xbench.py \
    --qubits-tfim 8,16,24,32,40,48,56,64 \
    --qubits-heisenberg 6,8,10,12,14 \
    --steps 10 --atol 1e-6 --iters 2 --out target/xbench

uv run --no-project --with matplotlib python3 benchmarks/cross-library/plot_xbench.py \
    --csv target/xbench/results.csv --out target/xbench/xbench.png \
    --title "Pauli propagation: ppvm-2 vs ppvm-1, PauliPropagation.jl, PauliStrings.jl, pauli-prop"
```

or `mise run xbench`. A subset, if you prefer: `--libs ppvm-1,ppvm-2` is the
v1→v2 comparison on its own.

## What is measured

Both workloads are **Heisenberg-picture propagation of an observable through an
explicit first-order Trotter product of Pauli rotations**, with a
coefficient-magnitude truncation after every gate. Everything below is fixed
across all five engines — the gate list, its order, the angles, the truncation
rule, and the readout.

### `tfim` — transverse-field Ising, magnetization

`H = J Σᵢ ZᵢZᵢ₊₁ + h Σᵢ Xᵢ` on an open chain. One step, in this order:

1. `RX(2h·dt)` on site `0, 1, …, n−1`
2. `RZZ(2J·dt)` on bond `(0,1), (1,2), …, (n−2,n−1)`

Observable `O = Σᵢ Zᵢ`; readout `⟨0…0|O(t)|0…0⟩`, i.e. the sum of the
coefficients of the X-free terms.

### `heisenberg` — isotropic Heisenberg + field, autocorrelator

`H = J Σᵢ (XᵢXᵢ₊₁ + YᵢYᵢ₊₁ + ZᵢZᵢ₊₁) + h Σᵢ Zᵢ`. One step, in this order:

1. `RXX`, `RYY`, `RZZ` at `2J·dt` on bond `(0,1)`, then `(1,2)`, …
2. `RZ(2h·dt)` on site `0, 1, …, n−1`

Observable `O = Z₀`; readout the autocorrelator `S(t) = tr[Z₀·O(t)]/2ⁿ`, which
is just the coefficient of `Z₀` (the Paulis are orthonormal under that pairing).

`θ = 2·c·dt` for a Hamiltonian term `c·G` is the convention all five engines use
for `exp(iθ/2·G)·P·exp(−iθ/2·G)`, so the propagated operators are identical, not
merely similar.

## The parameter contract

Every runner reads the same environment variables and writes the same CSV, so
they can be run directly as well as through the driver:

| variable | meaning |
|---|---|
| `MODEL` | `tfim` or `heisenberg` |
| `QUBITS` | comma-separated widths |
| `STEPS` | Trotter steps |
| `DT`, `JCOUP`, `HFIELD` | `dt`, `J`, `h` |
| `ATOL` | truncation threshold on `|c|` |
| `ITERS` | timed repeats; the **minimum** is reported |
| `DUMP` | print the propagated support instead of timing it |
| `MAX_TERMS` | `pauli-prop` only — its mandatory cap (see below) |

CSV columns: `model,library,qubits,steps,dt,atol,time_s,terms,observable`.

```bash
MODEL=heisenberg QUBITS=8,10,12 STEPS=10 ATOL=1e-6 ITERS=2 \
  cargo run --release -p ppvm-pauli-sum-2 --example xbench
```

## Validation — why the driver refuses to time first

`run_xbench.py` runs every engine with `DUMP=1` at `n=4, steps=3, atol=1e-14`
and diffs the whole propagated support against `ppvm-2`'s, term for term. It
aborts on any missing term, extra term, or coefficient difference above
`--validate-tol` (1e-10 by default).

This is not ceremony. Two real bugs in this harness produced numbers that looked
completely reasonable and were wrong:

* **The `pauli-prop` circuit was reversed.** In the Heisenberg frame it
  conjugates from the *end* of the instruction list backwards, so the spec's
  gate order needs `reversed()` on append. Appending forward propagates a
  different operator — 108 terms against TFIM's 124, 61 against Heisenberg's 64,
  coefficients off by up to 0.1 — while the readout still matched to 9 digits.
* **Duplicate Paulis were silently collapsing.** `propagate_through_circuit`
  returns a `SparsePauliOp` that may list the same Pauli more than once, so
  `len(op)` is a row count, not a support size, and a `{label: coeff}` dict
  comprehension keeps the last duplicate instead of summing them.

An observable is one scalar and can agree by luck or by cancellation. A
term-for-term diff cannot.

## Known differences between the engines

These are real and are the reason the plot shows the **workload size** next to
the runtime. Read them before quoting a ratio.

* **`pauli-prop`'s truncation is a different rule.** Its `atol` documents as
  "terms with coeff magnitudes less than this will *not be added* to the
  operator", i.e. it prunes at branch-creation time, where the other four
  accumulate first and drop the merged coefficient. The other four then agree
  on the support **exactly**, at every width, on both models; `pauli-prop` does
  not. On the sweep below it runs 11–36 % *under* the reference support on TFIM
  (the gap widening with `n`) and 4–15 % *over* it on Heisenberg. The driver
  prints a `note:` line for any engine more than 2 % off, and the plot's third
  panel is there so the ratio is never read without it.
* **`pauli-prop` also has a mandatory `max_terms` cap** with pre-allocation. The
  runner defaults it to `2²²` and **fails the run** if the support ever reaches
  it, rather than quietly reporting a differently-truncated number.
* **PauliStrings.jl stores `im^{#Y}` inside the coefficient** (its `Matrix`
  convention). The dump goes through `op_to_strings`, which puts its
  coefficients on the same real footing as everyone else's.
* **PauliStrings.jl's own front door is `evolve(H, O, tspan; method=Trotter())`,
  which we do not use** — it derives the gate list from the Hamiltonian's
  internal string order, which is not the spec's order. The runner builds the
  `TrotterGate` vector by hand instead.
* **Julia is timed after a warm-up run** so the reported time excludes JIT.
* **Circuit construction is outside the timed region for every engine.**
* Everything is single-threaded (`julia -t1`; `pauli-prop` is single-threaded by
  design; both `ppvm` versions here use no Rayon).

## A run

Darwin, Apple silicon, single-threaded, `--steps 10 --dt 0.1 --j 1 --h 1
--atol 1e-6 --iters 2`; runtime relative to `ppvm-2` (lower is `ppvm-2` winning
by more). `terms` is the shared support the four
accumulate-then-truncate engines all reach exactly.

`tfim`, widths 8…64:

| n | terms | ppvm-2 | ppvm-1 | PauliPropagation.jl | PauliStrings.jl | pauli-prop |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 4 701 | 0.0025 s | 1.17× | 2.52× | 6.25× | 7.06× |
| 24 | 34 353 | 0.0215 s | 1.13× | 6.22× | 8.73× | 6.30× |
| 40 | 64 001 | 0.0543 s | 1.10× | 6.49× | 9.74× | 6.27× |
| 64 | 108 473 | 0.1218 s | 1.11× | 10.38× | 11.44× | 6.72× |

`heisenberg`, widths 6…14:

| n | terms | ppvm-2 | ppvm-1 | PauliPropagation.jl | PauliStrings.jl | pauli-prop |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 16 324 | 0.0382 s | 1.25× | 4.11× | 7.92× | 2.62× |
| 10 | 225 353 | 0.3782 s | 1.25× | 3.80× | 8.46× | 2.05× |
| 12 | 1 174 849 | 1.9868 s | 1.27× | 3.13× | 8.06× | 1.52× |
| 14 | 2 915 879 | 4.3917 s | 1.27× | 2.87× | 7.95× | 1.18× |

### v1 → v2

`ppvm-1` is the same repository before the `traits-2` refactor, run on a
**matched** configuration — `HashMap` + `[u8; N]` + `FxHasher` +
`CoefficientThreshold` on both sides — so the column isolates the engine, not
the storage backend. `-2` is **10–15 % faster on TFIM and 22–27 % faster on
Heisenberg**, and the two agree term-for-term at every validated point.

That is a useful counterweight to `docs/performance-report.md`, whose old-vs-`-2`
matrix is nanobenchmark-heavy and where individual rows are noisy enough that
several needed `-Cllvm-args` layout controls to adjudicate. Here the unit is a
whole workload with a million-term support, and the gap is stable across widths.

Read `pauli-prop`'s column against the workload note above — on TFIM at `n=64`
it is propagating 36 % fewer terms than the other three, and on Heisenberg at
`n=10`–`12` about 15 % more.

## Scaling note

The TFIM support grows roughly linearly in `n` at fixed depth, so that sweep
reaches 64 qubits cheaply. The Heisenberg support grows much faster — three
non-commuting bond rotations per bond per step — and the `Z₀` autocorrelator
itself converges once `n` exceeds the light cone (identical from `n≈10` at
`steps=10`), so widths past that measure scaling rather than new physics.

## Files

- `../../crates/ppvm-pauli-sum-2/examples/xbench.rs` — the `ppvm-2` runner.
- `../../crates/ppvm-pauli-sum/examples/xbench_v1.rs` — the `ppvm-1` runner.
- `../../julia-benchmarks/benches/xbench_pp.jl` — PauliPropagation.jl.
- `../../julia-benchmarks/benches/xbench_ps.jl` — PauliStrings.jl.
- `xbench_qiskit.py` — `pauli-prop`.
- `run_xbench.py` — validation, driving, CSV merge, summary table.
- `plot_xbench.py` — the figure.

[pp]: https://github.com/MSRudolph/PauliPropagation.jl
[ps]: https://github.com/nicolasloizeau/PauliStrings.jl
[qk]: https://github.com/Qiskit/pauli-prop
