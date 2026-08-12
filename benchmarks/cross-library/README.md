# Cross-library Pauli-propagation benchmark

`ppvm` against the four other single-threaded Pauli-propagation engines we know
of, on two workloads:

| library | language | entry point used |
|---|---|---|
| **`ppvm`** (this repo) | Rust | `PauliSum<config::fxhash::Byte<N, f64, CoefficientThreshold>>` |
| [PauliPropagation.jl][pp] | Julia | `propagate(PauliRotation(...), psum; min_abs_coeff)` |
| [PauliStrings.jl][ps] | Julia | `trotter_step!(O, gates; truncation, truncate_every)` |
| [pauli-prop][qk] (Qiskit) | Rust-accelerated Python | `propagate_through_circuit(op, qc, max_terms, atol, frame="h")` |
| [monoprop][mp] (Algorithmiq) | C++ with Python bindings | `PauliPropagator.from_circuit(circuit, op, cutoff, lower_atol)` |

```bash
# Everything, with the term-for-term agreement check first.
uv run --no-project python3 benchmarks/cross-library/run_xbench.py \
    --qubits-tfim 8,16,24,32,40,48,56,64 \
    --qubits-heisenberg 6,8,10,12,14 \
    --steps 10 --atol 1e-6 --iters 2 --out target/xbench

uv run --no-project --with matplotlib python3 benchmarks/cross-library/plot_xbench.py \
    --csv target/xbench/results.csv --out target/xbench/xbench.png \
    --title "Pauli propagation: ppvm vs PauliPropagation.jl, PauliStrings.jl, pauli-prop, monoprop"
```

A subset, if you prefer: `--libs ppvm,monoprop`.

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

`θ = 2·c·dt` for a Hamiltonian term `c·G` is the convention four of the five
engines use directly for `exp(iθ/2·G)·P·exp(−iθ/2·G)`, so the propagated
operators are identical, not merely similar. monoprop is the exception and needs
a conversion; see below.

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
| `monoprop_NUM_THREADS`, `monoprop_PARTITIONS` | `monoprop` only — the thread cap (see below) |

CSV columns: `model,library,qubits,steps,dt,atol,time_s,terms,observable`.

```bash
MODEL=heisenberg QUBITS=8,10,12 STEPS=10 ATOL=1e-6 ITERS=2 \
  cargo run --release -p ppvm-pauli-sum --example xbench
```

## Validation — why the driver refuses to time first

`run_xbench.py` runs every engine with `DUMP=1` at `n=4, steps=3, atol=1e-14`
and diffs the whole propagated support against `ppvm`'s, term for term. It
aborts on any missing term, extra term, or coefficient difference above
`--validate-tol` (1e-10 by default).

This is not ceremony. Three real bugs in this harness produced numbers that
looked completely reasonable and were wrong:

* **The `pauli-prop` circuit was reversed.** In the Heisenberg frame it
  conjugates from the *end* of the instruction list backwards, so the spec's
  gate order needs `reversed()` on append. Appending forward propagates a
  different operator — 108 terms against TFIM's 124, 61 against Heisenberg's 64,
  coefficients off by up to 0.1 — while the readout still matched to 9 digits.
* **Duplicate Paulis were silently collapsing.** `propagate_through_circuit`
  returns a `SparsePauliOp` that may list the same Pauli more than once, so
  `len(op)` is a row count, not a support size, and a `{label: coeff}` dict
  comprehension keeps the last duplicate instead of summing them.
* **monoprop's angle convention is neither the spec's nor Qiskit's.** Its
  `ExpGate` applies `exp(+iθH)`, so the spec's `exp(−iθ_spec/2·G)` needs
  `θ = −θ_spec/2`, *and* the gate list needs reversing like `pauli-prop`'s. All
  eight sign/order combinations were tried against the reference: forward order
  loses terms outright (109 against TFIM's 124), and every wrong angle keeps the
  right support while moving coefficients by up to 1.3. Only `θ = −θ_spec/2` on
  a reversed list lands within 5e-13.

An observable is one scalar and can agree by luck or by cancellation. A
term-for-term diff cannot.

## monoprop is parallel unless you stop it

This is the one caveat that changes a headline number rather than a decimal, so
it gets its own section.

monoprop takes **one serial partition per physical core** when left alone: with a
single MPI rank, `resolve_partition_count_` reads the core count and fans out.
The PyPI wheels are built without MPI (`monoprop.has_mpi == False`), which is
easy to misread as "therefore serial" — it is not. On the 14-core machine below,
Heisenberg at `n=12` measures:

| | wall | CPU (user+sys) | CPU/wall |
|---|---:|---:|---:|
| default | 0.66 s | 6.28 s | 9.48× |
| `monoprop_NUM_THREADS=1 monoprop_PARTITIONS=off` | 1.97 s | 1.97 s | 1.00× |

So an uncapped monoprop reports a **3× faster** wall time than the serial one,
against engines that never had the option. Both variables are read once into a
cached C++ static, so they must be in the environment before the first
propagation — the runner sets them itself as well as receiving them from the
driver.

Because a stale environment variable fails silently, `xbench_monoprop.py`
**measures its own CPU/wall ratio around the timed region and exits non-zero** if
it exceeds `CPU_WALL_MAX` (1.5). Every monoprop row in the run below reported
`cpu/wall=1.00`.

## Known differences between the engines

These are real and are the reason the plot shows the **workload size** next to
the runtime. Read them before quoting a ratio.

* **`pauli-prop` and `monoprop` truncate by a different rule.** Both prune at
  branch-creation time, where the other three accumulate first and drop the
  merged coefficient. `ppvm`, PauliPropagation.jl and PauliStrings.jl then agree
  on the support **exactly**, at every width, on both models; the other two do
  not. On the sweep below `pauli-prop` runs 11–36 % *under* the reference support
  on TFIM (the gap widening with `n`) and 4–15 % *over* it on Heisenberg;
  `monoprop` is within 1 % on TFIM and 4–12 % over on Heisenberg. The driver
  prints a `note:` line for any engine more than 2 % off, and the plot's third
  panel is there so the ratio is never read without it.
* **`monoprop` tracks more rows than it reports.** It retains monomials whose
  coefficient has cancelled to exactly zero. Its `size()` at Heisenberg `n=14` is
  7 355 928 rows against the 3 204 697 terms above threshold — 2.3× — while on
  TFIM the two are within 1 %. The `terms` column is the above-threshold support,
  which is what the other four engines mean by it; `size()` goes to stderr next
  to each row.
* **`pauli-prop` also has a mandatory `max_terms` cap** with pre-allocation. The
  runner defaults it to `2²²` and **fails the run** if the support ever reaches
  it, rather than quietly reporting a differently-truncated number.
* **`monoprop` has a mandatory Pauli-weight `cutoff`**, which the others have no
  analogue of. The runner sets it to `n` — the whole register — so it never binds
  and `lower_atol` is the only truncation in play.
* **PauliStrings.jl stores `im^{#Y}` inside the coefficient** (its `Matrix`
  convention). The dump goes through `op_to_strings`, which puts its
  coefficients on the same real footing as everyone else's.
* **PauliStrings.jl's own front door is `evolve(H, O, tspan; method=Trotter())`,
  which we do not use** — it derives the gate list from the Hamiltonian's
  internal string order, which is not the spec's order. The runner builds the
  `TrotterGate` vector by hand instead.
* **Julia is timed after a warm-up run** so the reported time excludes JIT.
* **Circuit construction is outside the timed region for every engine.**
  monoprop's `propagate` re-expands its gate list internally on each call, which
  no flag hoists out; measured at 0.0–2.2 % of its total at these widths, so it
  is left in rather than worked around.
* Everything is single-threaded (`julia -t1`; `pauli-prop` is single-threaded by
  design; `monoprop` is capped as above; `ppvm` here uses no Rayon).

## A run

Apple M4 Pro (10 P + 4 E cores), macOS 26.6, rustc 1.96.0, Julia 1.12.6,
monoprop 0.8.0, PauliStrings.jl 1.10.1; single-threaded, `--steps 10 --dt 0.1
--j 1 --h 1 --atol 1e-6 --iters 2`. Runtime relative to `ppvm`; **below 1.00×
means that engine beat `ppvm`**. `terms` is the shared support the three
accumulate-then-truncate engines all reach exactly.

`tfim`, widths 8…64:

| n | terms | ppvm | PauliPropagation.jl | PauliStrings.jl | pauli-prop | monoprop |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 4 701 | 0.0028 s | 2.42× | 4.90× | 11.13× | 0.98× |
| 16 | 19 529 | 0.0122 s | 3.10× | 6.26× | 5.53× | 0.69× |
| 24 | 34 353 | 0.0260 s | 5.78× | 8.05× | 5.32× | 0.54× |
| 32 | 49 177 | 0.0431 s | 4.98× | 8.39× | 7.05× | 0.45× |
| 40 | 64 001 | 0.0657 s | 6.33× | 9.20× | 5.28× | 0.43× |
| 48 | 78 825 | 0.0899 s | 7.32× | 9.02× | 5.37× | 0.40× |
| 56 | 93 649 | 0.1188 s | 10.46× | 9.84× | 5.33× | 0.35× |
| 64 | 108 473 | 0.1429 s | 11.26× | 13.18× | 5.71× | 0.33× |

`heisenberg`, widths 6…14:

| n | terms | ppvm | PauliPropagation.jl | PauliStrings.jl | pauli-prop | monoprop |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 1 022 | 0.0028 s | 4.28× | 6.77× | 5.02× | 1.03× |
| 8 | 16 324 | 0.0505 s | 3.43× | 6.26× | 2.02× | 0.83× |
| 10 | 225 353 | 0.5290 s | 3.11× | 6.79× | 1.51× | 0.77× |
| 12 | 1 174 849 | 2.7528 s | 2.59× | 6.50× | 1.13× | 0.67× |
| 14 | 2 915 879 | 6.0906 s | 2.23× | 7.09× | 0.88× | 0.54× |

### Reading it

**`ppvm` beats both Julia engines everywhere**, by 2.2–11.3× against
PauliPropagation.jl and 4.9–13.2× against PauliStrings.jl, and the margin widens
with `n` on TFIM. Those three carry an identical support term-for-term, so those
are clean ratios.

**`ppvm` loses to `monoprop`** on both models — by up to 3.0× on TFIM at `n=64`
and 1.9× on Heisenberg at `n=14`, with the gap widening in `n` on both. On TFIM
the two carry the same support to within 1 %, so that column is a clean loss. On
Heisenberg `monoprop` is doing 4–12 % *more* above-threshold work than `ppvm` and
still finishing sooner, so the gap there is if anything understated — though it
is also tracking 2.3× as many rows in total, so the two engines are not making
the same space/time trade.

**`pauli-prop` is 5.3–11.1× slower on TFIM but competitive on Heisenberg**,
reaching 0.88× at `n=14`. Read its column against the workload note: on TFIM at
`n=64` it is propagating 36 % fewer terms than the reference three, and on
Heisenberg `n=10`–`12` about 15 % more.

## Scaling note

The TFIM support grows roughly linearly in `n` at fixed depth, so that sweep
reaches 64 qubits cheaply. The Heisenberg support grows much faster — three
non-commuting bond rotations per bond per step — and the `Z₀` autocorrelator
itself converges once `n` exceeds the light cone (identical from `n≈10` at
`steps=10`), so widths past that measure scaling rather than new physics.

## Files

- `../../crates/ppvm-pauli-sum/examples/xbench.rs` — the `ppvm` runner.
- `../../julia-benchmarks/benches/xbench_pp.jl` — PauliPropagation.jl.
- `../../julia-benchmarks/benches/xbench_ps.jl` — PauliStrings.jl.
- `xbench_qiskit.py` — `pauli-prop`.
- `xbench_monoprop.py` — `monoprop`, including the thread-cap assertion.
- `run_xbench.py` — validation, driving, CSV merge, summary table.
- `plot_xbench.py` — the figure.

[pp]: https://github.com/MSRudolph/PauliPropagation.jl
[ps]: https://github.com/nicolasloizeau/PauliStrings.jl
[qk]: https://github.com/Qiskit/pauli-prop
[mp]: https://github.com/Algorithmiq/monoprop
