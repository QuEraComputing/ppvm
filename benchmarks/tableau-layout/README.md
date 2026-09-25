# Tableau storage study

This harness compares the v2 tableau from ppvm PR 204, QuantumClifford.jl,
and a Rust storage experiment. It checks logical results before accepting
timings. The production ppvm storage is unchanged.

Pinned implementations:

- ppvm PR 204: `3befaf1b597033138d6c1bab394b967b09cd57fb`.
- QuantumClifford.jl: `1e553b25fdf12c67b55ccedfbc2f48f4dea78fcc`.

See [REPORT.md](REPORT.md) for findings and recorded results.

## Run

Use Rust, Julia 1.12 and `uv`. Create an isolated QuantumClifford checkout:

```bash
git -C /path/to/QuantumClifford.jl worktree add \
  -b Krastanov/codex/qc-layout-study /path/to/worktree-qc-layout-study \
  1e553b25fdf12c67b55ccedfbc2f48f4dea78fcc
uv run --no-project python benchmarks/tableau-layout/run.py \
  --quantumclifford /path/to/worktree-qc-layout-study
```

The command prepares fixtures, builds both Rust workers in release mode,
creates an isolated Julia environment, validates arbitrary phased Pauli
matrices, runs all workers sequentially, and writes raw and summarized CSV.
New Julia environments reuse the dependency versions in the recorded Manifest,
with local package paths adjusted to the supplied checkout.
The default sweep has 16 sizes from 8 to 2048 qubits, including 31/32/33,
63/64/65, 127/128/129 and 255/256/257. Each configuration has five samples
in each of three independent launches. Each sample has at least 10 ms of
calibrated work, subject to execution-time variation.

For a smoke run:

```bash
uv run --no-project python benchmarks/tableau-layout/run.py \
  --quantumclifford /path/to/worktree-qc-layout-study \
  --sizes 8,65 --samples 1 --launches 1 --min-ms 0.1 \
  --out target/tableau-layout-smoke
```

`--workers native,candidate` selects Rust only. `--julia-project` reuses
an existing environment; the runner checks which QuantumClifford checkout
it loads. `--skip-build` is useful during development; only use it when
the existing binaries were built with the intended compiler flags.
By default, Rust and Julia target the host CPU and Julia uses one thread.
Explicit `RUSTFLAGS` and `JULIA_CPU_TARGET` override CPU settings.
Keep builds, tests, and other heavy work separate from measured runs.

Focused Rust checks and the normal repository test entry point:

```bash
cargo test -p ppvm-tableau-2 --example candidate_storage
cargo test -p ppvm-tableau-2
cargo test --workspace
```

## What each worker measures

| Worker | Layouts and words | Scope |
|---|---|---|
| `tableau_layout.rs` | Native `u64`, row, column, column with batch transpose | Imports the actual private storage module; no duplicated ppvm kernel |
| Same worker, `ppvm-public` | Native canonical column layout | Actual `Tableau` API, including inverse signs and hash invalidation |
| `candidate_storage.rs` | `u8/u16/u32/u64/u128`, qubit bits or generator bits | Experimental padded X/Z vectors with two packed phase planes |
| `quantumclifford.jl` | `UInt8/16/32/64/128`, `fastrow` or `fastcolumn` | Actual QuantumClifford operations on 2n rows |

The candidate is a storage experiment, not a complete simulator backend.
It has no inverse cache, measurement support, quadrant split or extra SIMD
alignment. Both Rust packing axes pad to the selected word width.
Its contiguous-axis operations use whole words; its orthogonal-axis operations
use scalar bit access. This is a controlled storage experiment, not a claim
that these are the fastest possible Rust kernels for every layout.
QuantumClifford `fastcolumn` transposes **words**, while the candidate's
`generator_bits` and ppvm's columns transpose the **bit packing axis**.
QuantumClifford `UInt128`/`fastrow` multiplication is unsupported by SIMD.jl;
the runner requires exactly that exclusion and records it in `unsupported.json`.

## Operations and validation

Each timed iteration sweeps the whole tableau:

- `comm`: symplectic parity of row `r` with `(r+1) % (2n)`, for all rows.
- `mul`: sequential `row[r] = row[r] * row[(r+1) % (2n)]`, including phases.
- `h`, `s`, `cnot`: each qubit, with CNOT target `(q+1) % n`.
- `circuit`: H, S, CNOT at each qubit, in that order, totaling 3n gates.
- `transpose_roundtrip`: ppvm's two physical bit transposes, measured together.

The circuit is a sequence of gates, not dense Clifford matrix multiplication.
Timings are divided by the number of logical row operations or gates.
Transposition is reported per complete round trip.

The Python reference generates valid dense frames from a deterministic
Clifford circuit. A sidecar `.gates` file prepares exactly the same frame
through the ppvm public API. Arbitrary matrices at nine additional sizes
exercise odd phases and both commutation parities. The timed valid frames'
adjacent rows all commute. The arbitrary-input pass checks cases that this
benchmark pattern cannot cover.

Each worker computes a fresh-fixture, one-sweep digest of every logical
X/Z bit and phase (or each commutation result). The runner requires agreement
with the independent Python reference and exact coverage of all expected
configurations and samples. The reference checks all 16 single-site Pauli
products against complex matrices. The Rust candidate also compares its
complete output against scalar Pauli tables at word boundaries.

Fixture construction, initial layout conversion, cloning, compilation,
and checksums are outside the timers. Every sample starts from the same
fixture; repeated sweeps evolve it within the timed batch. These are warm
sweeps, not cold-cache or fixed-pair latency measurements.

## Outputs

`results.csv` contains every sample. `summary.csv` reports the median of
per-launch medians in ns per logical operation, plus sample and launch
ranges. `metadata.json` records source commits, hashes, CPU, toolchains,
flags, and actual loaded Julia package. Julia's Project and Manifest are
retained with results. Worker logs preserve unsupported cases and errors.
`expected.json` and `validation-expected.json` record reference digests.

To generate the size and word-width figures from a summary:

```bash
uv run --no-project --with matplotlib python benchmarks/tableau-layout/plot.py \
  target/tableau-layout
```

Library comparisons also include kernel differences: native row multiplication
copies its source to scratch, while the candidate and QuantumClifford borrow
it; native commutation uses two population-count reductions. Public ppvm
calls include simulator bookkeeping. These measurements do not isolate
padding as the only changed variable.
