# Old-vs-`-2` regression report

`mise run perf-report` screens every `ppvm-conformance-2` benchmark target and
writes `target/perf-report/report.md`. Each of those targets links **both**
engines into one binary, so a `<id>/old` + `<id>/new` pair is a same-build
ratio and the executable-layout bias cancels inside it. Keep this around until
the `-2` crates replace their legacy twins outright — it is the standing check
that the refactor has not lost ground.

Ratios are `new / old`: below 0.97 is an improvement, 0.97–1.03 is parity, and
above 1.03 is a regression. A regression only *blocks* when the slowest of
several processes also stays above the gate, which is what separates a real
one from an executable-layout artifact — so a one-launch screening run reports
but never fails.

```bash
# Screen everything (long; run it on an otherwise idle machine).
mise run perf-report

# Narrow to one target, or one Criterion filter.
mise run perf-report -- --bench tableau_surface_bench
mise run perf-report -- --filter noise

# Confirm suspects in fresh processes. This is the invocation that can fail.
mise run perf-report -- --filter clifford --launches 4

# Re-summarize an existing raw log without re-benchmarking.
mise run perf-report -- --reuse target/perf-report/raw.txt
```

Outputs land in `target/perf-report/`: `raw.txt` (concatenated Criterion
output), `pairs.tsv` (every pair, slowest first), and `report.md`.

Accepted regressions go in `perf-allowlist.txt`, one `fnmatch` pattern per
line with a comment explaining the acceptance. It is empty on purpose; the
standing policy is in the "Perf-drift allowlist" section of `docs/log.md`.

- `perf_regression_report.py` — the runner, classifier, and report writer.
- `summarize_criterion_output.py` — the pairing/median parser it reuses; also
  usable standalone against a hand-collected log.

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
  cargo run --release -p ppvm-pauli-sum --example trotter_qubit_sweep \
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

- `../crates/ppvm-pauli-sum/examples/trotter_qubit_sweep.rs` — ppvm sweep,
  replicating the Python storage-tier dispatch for `[u8; N]` storage.
- `../julia-benchmarks/benches/trotter_sweep.jl` — PauliPropagation.jl sweep.
- `plot_tfim_sweep.py` — renders the log-y comparison from the two CSVs.

[pp]: https://github.com/MSRudolph/PauliPropagation.jl

# Branch-coalesce scaling: sort-merge vs FxHashMap

Follow-up study for PR #154, which replaced the `FxHashMap` coalesce in the
T-gate hot path (`GeneralizedTableau::branch_with_coefficients`) with a
sort-merge and measured ~10× on `cultivation_d5`. This harness asks whether
that win **persists as the branch count `m` grows**, and **where the hash
coalesce wins again**. Because #154 deleted the hash path from the default
build, the bench reimplements *both* coalesce routines (faithful ports, asserted
equivalent at start-up) and drives them with identical real inputs at
`m = 2^k` (`k` branching T gates on an 80-qubit, `u128`-indexed tableau).

Two collision regimes:

- **doubling** — the next T flips a fresh index bit (output `2m`, zero merges);
  the canonical per-T-gate cost. Sort-merge wins throughout and the gap
  *widens* with scale (≈3.8× at `m = 2^20`).
- **merge** — the next T flips a bit the set is already closed under (output
  `m`, all collisions); the flavour of the measurement case-a path. The hash
  coalesce overtakes sort-merge for `m ≳ 2048` (the dense-collision regime is
  where probing's free coalesce-on-insert beats paying for a full sort).

## Reproduce

```bash
# 1. Run the bench (writes target/criterion/branch-coalesce-*/...).
#    Default sweep tops out at m = 2^20; bump with PPVM_BRANCH_MAX_EXP.
cargo bench -p ppvm-tableau --bench branch-coalesce-scaling
# PPVM_BRANCH_MAX_EXP=22 cargo bench -p ppvm-tableau --bench branch-coalesce-scaling

# 2. Plot (reads criterion's estimates.json directly — no CSV step).
uv run --with matplotlib python benchmarks/plot_branch_coalesce.py \
  --out /tmp/branch_coalesce_scaling.png
```

- `../crates/ppvm-tableau/benches/branch-coalesce-scaling.rs` — the A/B bench.
- `plot_branch_coalesce.py` — left panel: time vs `m` (log-log); right panel:
  sort-merge speedup `t_hash / t_sortmerge` vs `m`, with the crossover line and
  the "hash wins" band.
