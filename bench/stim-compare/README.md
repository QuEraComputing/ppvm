# PR 204 tableau vs Stim

Experiments comparing the `ppvm-tableau-2` generalized tableau (PR 204,
`codex/traits-2-impl`) with Stim's `TableauSimulator` on pure-Clifford circuits,
and the measurement changes on `david/batch-mr-reset` that came out of them.

All numbers: Apple M-series (arm64), single thread, release builds, time per
shot of the circuit execution only (parsing and simulator construction are
outside the timer). On arm64 Stim uses plain 64-bit words (`bitword<64>`; it has
no NEON tableau path), so none of these comparisons involve SIMD.

## Setup

- **ppvm**: `stim-compare` (this crate, its own `[workspace]` so the C++ dependency
  stays out of ppvm's). Runs `execute_validated_with_rng` on a fresh
  `GeneralizedTableau<U2048>` per shot.
- **Stim reference**: `stim_bench.py`, official PyPI `stim` 1.15.0 from
  `ppvm-python/.venv`, one `TableauSimulator.do_circuit` call per shot. Use this,
  not the Rust side's `stim` crate (0.4.5): that crate compiles Stim without
  `NDEBUG`. In practice the difference was small (18.8 vs 18.2 ms on
  `surface_d30`), but the PyPI build is what Stim users run.
- **surface_d30**: `crates/ppvm-stim/examples/surface_d30.stim` (1889 qubits,
  27 870 measurements, `REPEAT 29`). The harness also derives variants by
  deleting instructions: `full-r1` (one repeated round), `no-noise`,
  `no-measure` (drops `M`/`MR`/`R` and the records that reference them),
  `gates-only`.

```bash
STIM_RS_BUILD_FROM_SOURCE=1 cargo build --release
./target/release/stim-compare 20 ppvm [variant]       # surface_d30 variants
./target/release/stim-compare file <circuit.stim> 50  # any circuit
../../ppvm-python/.venv/bin/python stim_bench.py 20 [variant]
../../ppvm-python/.venv/bin/python stim_bench.py file <circuit.stim> 50
../../ppvm-python/.venv/bin/python sweep.py 3          # regenerates circuits/, writes sweep.csv
../../ppvm-python/.venv/bin/python attribute.py <bin-dir> <commit>...  # writes attribution.csv
```

`attribute.py` expects one harness binary per commit in `<bin-dir>`, named by
commit, each built from this crate with the workspace checked out at that commit.

## Starting point: why PR 204 was 2× slower than Stim

On `surface_d30` PR 204 took 37.6 ms/shot against Stim's 18.2 ms. Removing the
measurements made ppvm *faster* than Stim (10.1 vs 12.1 ms), so gates, noise
and SIMD were not the cause. The whole gap was measurement and reset, and
mostly the ~900 random ones (≈27 µs each against Stim's ≈7 µs).

PR 204 already had Stim's complexity: deterministic measurements are O(n/W)
through the inverse tableau, random ones O(n²/W) in the worst case. The cost
was memory access. Profiled with samply:

- `MR`/`R` ran one qubit at a time, column-major, so every random projection
  gathered single generators bit by bit across all 1889 columns
  (`project_inverse` → `gather_row`, ≈11.5 ms) plus a strided per-qubit pass in
  the forward projection (≈6 ms).
- The batched `M` took a row guard up front, so even its deterministic targets
  read their columns strided (`gather_column`, ≈6 ms).

## Attribution

Each commit on `david/batch-mr-reset`, timed in one session, alternating all
binaries and Stim over 3 rounds (median of round medians; PR 204 is `3befaf1b`).
Stim = official 1.15.0. Raw data: `attribution.csv`.

| Workload | PR 204 | `2fa24ce4` | `d7d23c0e` | `4e8aeb5a` | `b4ea1ce9` | `295aaa8b` | `cfd25143` | Stim |
|---|---|---|---|---|---|---|---|---|
| surface_d30 full | 38.46 ms | 38.42 | **26.75** | **22.84** | 20.78 | **16.85** | 16.89 | 18.66 |
| surface_d30 1 round | 26.07 ms | 25.94 | **13.96** | **9.62** | 9.09 | **6.19** | 6.24 | 6.79 |
| surface_d30 no measure | 10.19 ms | 10.22 | 10.25 | 10.26 | 10.24 | 10.26 | 10.22 | 12.32 |
| surface d7 (118 q) | 160.4 µs | 161.3 | 143.6 | 139.7 | **113.0** | **86.2** | 84.7 | 61.7 |
| surface d11 (274 q) | 740.3 µs | 740.8 | 621.4 | 591.0 | **484.5** | **376.5** | 374.6 | 238.4 |
| surface d19 (778 q) | 6.41 ms | 6.40 | **4.61** | 4.09 | 3.63 | **2.77** | 2.79 | 2.49 |
| repetition d75 (149 q) | 988.4 µs | 988.8 | 996.6 | 1060.0 | **682.2** | **505.1** | 504.6 | 207.9 |
| color d31 (1081 q) | 60.15 ms | 60.20 | 53.02 | 52.70 | 51.40 | 50.85 | **8.18** | 7.21 |
| color d43 (2080 q) | 349.93 ms | 348.72 | 330.32 | 327.39 | 321.34 | 283.97 | **32.09** | 32.36 |

Every commit leaves outcomes and RNG draws unchanged on these circuits; for
`295aaa8b`, `surface_d30` records were checked bit-identical against
`b4ea1ce9` over 22 seeded shots. As agreed, seeded results *can* change for
noisy `MR(p)` and multi-amplitude states (draw order), never the distribution.

### `2fa24ce4` feat(tableau-2): Stim-style `measure_batch`

API only: `GeneralizedTableau::measure_batch` checks which targets are random
in the canonical orientation (contiguous), collapses those under one row guard,
then measures the deterministic ones column-major — Stim's `collapse_z`. Nothing
calls it yet, so no change.

### `d7d23c0e` perf(stim): batch `M`, `MR`, `R` — surface_d30 38.4 → 26.8 ms

The traits-2 executor routes noise-free `M`, `MR` and `R`/`RZ` through
`measure_batch`, applying the `X` resets afterwards (gates need column-major;
`X_q` commutes with `Z_p`, so deferring is exact; repeated targets keep the old
loop). Random measurements now amortize one transpose per instruction and
contiguous projections; all-deterministic instructions never transpose. This is
the largest single win on `surface_d30` (−11.7 ms; −12.0 ms of the 1-round
variant).

**Dead end, not committed:** batching under the *eager* row guard (transpose
every `MR`/`R` regardless) made `surface_d30` 123.6 ms/shot. Under the guard,
every measurement's column reads become strided, including the 29 all-
deterministic rounds. Checking determinism first, as Stim does, is what makes
batching pay.

### `4e8aeb5a` perf(tableau-2): gather columns once — 26.8 → 22.8 ms

A random measurement read the same two X columns at the measured qubit three
times (decomposition, `project_inverse`, `project_row_major`), each a strided
pass under the guard. Now gathered once and passed down. −3.9 ms on
`surface_d30`. `repetition_d75` (deterministic measurements only) is ~6%
slower with this commit, reproduced in two separate attribution runs; not
investigated (b4ea1ce9 more than recovers it).

### `b4ea1ce9` perf(tableau-2): reuse buffers, skip unused masks — 22.8 → 20.8 ms

Column buffers live in `MeasureScratch` instead of being allocated per
measurement, and a deterministic outcome on a state whose amplitudes are all at
index 0 (every Clifford run) skips widening the masks into the 2048-bit branch
index (one full-width shift and OR per set bit). Small on `surface_d30`
(−2.1 ms) but the biggest win for many-measurement, small-n circuits:
`repetition_d75` −36%, surface d7 −19%.

### `295aaa8b` perf(tableau-2): Stim-style collapse on stabilizer states — 20.8 → 16.9 ms

On a stabilizer state (one amplitude at index 0, inverse signs valid)
`measure_batch` now collapses like Stim's `collapse_qubit_z`: CX appends from
the pivot to the other anticommuting stabilizers, `S` if the accumulated
destabilizer anticommutes with `Z`, `H`, `X` if the outcome needs it. The
destabilizer column is never read and no other destabilizer is multiplied; a
deterministic outcome is one inverse-sign read. The new stabilizer is `±Z`
times other stabilizers rather than `±Z` itself, so the frame differs from the
textbook projection's while describing the same state. This is what puts ppvm
ahead of Stim on `surface_d30` (0.90×) and gives a further −24% on surface d7.

### `cfd25143` perf(stim): batch `MX`, `MY`, `MRX`, `MRY`, noisy `M` — color d43 284.0 → 32.1 ms

These still ran per qubit (`H`, one measurement, `H`) on the column-major
strided path, which profiling showed was 78% of color d31: Stim's color-code
generator ends with `MX`. The executor now rotates every target onto Z at once
around one `measure_batch` (`measure_in_basis`; a repeated target keeps the
per-target order), and `StimTableau::measure_noisy_many` flips each record
after the batch. Color d31 −84%, d43 −89%, which takes the large color codes
level with Stim; circuits without these instructions are unchanged.

## Size sweep

`sweep.py` generates Stim memory circuits (`rounds = d`, all four noise
parameters 0.001): rotated surface code, repetition code, and color code
(`C_XYZ`, which ppvm-stim rejects, rewritten as `H` then `SQRT_X_DAG` for both
simulators), plus clifft-bench's `pure_surface_d7_r7`. Branch = `cfd25143`.
Raw data: `sweep.csv`, `sweep.log`. Ratios are time / Stim time.

| Surface | d=3 | 5 | 7 | 9 | 11 | 15 | 19 | 23 | 27 | 31 |
|---|---|---|---|---|---|---|---|---|---|---|
| qubits | 26 | 64 | 118 | 188 | 274 | 494 | 778 | 1126 | 1538 | 2014 |
| branch / Stim | 1.07 | 1.25 | 1.39 | 1.46 | 1.59 | 1.24 | 1.11 | 1.07 | 0.96 | 0.91 |
| PR 204 / Stim | 1.81 | 2.26 | 2.63 | 2.92 | 3.11 | 2.61 | 2.57 | 2.47 | 2.22 | 2.07 |

| Repetition | d=3 | 9 | 25 | 75 | 225 | 675 |
|---|---|---|---|---|---|---|
| qubits | 5 | 17 | 49 | 149 | 449 | 1349 |
| branch / Stim | 0.69 | 1.41 | 2.35 | 2.42 | 1.85 | 1.57 |
| PR 204 / Stim | 1.41 | 2.85 | 4.54 | 4.75 | 3.40 | 2.30 |

| Color | d=3 | 5 | 7 | 9 | 13 | 17 | 21 | 25 | 31 | 37 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| qubits | 10 | 28 | 55 | 91 | 190 | 325 | 496 | 703 | 1081 | 1540 | 2080 |
| branch / Stim | 0.76 | 1.56 | 2.17 | 2.00 | 1.73 | 1.66 | 1.36 | 1.19 | 1.14 | 1.06 | 0.98 |
| `295aaa8b` / Stim | 0.78 | 1.94 | 2.85 | 1.99 | 3.18 | 3.59 | 1.34 | 5.75 | 6.99 | 8.00 | 8.69 |
| PR 204 / Stim | 1.26 | 2.83 | 4.17 | 3.70 | 4.53 | 5.18 | 3.04 | 7.15 | 8.36 | 9.29 | 10.74 |

The `295aaa8b` color row is from the previous sweep: its non-monotonic jumps
(d21 at 1.34× between d17 at 3.59× and d25 at 5.75×) were the unbatched `MX`.

clifft-bench `pure_surface_d7_r7`: branch 84.9 µs, PR 204 160.8 µs, Stim 61.5 µs.
Sanity check: mean `1` outcomes per shot agree with Stim's compiled sampler
within ≈2 standard errors on repetition d9, surface d5, color d5 and color d9.
Both color circuits sit ≈2σ *above* Stim (also before `cfd25143`), which is
worth a dedicated statistical test.

## Open gaps

1. **Repetition code: `CX` dominates.** At d675, `CX` is 67% of samples (forward
   `cnot` 32%, inverse-sign `inv_pair_phase` 28%, `gate2_mut` 7%), while on
   `surface_d30` ppvm's gates beat Stim's. Needs a gate-level comparison.
2. **Small and mid-size circuits (≈30–500 qubits) are 1.2–2.4× slower** on every
   family. Not profiled yet; fixed per-instruction and per-measurement costs are
   the likely suspects.
3. `4e8aeb5a`'s ~6% slowdown on deterministic-only measurement
   (`repetition_d75`).
4. Color codes sit ≈2σ above Stim's mean `1` count (see the sweep).
