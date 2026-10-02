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
../../ppvm-python/.venv/bin/python nonclifford.py <bin-dir> <commit>... # ppvm only; writes nonclifford.csv
../../ppvm-python/.venv/bin/python gates.py 50                         # per-gate cost; see gates.log
```

`file` mode takes an optional fifth argument, the qubit count; with `only = ppvm`
it then skips Stim's parse, so programs with `T` gates and rotations run.

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
This is the third such run; per-commit numbers agree with the earlier two
within a few percent.
Stim = official 1.15.0. Raw data: `attribution.csv`.

| Workload | PR 204 | `2fa24ce4` | `d7d23c0e` | `4e8aeb5a` | `b4ea1ce9` | `295aaa8b` | `cfd25143` | `6e8bc36a` | Stim |
|---|---|---|---|---|---|---|---|---|---|
| surface_d30 full | 37.86 ms | 37.61 | **26.29** | **22.38** | 20.47 | **16.72** | 16.62 | **15.49** | 18.26 |
| surface_d30 1 round | 25.45 ms | 25.42 | **13.65** | **9.57** | 9.08 | **6.10** | 6.32 | 6.07 | 6.65 |
| surface_d30 no measure | 9.99 ms | 9.99 | 10.02 | 10.02 | 10.05 | 10.01 | 10.00 | **8.90** | 12.06 |
| surface d7 (118 q) | 160.1 µs | 160.2 | 140.8 | 138.1 | **110.9** | **84.4** | 83.8 | **72.2** | 61.4 |
| surface d11 (274 q) | 734.0 µs | 733.2 | 610.9 | 582.2 | **482.1** | **371.6** | 374.5 | **326.3** | 234.9 |
| surface d19 (778 q) | 6.32 ms | 6.32 | **4.51** | 4.05 | 3.58 | **2.75** | 2.76 | **2.38** | 2.45 |
| repetition d75 (149 q) | 972.0 µs | 975.7 | 983.7 | 1050.0 | **675.5** | **499.0** | 502.5 | **445.1** | 207.1 |
| color d31 (1081 q) | 59.27 ms | 59.51 | 52.35 | 52.18 | 50.46 | 50.02 | **8.11** | 7.13 | 7.08 |
| color d43 (2080 q) | 343.40 ms | 342.69 | 323.06 | 320.34 | 314.09 | 278.39 | **31.42** | 28.95 | 31.70 |

Every commit leaves outcomes and RNG draws unchanged on these circuits; for
`295aaa8b`, `surface_d30` records were checked bit-identical against
`b4ea1ce9` over 22 seeded shots. As agreed, seeded results *can* change for
noisy `MR(p)` and multi-amplitude states (draw order), never the distribution.

### `2fa24ce4` feat(tableau-2): Stim-style `measure_batch`

API only: `GeneralizedTableau::measure_batch` checks which targets are random
in the canonical orientation (contiguous), collapses those under one row guard,
then measures the deterministic ones column-major — Stim's `collapse_z`. Nothing
calls it yet, so no change.

### `d7d23c0e` perf(stim): batch `M`, `MR`, `R` — surface_d30 37.6 → 26.3 ms

The traits-2 executor routes noise-free `M`, `MR` and `R`/`RZ` through
`measure_batch`, applying the `X` resets afterwards (gates need column-major;
`X_q` commutes with `Z_p`, so deferring is exact; repeated targets keep the old
loop). Random measurements now amortize one transpose per instruction and
contiguous projections; all-deterministic instructions never transpose. This is
the largest single win on `surface_d30` (−11.3 ms; −11.8 ms of the 1-round
variant).

**Dead end, not committed:** batching under the *eager* row guard (transpose
every `MR`/`R` regardless) made `surface_d30` 123.6 ms/shot. Under the guard,
every measurement's column reads become strided, including the 29 all-
deterministic rounds. Checking determinism first, as Stim does, is what makes
batching pay.

### `4e8aeb5a` perf(tableau-2): gather columns once — 26.3 → 22.4 ms

A random measurement read the same two X columns at the measured qubit three
times (decomposition, `project_inverse`, `project_row_major`), each a strided
pass under the guard. Now gathered once and passed down. −3.9 ms on
`surface_d30`. `repetition_d75` (deterministic measurements only) is ~6%
slower with this commit, reproduced in three attribution runs; not
investigated (b4ea1ce9 more than recovers it).

### `b4ea1ce9` perf(tableau-2): reuse buffers, skip unused masks — 22.4 → 20.5 ms

Column buffers live in `MeasureScratch` instead of being allocated per
measurement, and a deterministic outcome on a state whose amplitudes are all at
index 0 (every Clifford run) skips widening the masks into the 2048-bit branch
index (one full-width shift and OR per set bit). Small on `surface_d30`
(−1.9 ms) but the biggest win for many-measurement, small-n circuits:
`repetition_d75` −36%, surface d7 −20%.

### `295aaa8b` perf(tableau-2): Stim-style collapse on stabilizer states — 20.5 → 16.7 ms

On a stabilizer state (one amplitude at index 0, inverse signs valid)
`measure_batch` now collapses like Stim's `collapse_qubit_z`: CX appends from
the pivot to the other anticommuting stabilizers, `S` if the accumulated
destabilizer anticommutes with `Z`, `H`, `X` if the outcome needs it. The
destabilizer column is never read and no other destabilizer is multiplied; a
deterministic outcome is one inverse-sign read. The new stabilizer is `±Z`
times other stabilizers rather than `±Z` itself, so the frame differs from the
textbook projection's while describing the same state. This is what puts ppvm
ahead of Stim on `surface_d30` (0.92×) and gives a further −24% on surface d7.

### `cfd25143` perf(stim): batch `MX`, `MY`, `MRX`, `MRY`, noisy `M` — color d43 278.4 → 31.4 ms

These still ran per qubit (`H`, one measurement, `H`) on the column-major
strided path, which profiling showed was 78% of color d31: Stim's color-code
generator ends with `MX`. The executor now rotates every target onto Z at once
around one `measure_batch` (`measure_in_basis`; a repeated target keeps the
per-target order), and `StimTableau::measure_noisy_many` flips each record
after the batch. Color d31 −84%, d43 −89%, which takes the large color codes
level with Stim; circuits without these instructions are unchanged.

### `6e8bc36a` perf(tableau-2): trim gate kernels to live words — every circuit −4% to −14%

A column's stride is rounded up to a whole 4-word block, so below 256 qubits
every gate kernel and the inverse-sign row product swept 4 words where only
`n.div_ceil(64)` can hold a set bit. The padding is zero and every kernel maps
zero words to zero words, so `gate1_mut` / `gate2_mut` and `inv_pair_phase` now
borrow only the live words, and the per-gate `get_disjoint_mut` overlap check
becomes a debug check on ranges that are disjoint by layout. `CX` at n ≤ 128
goes from ~25 to ~15.5 ns (see "Per-gate cost"), surface d7 −14%, surface d19
−14% (now ahead of Stim), `surface_d30` −7%, and the gates-only `no measure`
variant −11%. Tableaus are unchanged: mean outcomes over 200 seeded shots are
identical to `cfd25143` on five Clifford and three non-Clifford circuits.

## Size sweep

`sweep.py` generates Stim memory circuits (`rounds = d`, all four noise
parameters 0.001): rotated surface code, repetition code, and color code
(`C_XYZ`, which ppvm-stim rejects, rewritten as `H` then `SQRT_X_DAG` for both
simulators), plus clifft-bench's `pure_surface_d7_r7`. Branch = `6e8bc36a`.
Raw data: `sweep.csv`, `sweep.log`. Ratios are time / Stim time.

| Surface | d=3 | 5 | 7 | 9 | 11 | 15 | 19 | 23 | 27 | 31 |
|---|---|---|---|---|---|---|---|---|---|---|
| qubits | 26 | 64 | 118 | 188 | 274 | 494 | 778 | 1126 | 1538 | 2014 |
| branch / Stim | 0.99 | 1.06 | 1.17 | 1.33 | 1.39 | 1.09 | 0.97 | 0.96 | 0.88 | 0.86 |
| PR 204 / Stim | 1.87 | 2.32 | 2.61 | 2.93 | 3.14 | 2.60 | 2.57 | 2.47 | 2.23 | 2.07 |

| Repetition | d=3 | 9 | 25 | 75 | 225 | 675 |
|---|---|---|---|---|---|---|
| qubits | 5 | 17 | 49 | 149 | 449 | 1349 |
| branch / Stim | 0.60 | 1.16 | 1.88 | 2.15 | 1.54 | 1.38 |
| PR 204 / Stim | 1.36 | 2.88 | 4.56 | 4.75 | 3.42 | 2.29 |

| Color | d=3 | 5 | 7 | 9 | 13 | 17 | 21 | 25 | 31 | 37 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| qubits | 10 | 28 | 55 | 91 | 190 | 325 | 496 | 703 | 1081 | 1540 | 2080 |
| branch / Stim | 0.68 | 1.30 | 1.77 | 1.68 | 1.57 | 1.44 | 1.20 | 1.11 | 1.00 | 0.97 | 0.92 |
| PR 204 / Stim | 1.23 | 2.82 | 4.14 | 3.67 | 4.52 | 5.09 | 3.06 | 7.15 | 8.37 | 9.25 | 10.71 |

Before `cfd25143` the color row was non-monotonic (d21 at 1.34× between d17 at
3.59× and d25 at 5.75×, `295aaa8b`): that was the unbatched `MX`.

clifft-bench `pure_surface_d7_r7`: branch 72.7 µs, Stim 61.7 µs (1.18×).
Sanity check: mean `1` outcomes per shot agree with Stim's compiled sampler
within ≈2 standard errors on repetition d9, surface d5, color d5 and color d9.
Both color circuits sit ≈2σ *above* Stim (since before this work), which is
worth a dedicated statistical test.

## Non-Clifford regression check

The fast collapse only fires on stabilizer states, but `4e8aeb5a` and
`b4ea1ce9` touch the general measurement path too. `nonclifford.py` times PR 204
against `cfd25143` (3 alternating rounds) on ppvm's `cultivation_d5` and
clifft-bench's non-Clifford circuits, with clifft's half-turn `R_X(a)` / `U3`
rewritten as ppvm's `I[R_X(theta=a*pi)]` / `I[U3(...)]` tags for both builds.
Raw data: `nonclifford.csv`, `nonclifford.log`.

| Program | qubits | PR 204 | `cfd25143` | ratio |
|---|---|---|---|---|
| cultivation_d5 | 42 | 5.98 ms | 5.99 ms | 1.00× |
| msc d3 | 15 | 83.8 µs | 82.0 µs | 0.98× |
| msc d5 | 42 | 5.99 ms | 6.01 ms | 1.00× |
| distillation | 85 | 201.9 µs | 193.7 µs | 0.96× |
| coherent d3 r1 | 26 | 285.3 µs | 280.6 µs | 0.98× |
| coherent d3 r3 | 26 | 4.29 ms | 4.30 ms | 1.00× |
| quantum volume q10 | 10 | 140.1 ms | 140.6 ms | 1.00× |

No regressions. (`coherent_d5` and `quantum_volume_q20` take seconds to minutes
per shot and were left out.)

## Per-gate cost

Profiles of surface d11, repetition d75 and color d9 put `CX` at 39–55% of the
time at small n, with measurement down to 22–25%. `gates.py` isolates it:
`REPEAT 200` of a `CX` brickwork or an `H` layer, no noise or measurement.
Raw data: `gates.log`.

| n | ppvm `CX` before | after (`6e8bc36a`) | Stim `CX` | after / Stim | ppvm `H` before | after | Stim `H` | after / Stim |
|---|---|---|---|---|---|---|---|---|
| 32 | 25.0 ns | 15.6 | 7.4 | 2.11× | 8.0 ns | 4.6 | 3.6 | 1.26× |
| 64 | 24.7 | 15.6 | 7.2 | 2.18× | 7.9 | 4.5 | 3.5 | 1.28× |
| 128 | 23.8 | 15.4 | 9.1 | 1.70× | 7.3 | 4.4 | 5.2 | 0.85× |
| 274 | 27.3 | 19.4 | 14.0 | 1.38× | 10.2 | 8.4 | 7.4 | 1.14× |
| 512 | 27.8 | 20.2 | 25.5 | 0.79× | 11.3 | 10.6 | 9.8 | 1.09× |
| 1024 | 40.7 | 34.5 | 39.6 | 0.87× | 15.6 | 12.9 | 14.6 | 0.89× |
| 2048 | 72.5 | 65.2 | 69.5 | 0.94× | 22.1 | 21.7 | 18.9 | 1.15× |

"Before" is `cfd25143`. Before `6e8bc36a` ppvm's `CX` had a ~25 ns floor up to
~512 qubits: a fixed per-gate cost (the 4-word stride padding, swept twice per
`CX`, plus the per-gate overlap check). What remains at small n is mostly
structural: every gate updates the forward generator phases *and* the inverse
signs (`inv_pair_phase`, two row-phase products per `CX`), where Stim keeps only
its inverse tableau. The forward phases carry the multi-amplitude state
(`odd_phase_destabilizer_mask` in branching, case-a merge and expectations) and
the fallback when the inverse goes stale, so dropping either side is a redesign.

## Open gaps

1. **Per-gate cost at small n.** `CX` is ~15.5 ns against Stim's 7–9 ns up to
   ~128 qubits, mostly the dual phase bookkeeping (see "Per-gate cost"). It is
   most of why small and mid-size circuits (≈30–500 qubits) are still
   1.1–2.2× slower, and why the repetition code, almost all `CX`, trails at every
   size.
2. **Measurement overhead at small n**: `measure_batch_one`, the per-target
   2048-bit "all amplitudes at index 0" check (`memcmp`), the determinism
   pre-scan, `has_repeats`.
3. `4e8aeb5a`'s ~6% slowdown on deterministic-only measurement
   (`repetition_d75`); more than recovered by `b4ea1ce9`.
4. Color codes sit ≈2σ above Stim's mean `1` count (see the sweep).
