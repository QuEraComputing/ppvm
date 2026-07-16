# Log for 2026-07-01-expm-pc-step

## 2026-07-01

## Architecture notes (profile before optimizing)

Metric: total wall for 15 `pc_step_arr` steps on the real-space XY two-leg
ladder, N=10 (L=5), gamma=0, dt=0.05, drop_tol=1e-3, localized-Z seed.
Baseline = 3.28 s, RSS 250 MB, peak basis 15414. Conservation Sum a_q = 1.

Phase breakdown (pc_step_timed, steady-state, N=10, drop=1e-3):
- expm  (matrix-exponential apply): **96.8%**  (expm2 70.6% + expm1 26.1%)
- leakage: 3.1%
- generator: 0.0% (matrix-free path builds no CSR)
- expand: 0.1%

=> The entire per-step cost is `mf_expm::expm_apply_mf`. Its `MfOp::dot` /
`spmv_matrix_free` RECOMPUTES the generator action (spec.action per column)
on every Krylov/Taylor matvec. The orbit-rep path already caches the action
once per expm call (`build_orbit_rep_cols` + `OrbitRepCscOp`); the real-space
path does not. Hypothesis 1: cache the action (CSC) once per expm call and
reuse across matvecs -> large expm speedup. Known tradeoff: matrix-full was
~2.8x faster / ~1.43x more memory in an earlier ad-hoc test.

## !! MEASUREMENT BUG (invalidates iterations 1-2 first measurements)

`maturin develop` was run from `crates/ppvm-python-native` (which builds the
crate as its OWN package). The module actually loaded by `import ppvm._core`
is `ppvm-python/src/ppvm/_core.abi3.so`, built via the maturin config in
`ppvm-python/pyproject.toml` — you MUST run maturin from `ppvm-python/`. So
every "measurement" up to and including the first cache-action and loose-table
runs used a STALE .so and is meaningless (baseline 3.28s, cache-action 3.34s,
loose 3.2s were all the same old binary). Correct build command:
`cd ppvm-python && VIRTUAL_ENV=.venv uvx --from maturin maturin develop --release --uv`
Verify freshness: `strings src/ppvm/_core.abi3.so | grep <marker>`.

## Correct profile (freshly-built binary, N=10, drop=1e-3)

Per expm call, the APPLY (Taylor matvecs) dominates: on a 177k-term corrector
apply=355ms (loose, m=15) / 1048ms (double, m=26), vs per-col norm pass 37ms,
setup 10ms. => the SpMV count AND per-SpMV cost both matter; the old
fully-matrix-free op recomputed the Pauli-commutator action on all columns on
EVERY matvec.

## Iteration 1 (as first recorded): cache-action (CSC) — WRONGLY DISCARDED

Ported orbit-rep's `build_orbit_rep_cols` + `OrbitRepCscOp` to the real f64
path in `expm_apply_mf` (build the in-basis action once per call, reuse across
the m*·s matvecs). commit f93a230e (reverted).

Result: wall 3.343 s (baseline 3.275), RSS 264 MB (baseline 250). **No wall
win, slightly more RAM.** Checksum 1.0 (correct).

Diagnosis: caching the action does not help because the action *recompute* is
NOT the bottleneck. Each matvec, matrix-free does compute_action_terms (produce
(row,coeff)) + scatter; cached does read-cols (row,coeff) + scatter — the two
produce the same term stream, so the scatter (memory-bound accumulation over
~basis·O(N) nnz) dominates equally. At dt=0.05 the true cost is the *count* of
scaling-squaring matvecs (m*·s per expm, x2 expm/step), not per-matvec action.
This is algorithmic (Al-Mohy-Higham), not a caching target.

## Root cause + next lever

The expm computes exp to double-precision backward error (`THETA` table is for
u=2^-53 ~ 1e-16) while the basis is truncated at drop_tol=1e-3. That is ~13
orders of over-accuracy. matvec count s = ceil(dt*||A||/theta_m); a looser
theta table (larger theta_m) directly cuts s -> fewer matvecs. Lever 2: use a
single-precision / 1e-6-matched theta table in the mf_expm path. Est ~1.5-2x
expm speedup with no net accuracy loss (still >> truncation). Bigger lever:
Krylov/Lanczos expm (fewer matvecs than Taylor scaling-squaring) — larger
rewrite, deferred.

## Iterations 1+2 combined — KEPT (commit 8eab9190) — ~10.7x

Once the stale-binary bug was found and fixed, BOTH levers are real and
compose:

- cache-action (CSC): build the per-column action once per call
  (`build_mf_cols`) + `CachedCscOp`; the same pass also yields mu/1-norm
  (dropped the separate norm pass). Per-matvec apply 355ms -> 25ms (~14x) —
  the earlier "discard" was purely the stale binary.
- loose table: `select_ms_loose` (tol=1e-6, drop_tol>=1e-4) cuts m ~24 -> 15.

Clean measurements (correct build), N=10 XY ladder, drop=1e-3, min of runs:
- double + no cache (true baseline): 7.15 s, 250 MB
- loose  + no cache:                 5.49 s, 284 MB
- cache + double:                    0.69 s, ~300 MB
- cache + loose (SHIPPED):           0.67 s, ~310 MB

=> ~10.7x wall vs baseline for ~+24% RAM (the CSC cache). Conservation
Sum a_q = 1 and checksum unchanged; all 9 ppvm-lindblad tests pass (drop=0
keeps the double table + cache, still bit-exact vs orbit-rep/merged refs).

Remaining bottleneck is now `build_mf_cols` (the one action pass, ~40ms on the
177k corrector) + the CSC scatter. Next levers if more is wanted: reuse the
cache across predictor+corrector (same basis prefix), or the Krylov rewrite.
