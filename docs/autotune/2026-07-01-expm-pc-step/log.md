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
