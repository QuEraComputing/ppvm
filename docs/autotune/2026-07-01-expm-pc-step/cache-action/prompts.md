# Approach: cache-action for the real-space expm

Hypothesis: the real-space `mf_expm::expm_apply_mf` recomputes the generator
action (`spec.action` per column) on every Krylov/Taylor matvec (profile:
expm = 97% of pc_step). Build the in-basis generator action ONCE per expm call
as CSC columns (rows + real coeffs), then do cheap cached matvecs — exactly
what the orbit-rep path already does (`build_orbit_rep_cols` + `OrbitRepCscOp`
in orbit_rep.rs), but for the real f64 path (no phases/characters).

Target metric: 15 pc_step_arr steps, N=10 XY ladder, drop=1e-3. Baseline 3.28s.
Track wall (primary) AND rss_mb (memory tradeoff — caching costs memory).
