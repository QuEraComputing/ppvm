# Log for 2026-07-02-expm-ladder

Target: expm (adaptive Pauli-Lindblad pc_step) on the two-leg XY ladder, AFTER
the 2026-07-01-expm-pc-step optimization (cache-action + tol-matched table,
~10.7x). Metric (/tmp/metric_expm.py): N=10, 15 steps, drop=1e-3. Baseline
(min of 3): wall 0.457s, peak_basis 15398, RSS 306MB, checksum 1.0.

## Profile (pc_step_timed phase split)
- expm apply: **77%** (expm2 60% + expm1 17%)
- leakage: **22%** (was ~3% pre-optimization; now significant since expm got 10x)
- generator: 0% (matrix-free), expand 0.7%

## Hypotheses / next iterations
1. Share the `compute_action_terms(p)` computation between the expm CSC build
   (mf_expm::build_mf_cols) and the immediately-following leakage on the SAME
   basis (pc_step_inner: expm1-build on B1, then leakage2 on B1 -- consecutive,
   same basis). The action(p) depends only on p, not coefficients, so it is
   recomputed. CAVEAT: build_mf_cols currently DISCARDS out-of-basis outputs
   (filters index.get(w).is_some()) -- exactly the terms leakage needs -- so
   sharing requires retaining full action outputs (memory cost) and threading
   the cache through pc_step_inner. Structural, medium risk; deferred pending
   the Trotter iteration.
2. Parallelism of leakage vs expm build (check if leakage is already parallel).
