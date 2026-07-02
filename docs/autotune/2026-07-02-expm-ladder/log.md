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

## Iteration 1: reuse predictor action-cache in leakage2 — DISCARD (memory)
Implemented (cherry-pick 191e92aa->2411240f): build_mf_cols also retains
out-of-basis action outputs (oob: Vec<Vec<(Word,f64)>>), expm_step returns it,
and a new leakage_from_action_cache consumes it for leak2 instead of recomputing
compute_action_terms on B1. Correct: checksum/sum_a = 1.0, 9/9 crate tests pass
(incl. tight pc_step agreement tests).
Measured N=12 (under Zoom load, so walls unreliable): RSS 1208-1210MB vs baseline
522-527MB -- ~2.3x. RSS is stable run-to-run (load-independent), so this verdict
is robust despite the noisy walls. Wall looked neutral-to-slower but that is load.
Root cause of the memory blow-up: build_mf_cols now ALWAYS builds the oob Word
cache, including for the CORRECTOR expm2 on the doubly-enriched basis B2 (largest
basis of the step), whose cache is then discarded -- pure waste. Even fixing that
(flag to skip oob on the corrector), the predictor's oob cache on B1 still adds
memory for a wall gain bounded by leakage=17% (realistically <=10%).
DECISION: DISCARD. Memory is expm's established weak point (see
[[../2026-07-01-expm-memory]]); ~2.3x RAM for a small, unmeasurable wall gain is
a net loss. If revisited: only build oob for the predictor, and store row-indices
into a preliminary B2 index instead of full Words to bound the memory.
