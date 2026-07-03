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

## CONCLUSION
No net change kept. Post the prior ~10.7x (2026-07-01-expm-pc-step), the profile
is 82% scaling-squaring matvecs (already cached + parallel + tolerance-matched)
and 17% leakage. The one structural idea (share the action pass into leakage2)
is correct but costs ~2.3x RAM (discarded). expm is near its practical optimum;
the earlier 10.7x per-step speedup stands as the real win.

## tau_add / K admission study (both models, K-scan, high + low precision)

Motivation: the doubly-enriched PC transient is ~11-16x the pruned basis, ~94%
pruned (load-independent). Q: does a rate-based admission filter (tau_add =
K*drop_tol/dt, PPVM_K_LEAKAGE) beat admit-all / a max_basis cap? Judged at
MATCHED accuracy vs exact ED. (One K on one model can't decide it.)

Ladder N=10 (dt=0.05, T=2), min-of-2, K-scan:
- loose (rel~1e-2): K=1 ~2x faster than admit-all (2.3s vs 4.7s), better rel.
- high  (rel~1e-3): K=1..3 ~1.2-1.4x faster (15s vs ~18-21s); never worse.
- K=5 OVER-filters (rel 6e-2 at drop=3e-4 vs admit-all 6.8e-2 only by using a
  tiny basis) -> bad at matched accuracy. K=3 borderline; K=1 best.
- max_basis cap: ~3.6x at LOOSE accuracy (rel~1.3e-2, 0.9s) BUT cannot reach
  high precision (the cap clips the basis below what rel<=1e-3 needs). It is a
  memory-bound tool, not an accuracy knob.

Long-range XY chain N=10 (alpha=1, gamma=0), min-of-2:
- same pattern: K=1 ~1.1-1.34x faster than admit-all at matched accuracy
  (rel~1.2e-2: 26s vs 35s; rel~2.8e-4: 46s vs 52s). K=3/K=5 over-filter.

Conclusions:
1. Removing the admission filter (running admit-all) was a modest-to-moderate
   LOSS: a small K (~1) recovers ~1.1-2x at matched accuracy on BOTH models,
   MORE at loose accuracy, and never hurts. K auto-scales (tau_add∝drop_tol), so
   at tight drop it filters little and degrades gracefully to admit-all.
2. K is regime-dependent. K=1 is the sweet spot for these UNITARY (gamma=0)
   cases; K>=3 over-filters. The notes' K≈5 was calibrated on a DISSIPATIVE
   (gamma=1) long-range case -- dissipation damps high-weight operators and
   changes the leakage structure, so its optimal K is higher. Judging tau_add
   from K=5 on a unitary ladder (as first done) was unfair.
3. max_basis is a hard memory bound (useful, and good at loose accuracy) but NOT
   a substitute for K at high precision (it clips accuracy). K works at all
   precisions.
Recommended: expose K (done, PPVM_K_LEAKAGE, default 0=off). A default of K≈1
would be a Pareto improvement for unitary runs; dissipative runs may want higher.
Default left at 0 pending a decision (K is regime-dependent).

## RAM (peak RSS) is the real metric — admission control is NOT a RAM lever

Re-judged the K/admit-all comparison on PEAK RSS (max throughout runtime) at
MATCHED accuracy (RSS is load-independent, so robust).

N=10 high precision (rel~1e-3): peak RSS ~FLAT across K (~390MB ladder, ~510MB
lrxy) — N=10 saturates (~260k of 4^10 reachable), so no transient to filter.
N=12 (non-saturated), drop=1e-4:
    K=0 : rel 7.0e-3, RSS 1772MB   (peak_basis 1.60M)
    K=1 : rel 5.2e-3, RSS 1715MB   (peak_basis 1.65M)   <- matched accuracy, ~same RAM (3%)
    K=3 : rel 2.0e-2, RSS  997MB   (peak_basis 0.87M)   <- 1.8x less RAM but 3x WORSE rel

Conclusion:
- At MATCHED accuracy, the admission filter does NOT reduce peak RAM (K=0 vs K=1:
  1772 vs 1715MB). K is a WALL optimization (~1.3-2x) only. Higher K (K=3) cuts
  RAM ONLY by cutting accuracy — moving DOWN the RAM<->accuracy Pareto, not
  beating it.
- Peak RAM is fundamentally set by the retained basis the accuracy target
  requires (+ its action cache + the doubly-enriched corrector transient B2,
  which at high precision is close to the retained basis and barely filterable).
- The direct knob for a RAM BUDGET is max_basis (hard bound), at a known
  accuracy cost (capping below what the accuracy needs degrades rel steeply).
- => For RAM specifically, removing tau_add and keeping max_basis was FINE:
  max_basis IS the RAM knob; tau_add/K would not lower peak RAM at matched
  accuracy. (This tempers the earlier "removing tau_add was a loss" — that was a
  WALL statement; for RAM it is not a loss.)
- The only way to cut peak RAM at FIXED accuracy is a code change to the peak
  itself: reduce per-term footprint or stream the corrector's action instead of
  caching all of B2 (trades back toward the recompute cost the cache removed).

## CORRECTION: tau_add and max_basis limit RAM by DIFFERENT mechanisms

(User's point.) The earlier "admission control is not a RAM lever" was measured
with tau_add(K) only, both at max_basis=inf -- which missed a real mechanistic
difference:
- max_basis room-cap is applied PER CHUNK inside leakage_with_prune
  (merged.retain(top-room)) AND caps B2 via add_leakage_capped -> it bounds the
  LIVE transient (leakage map + doubly-enriched B2) during accumulation.
- tau_add (as implemented) filters only at the END -> the merged map still
  accumulates ALL candidates first, so it shrinks the RESULT but NOT the peak
  transient memory.

N=12 ladder, drop=1e-4, MATCHED accuracy (rel~7e-3), peak RSS:
    admit-all (mb=inf)      : 1795 MB
    tau_add K=1 (end-filter): 1715 MB  (-4%)
    max_basis = 1.8M cap    : 1611 MB  (-10%, rel 6.7e-3, cap just above retained
                                        ~1.6M so accuracy preserved; peak_basis
                                        1.52M shows the cap clipping the transient)
    (mb = 2.0M/2.5M ~ 1759/1838 MB: cap above the overshoot -> no effect)

=> max_basis is the better RAM knob: it clips the transient overshoot (B2 >
retained) during accumulation. The edge is ~10% at N=12 high precision (small
overshoot there) and LARGER at looser accuracy (N=10 rel~1e-2: 233 vs 306MB,
~1.3x) where B2 overshoots the retained basis much more (up to ~16x). tau_add
COULD match this if it filtered PER CHUNK (bounding the live map) rather than at
the end -- a code refinement. For peak RAM, set max_basis ~1.1-1.3x the expected
retained basis: preserves accuracy, clips the transient, hard-bounds RAM.
