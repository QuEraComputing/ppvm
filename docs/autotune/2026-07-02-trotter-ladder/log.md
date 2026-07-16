# Log for 2026-07-02-trotter-ladder

Target: Trotter (PauliSum rxx/ryy/rzz gate propagation) on the two-leg XY ladder.
Metric (/tmp/metric_trotter.py): brickwork 2nd-order Strang, N=10 (L=5), dt=0.1,
8 steps, min_abs_coeff=1e-3, localized-Z seed. Baseline (min of 3, quiet build):
wall 0.777s, peak_basis 34030, RSS 173MB, M=1.0 (conserved).

BUILD REMINDER: rebuild native from ppvm-python/ (not crates/ppvm-python-native)
  cd ppvm-python && VIRTUAL_ENV=.venv uvx --from maturin maturin develop --release --uv

## Profile (Python-level phase split: gates with truncate=False vs .truncate())
- gate application (rxx/ryy_many via map_insert): **98%**
- truncation (drop below min_abs_coeff): 2%

=> The bottleneck is the 2-qubit rotation application, not truncation. Each gate
calls PauliSum::map_insert -> ACMapInsert::map_insert, whose HashMap impl
(ppvm-traits/src/map/hashmap.rs:157) is a SINGLE-THREADED `for (k,v) in
self.iter_mut()`. The closure is already `Fn + Sync + Send`. The expm path
parallelises the analogous action pass (mf_expm::build_mf_cols) with rayon.

## Iteration 1 (running): parallelise map_insert gate application
Hypothesis: rayon-parallelise the per-entry closure (independent `*v *= cos`
mutation + collected branch terms), merge branch terms into dest. Expected
~cores x on the gate phase (98% of step). Small-map sequential fallback.

## Iteration 1: parallelize map_insert — DISCARD (~30% slower)
Cherry-picked 3d5f0409 (subagent parallelized map_insert_vec/map_insert via
collect-mut-refs + par_iter + par_extend, threshold 1024). Result: wall 1.0s
(min) vs 0.777s baseline, RSS 230 vs 173MB. Correct (M=1.0), but SLOWER.
Root cause: a Trotter step applies ~120 gates (rxx+ryy over ~60 bonds x2
directions), and EACH gate calls map_insert over the full ~34k-term map. The
per-call `iter_mut().collect::<Vec<(&W,&mut C)>>()` (~0.5MB alloc) + rayon
dispatch is paid 120x/step; the per-term work (a few-ns commutation check) is
memory-bandwidth-bound, so parallel gain is small and the per-call overhead
dominates. Opposite of the expm path (ONE big par pass per call over the basis).
Anti-pattern: don't rayon-parallelize a cheap loop that's invoked many times per
step. Better directions (next): (a) skip non-overlapping terms via a
qubit-support index so a gate on bond (a,b) visits only terms touching a/b
instead of all 34k (algorithmic, cuts the O(gates x basis) work); (b) if
parallelizing, do it at a coarser grain with no per-gate collect (hashbrown
native par_iter_mut), only for very large maps.

## Iteration 1 REVISITED (per user: "regression an artifact of too-small a problem?") — KEEP
Re-tested parallel vs sequential at larger N, clean/quiet, min of 3:
- N=10 (L=5, 8 steps, peak 34k): seq 0.777s vs par 0.702-0.77s -> NEUTRAL.
  (The earlier "1.0s / +30%" was measured under machine load — seq was also
  inflated to 1.17s then; on a quiet machine parallel N=10 matches sequential.)
- N=12 (L=6, 6 steps): seq 18.4s (min) vs par 7.5s (min) -> ~2.5x FASTER.
The user was right: the regression was a small-problem + load artifact. At
N=12 each gate operates on a large intra-step map (gates deferred with
truncate=False grow the map before the single per-step truncate), so the
per-term parallel work amortizes the rayon dispatch overhead. Correct (M=1.0).
Cost: RSS +60% at N=12 (627 vs 396MB — per-thread accumulation buffers).
DECISION: KEEP (helps where Trotter is slow; neutral where it's fast).
Lesson: benchmark the autotune metric at a size representative of the SLOW
regime, not the fast one — a small metric hid a real 2.5x structural win.

## Iteration 1 FINAL VERDICT (clean quiet machine) — DISCARD
The "KEEP 2.5x" above was ALSO a measurement artifact. Re-measured parallel vs
sequential BOTH freshly-built on a quiet machine (min of 3-5):
  N=12  sequential 8.06s / 396MB   vs   parallel 7.25s / ~600MB
=> only ~1.1x faster for +50% RAM. The earlier 18.4s "sequential" was
load-inflated (~2.3x); the earlier 30%-slower "parallel" was also load. Both
prior verdicts were corrupted by single-run walls under variable machine load.
Clean physics: the gate loop calls map_insert ~120x/step (short parallel
sections) on memory-bandwidth-bound work, so parallel speedup is small and the
per-gate collect overhead nearly cancels it; the per-thread buffers cost +50%
RAM. DISCARD (reverted to sequential).
PROCESS LESSON (critical): NEVER decide keep/discard from single-run walls under
uncontrolled load. Require: quiet machine, both variants freshly built, min-of-N,
ideally back-to-back. Two flip-flops here came from ignoring this.

## CONCLUSION
No net change kept. Trotter gate application (98%) is memory-bandwidth-bound and
invoked ~120x/step; naive rayon parallelization gives only ~1.1x for +50% RAM
(discarded). The only real remaining lever is architectural — columnar Pauli
bit-plane storage to make the per-gate pass cache-friendly + SIMD-able before
parallelizing — deferred (major rewrite, uncertain payoff). Trotter is near its
practical optimum for the current HashMap-of-Paulis representation.
