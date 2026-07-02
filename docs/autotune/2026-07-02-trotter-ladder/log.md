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
