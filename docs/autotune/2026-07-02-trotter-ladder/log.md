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
