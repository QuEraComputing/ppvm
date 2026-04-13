# Log for 2026-04-12-measure-parallel

## 2026-04-12

### Architecture Notes

**Profiling breakdown for first measurement (max coefficients, 85 qubits):**

| Component | 256 coeffs | 1024 coeffs | 4096 coeffs | 16384 coeffs |
|-----------|-----------|------------|------------|-------------|
| decomp | 1.8% | 0.6% | 0.1% | <0.1% |
| clone (coeffs→HashMap) | 23.9% | 26.8% | 36.1% | 26.3% |
| overlap loop | 21.1% | 22.7% | 17.6% | 22.4% |
| update (drain+norm+filter) | 51.9% | 48.8% | 46.4% | 50.5% |

**Key findings:**
- The clone phase (copying SparseVector into FxHashMap for lookup) wastes ~25-36% of time.
- The update phase (drain B→A, compute norm, filter) dominates at ~50%.
- The overlap loop is ~20%, mostly HashMap lookups.
- decomp is negligible (<1%).

**Optimization strategy (ordered by expected impact):**
1. Eliminate the clone by using a HashMap-backed coefficient storage or building the coeff_map once and reusing it.
2. Parallelize/optimize the update phase (50% of time).
3. Parallelize the overlap loop (20% of time) — pure reduction, embarrassingly parallel.
4. Possibly fuse the overlap and update loops to avoid iterating twice.
- Iteration 1 (eliminate-clone): Replaced clone().into_iter() with mem::replace() in measure(). Avoids allocating a full copy of coefficients. Also inlined case_b trim logic. Results: t8 -5.8%, t10 -6.6%, t12 -4.6%. The improvement is consistent but modest — the clone was ~25% of first-measurement time but only ~10% of total measurement (85 sequential measurements). Next: optimize the update phase (drain+norm+filter, 50% of per-measurement time).
- Iteration 2 (case_b fast path): When stab_anticomm_bits==0, branch_index=idx (self-pairing). Skip HashMap construction entirely — compute overlap and filter directly on Vec. MSD-fused: 63.8→58.2µs (-8.7%), t8: 23.1→17.2µs (-25.5%), t10: 49.5→42.6µs (-13.9%), t12: 156.6→142.4µs (-9.1%). Case_b dominates for lower T-gate counts. Next: optimize case_a path (HashMap overhead for overlap + drain).
- Iteration 3 (avoid retain): Two-pass case_b: iterate by ref for overlap, then filter into coefficients directly. Avoids push-all-then-retain. MSD-fused: -1.3%, t10: -2.7%, t12: -2.0%. Small win — approaching diminishing returns for micro-optimizations.
- Iteration 4 (real-only overlap): Accumulate only real part of z_overlap. Case_b: skip phase 1,3 entries entirely (zero contribution). Case_a: direct Re computation with 2 muls instead of complex chain. t8: -3.3%, t10: -2.5%, t12: noise. Approaching diminishing returns.
- Iteration 5 (retain-based drain): Replaced collect-b_keys + individual-remove with HashMap::retain to partition A/B in a single pass. Huge win: t8 -8.6%, t10 -12.5%, t12 -15.9%. Individual HashMap removes were the bottleneck — retain avoids probe+shift per entry. Total from baseline: t8 -34.8%, t10 -27.3%, t12 -24.3%.
- Iteration 6 (paired overlap, DISCARDED): Tried iterating only A entries and computing both forward+reverse phase contributions per pair to halve HashMap lookups. Regression: t8 +7%, t10 +10%. The extra compute_phase_with_mask call per pair is more expensive than the saved HashMap lookup. Phase computation on u128 (count_ones) is not cheap.
- Iteration 7 (rayon measurement, DISCARDED): Tested rayon for case_a overlap (par_iter sum), case_b overlap, and B→A merge (parallel phase computation). All showed regressions: t14 seq 389µs vs rayon 403µs (+4%), t16 seq 1986µs vs rayon 2379µs (+20%). Root cause: per-element work (phase computation + complex multiply) is too cheap (~5ns) for rayon overhead. The sequential HashMap accumulation dominates and cannot be parallelized. Unlike T-gate branching where rayon helps (heavier per-element work), measurement's bottleneck is HashMap operations.
