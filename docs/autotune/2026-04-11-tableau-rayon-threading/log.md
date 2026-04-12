# Log for 2026-04-11-tableau-rayon-threading

## 2026-04-11

### Architecture Notes

**Profiling breakdown (85-qubit MSD circuit, ~135µs total):**
- Encoding (H+T+encode): 55.5µs (41.1%) — includes 5 T gates that branch coefficients
- Middle gates (Clifford): 38.7µs (28.7%) — all Clifford operations
- Measurement: 39.3µs (29.1%) — 85 sequential measurements

**T-gate branching scaling (85 qubits):**
- 1 T gate: 0.2µs (2 coefficients)
- 4 T gates: 2.1µs (16 coefficients)
- 8 T gates: 10.1µs (256 coefficients)
- 12 T gates: 105.1µs (4096 coefficients)
- 16 T gates: 1569.5µs (65536 coefficients)

**Single Clifford gate cost:**
- 32 qubits (64 rows): 43 ns/call
- 64 qubits (128 rows): 82 ns/call
- 128 qubits (256 rows): 161 ns/call

**Key insight:** Parallelizing individual Clifford gate row iteration is NOT viable — per-call overhead (~1µs rayon) far exceeds the gate cost (43-161ns). The primary parallelism target is `branch_with_coefficients` which dominates at high coefficient counts (exponential growth with T gates).

**Parallelism targets (by expected impact):**
1. `branch_with_coefficients` (data.rs) — coefficient loop, O(M) where M grows exponentially with T gates. Each coefficient independently computes phase+branch, but accumulates into shared HashMap. Strategy: parallel fold with local HashMaps, then merge.
2. Measurement overlap computation (measure.rs) — O(M) iteration computing z_overlap.
3. Measurement coefficient update (measure.rs) — processing b_keys entries.

**Design constraint:** All parallelism behind `features = ["rayon"]`. When off, zero regression — identical codepath.
## 2026-04-12
- Iteration 1: Added rayon feature with threshold-based parallel coefficient branching. No regression without feature. With feature enabled, rayon overhead dominates at current benchmark sizes (4096 coefficients). The Vec::collect + par_iter + HashMap fold/reduce pattern adds ~20µs per rayon call. Threshold set to 4096 items; below that sequential path is used. At 65536 coefficients (16 T gates), rayon is still ~2x slower than sequential due to HashMap merge overhead. The per-coefficient work is too cheap for rayon. Next steps: explore chunked accumulation or different data structures for the parallel path.
- Iteration 2: Replaced rayon fold/reduce-with-HashMap with parallel map → Vec + sequential accumulate. The HashMap merge in fold/reduce was the bottleneck (not the computation). New approach: rayon par_map computes all (key, value) pairs into a Vec, then sequential loop inserts into pre-sized HashMap. Results: 35% faster at 32K coeffs, 39% at 131K, 47% at 1M. Threshold raised to 16384. Profile breakdown showed accumulation is 70-85% of total time, computation only 15-25%.
