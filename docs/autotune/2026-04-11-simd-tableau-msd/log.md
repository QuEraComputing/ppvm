# Log for 2026-04-11-simd-tableau-msd

## Goal
Implement ARM NEON SIMD optimizations (behind optional feature flag) to improve generalized
tableau simulation of MSD on Apple M-series chips.

## Architecture Notes
- MSD uses 85 qubits (5×17), `Byte8F64<2>` config (128-bit storage via `BitArray<[u64; 2]>`), `u128` index type
- Tableau has 170 rows (2×85). Each Clifford gate iterates all 170 rows.
- PauliWord uses `[u64; 2]` for xbits/zbits — maps perfectly to NEON `uint64x2_t` (128-bit)
- Prior autotune achieved 32.3% improvement (181→122.5µs) via micro-optimizations
- Key hot paths: PhasedPauliWord MulAssign (fused loop), Clifford gate inner loops, compute_decomposition

## SIMD Targets (ARM NEON)
1. PhasedPauliWord MulAssign: fused loop over [u64; 2] with XOR, AND, OR, NOT, popcount → single NEON 128-bit ops
2. Clifford gate inner loops: 170 rows × bit operations per row → potential for batched NEON
3. compute_decomposition: anticommutes_at checks + Pauli multiplications

## Entries

## 2026-04-11
- neon-mulassign (discard): Replaced scalar loop in PhasedPauliWord MulAssign<&Self> with ARM NEON uint64x2_t intrinsics for [u64;2] storage. No measurable impact (123.7µs vs 123.9µs baseline). LLVM already auto-vectorizes the 2-element loop effectively — the explicit NEON intrinsics provide no benefit over compiler-generated code. For [u64;2], the loop body is only 2 iterations; the NEON version does 1 pass but with the same number of logical operations. Finding: explicit NEON for small fixed-size loops over u64 is redundant on aarch64 — the compiler generates equivalent SIMD already.
- circuit-fusion-match (discard): Fused all ~680 Clifford gates into a single loop with match-based dispatch. 3x SLOWER (370µs vs 126µs). The match dispatch adds ~115K indirect branches across all rows. Branch predictor works perfectly for the original tight per-gate loops (same branch taken 170 times) but poorly for the fused version with 13-variant enum dispatch. Finding: dynamic dispatch via match in the inner loop is catastrophic. Need batch same-type gates (no dispatch) or compile-time specialization.
- batch-same-type (discard): Batch same-type gates (sqrt_y_batch, cz_batch) iterating over indices per row. 1.66x slower than individual calls for 16x sqrt_y (4766ns vs 2877ns). MSD fused: 248µs vs 130µs baseline. Variable addr0 in inner loop prevents compiler from hoisting word_idx/bit_mask computation and constant-folding the bit access. Finding: with AoS layout, per-gate tight loops with constant target are optimal — batching makes each application slower. But combined bitmask approach (merge all same-word gates into single mask op) should work since it reduces to O(1) per word per row.
- combined-bitmask (KEEP): Merged all same-word, same-type single-qubit Clifford gates into combined bitmask operations. For sqrt_y/adj/sqrt_x/adj, compute combined_mask from all qubit indices, then: swap/XOR bits within mask, use popcount parity for phase. Reduces 16x170 iterations to 170 iterations with O(1) work per word per row. Key optimizations: (1) stack-allocated [u64; 8] mask array instead of Vec, (2) fast-path in GeneralizedTableau skipping loss filter when no qubits lost, (3) individual CZ calls (CZ batch was slower due to variable-address penalty), (4) inline const arrays in benchmark (no Vec allocations). Result: MSD 81.4µs vs 127.4µs baseline = 36% improvement. Combined with prior autotune: 181µs → 81.4µs = 55% total improvement.
- cz-block-pairs (KEEP): For CZ pairs with constant qubit offset and same u64 word, replace N individual CZ calls with a single word-level shift+XOR+popcount per row. 17 CZ pairs reduced from 3738ns to 829ns = 4.5x faster. Combined with bitmask batches: MSD 70.4µs vs 126.9µs baseline = 44.5% improvement. Total from 181µs = 61.1% improvement.
- cz-independent-pairs (discard, not committed): Implemented branchless CZ for independent same-word pairs with per-pair conditionals and deferred z_delta accumulation. No measurable improvement (69.7µs vs 70.4µs — noise). The per-pair conditional branches in the inner loop have similar cost to PhasedPauliWord delegation. The compiler already generates near-optimal code for individual CZ (~1ns per row). Only constant-offset CZ can be significantly optimized (via cz_block_pairs).
- PHASE SUMMARY: After 6 iterations (2 kept, 4 discarded), MSD improved from 126.9µs to 70.4µs (44.5%). From original 181µs: 61.1% total. Remaining breakdown: encoding CZ 30.5% (near-optimal), measurement+T 47.6% (sequential/data-dependent), between-block CZ 10.7% (block_pairs), encoding sqrt 11.2% (bitmask batch). Further micro-optimizations on Clifford gates unlikely to yield >1%. Structural changes needed for the next level: measurement path optimization (compute_decomposition, coefficient manipulation), or algorithmic changes to reduce measurement cost.
- cz-block-pairs-cross-word (KEEP): Generalized CZ block pairs to handle cross-word cases where controls are in word 0 and targets in word 1. Replaces 38 individual cross-word CZ calls with 4 batch calls. MSD 63.4µs vs 70.4µs = 10% improvement. From baseline: 63.4µs vs 125.9µs = 49.6%. From original 181µs = 65.0%.
- cz-pairs-branchless (discard, not committed): Direct word-level CZ for independent same-word pairs with per-pair bit extraction and conditional z_delta accumulation. 35% faster in isolation (4.0µs vs 6.2µs for 24 CZ), but <1% end-to-end improvement (62.0µs vs 62.9µs) — absolute savings (~2µs) too small relative to total. cz-lookup-table also discarded — PEXT-like bit extraction loop on ARM negates lookup table benefit.
- ESCALATION: After 8 iterations (3 kept, 5 discarded), micro-optimization well is dry. Remaining 63µs: encoding CZ ~22µs (near-optimal), measurement ~34µs (sequential data-dependent), between-block ~7µs (already optimized). Consecutive failures on encoding CZ (iterations 6, 8) confirm diminishing returns on Clifford gate optimization. Next phase should target measurement path or algorithmic improvements.
