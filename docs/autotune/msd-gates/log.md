# Autotune Log: msd-gates

## Goal
Improve MSD sampling performance by optimizing the microbenchmark of relevant gates
(sqrt_y, sqrt_x, cz, h, t, measure) used in the MSD circuit.

## Architecture Notes
- MSD uses 85 qubits (5 x 17), `Byte8F64<2>` config (128-bit storage), `u128` index type
- Tableau has 170 rows (2 x 85). Each Clifford gate iterates all 170 rows.
- PauliWord uses `BitArray<[u64; 2]>` for xbits/zbits — REHASH=false for tableau (PauliWordNoHash)
- T gate triggers branching: compute_decomposition O(n^2) + HashMap-based coefficient merge
- Measurement: O(n^2) decomposition + HashMap for overlap + normalize

## Entries

### direct-bit-ops (discard)
Replaced `impl_tableau_clifford!` macro calls for basic Clifford gates (h, x, y, z, cnot, cz)
with direct bit-manipulation implementations matching the CliffordExtensions style.
**Result:** H gate improved ~6% (107→101 ns) but CZ *regressed* ~28% (119→152 ns).
MSD end-to-end 3% slower (181→186 µs).
**Finding:** The compiler was already optimizing the layered PhasedPauliWord delegation chain well.
Manually inlining the same operations didn't help and may have inhibited compiler optimizations
(possibly register allocation or instruction scheduling). The CliffordExtensions direct style
isn't universally faster than the macro delegation. Don't assume bitvec overhead is a bottleneck.

### precompute-phase-mask (keep)
Added `odd_phase_destabilizer_mask()` to precompute a bitmask of destabilizers with odd phase.
New `compute_phase_with_mask()` replaces the O(n) loop with `(active & mask).count_ones()`.
Applied in `branch_with_coefficients`, `compute_coefficients_after_pauli_apply`, and `measure`.
**Result:** MSD **17% faster** (170µs vs 206µs back-to-back). The O(n) loop per coefficient
entry was a significant cost for MSD with 85 qubits and ~32 coefficients.
**Finding:** `compute_phase` was a hidden bottleneck. Precomputing masks that are constant
across the coefficient loop is highly effective. The measurement generalized microbench (32 qubits,
16 coefficients) showed no improvement — the optimization pays off proportionally to n×k.

### fast-deterministic-measure (discard)
Split measurement into two paths: Case B (Z is stabilizer) computes overlap directly from Vec
without HashMap allocation. Case A uses existing HashMap path.
**Result:** ~7µs SLOWER (161µs vs 152µs baseline). Code duplication increased binary size and
hurt instruction cache. With only ~32 coefficients, HashMap<u128, Complex64> construction is cheap.
**Finding:** Don't split hot paths to avoid small allocations. The compiler optimizes the unified
path better. HashMap with 32 entries costs ~200ns to build — negligible.

### avoid-sqrt-normalize (discard)
Replaced `v.abs() * v.abs()` with `v.re() * v.re() + v.im() * v.im()` in SparseVector::normalize.
**Result:** No measurable impact (153µs vs 154µs — noise). LLVM likely already optimizes
`sqrt(x)*sqrt(x)` → `x`. Even if not, 32 sqrts at ~5ns each = ~160ns — negligible.

### bulk-tableau-ops (keep)
Replaced O(n) bit-by-bit stabilizer reset in `update_tableau_according_to_outcome` with
`BitArray::ZERO` bulk assignment + single `set` call.
**Result:** **4% faster** (142µs vs 148µs). The O(85) individual `BitArray::set` calls were
measurably slower than 2 bulk zero operations (zeroing [u64; 2] is 2 stores).

### avoid-clone-measure (REVERTED — correctness bug)
Replaced `self.coefficients.clone().into_iter()` with `std::mem::replace`. Case B trimmed from HashMap.
**Result:** 1.3% faster but introduced a bug: HashMap iteration order is non-deterministic, so
Case B coefficient ordering differed between measurements. Test `test_measure_generalized_idempotent` failed.
**Finding:** HashMap iteration order is not deterministic in Rust. Code that rebuilds Vec from HashMap
changes coefficient ordering. The original clone approach preserves Vec ordering for Case B.
