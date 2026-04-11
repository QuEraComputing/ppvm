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
