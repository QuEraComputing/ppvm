# Approach: direct-bit-ops

## Hypothesis
The Tableau-level basic Clifford gates (h, cz, cnot, x, y, z, s) currently delegate
to `PhasedPauliWord` methods via the `impl_tableau_clifford!` macro. These methods
add unnecessary overhead per row:
1. A `get_lbit()` check that always returns false (PauliWordNoHash has no loss bits)
2. Redundant bit reads (PhasedPauliWord reads x/z bits for phase, then PauliWord reads them again for the transform)
3. A `rehash()` call that is a no-op (REHASH=false for tableau PauliWords)

The CliffordExtensions (sqrt_x, sqrt_y, etc.) already bypass this by directly
accessing `pw.word.xbits[addr0]` and `pw.word.xbits.set()`. But the basic Clifford
gates still go through the layered delegation.

By rewriting the basic Clifford gates on `Tableau<T>` to match the direct style
already used for CliffordExtensions, we eliminate per-row overhead across 170 rows.

The CZ gate is the hottest target: ~260 calls in the MSD benchmark, each iterating
170 rows. Currently each row does 6 bit reads + 2 writes through bitvec abstraction.

## Target metrics
- gates/single-qubit/h: baseline 107 ns
- gates/two-qubit/cz: baseline 119 ns
- msd/msd-0: baseline 181 µs

## Files to modify
- `crates/ppvm-tableau/src/gates/clifford.rs` — rewrite the `Clifford` impl for `Tableau<T>`
