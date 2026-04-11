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
