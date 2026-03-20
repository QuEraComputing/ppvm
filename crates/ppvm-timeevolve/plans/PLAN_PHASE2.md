# Performance Brainstorm: ppvm Lindblad Solver (Phase 2)

## Context

The ppvm-timeevolve solver implements DOPRI5 for the Lindblad master equation in the Pauli
operator basis. Tasks 11–15 are complete (allocations eliminated, FSAL, single-multiply
anticommutator, per-step truncation). The code already uses `ByteF64<N>` =
`HashMap<W, f64, FxBuildHasher>` — no new crates needed. Target: n=20+ qubits, single
solve, rayon acceptable.

**Cost per RHS at n=20 (dense rate matrix, |P| ≈ 500):**
- ~1600 Lindblad terms × 500 entries × 2.5 MulAssigns = ~2M MulAssigns per call
- 6 rhs_into calls/step (FSAL: k[0] reused from previous step) → bottleneck is the O(|Γ|×|P|) inner loop

---

## ODE Solver Choice

For the truncated Lindblad system, **DOPRI5 with adaptive step is the right choice**.
Lower-order methods require many more steps for the same tolerance:

| Method | RHS/step | Steps at ε=1e-6 | Total cost |
|--------|----------|-----------------|------------|
| Euler (fixed) | 1 | ~10^6 | catastrophic |
| RK4 (fixed) | 4 | ~1000 | 4000 |
| RK23 (adaptive 2(3)) | 3 | ~500 | 1500 |
| **DOPRI5 (adaptive 4(5))** | **6** | **~100** | **600** |

Truncation removes fast-decaying modes, so the effective system is not stiff after the
initial transient — explicit adaptive methods work well. **The right lever is to match
`rtol` to the truncation tolerance**: if Budget truncates at 500 entries (accuracy ~1e-3),
set `rtol = 1e-3` and DOPRI5 takes large steps automatically.

---

## Task 16 — Commutator form in `LindbladOp::apply`

**Mathematical basis.**
```
2·left·P·right − {a_kl, P}  =  [left, P]·right  +  left·[P, right]
```
The total contribution for state entry W_a is non-zero iff `left` or `right` anticommutes
with W_a. Define `p1 = comm_parity(left.word, W_a.word)`, `p2 = comm_parity(W_a.word,
right.word)`, `multiplicity = p1 + p2` ∈ {0, 1, 2}:

- multiplicity=0 (prob≈25% at large n): contribution is **zero** — skip entirely
- multiplicity=1 (prob≈50%): coefficient = `2·weight·re_phase(tmp.phase)·coeff_a`
- multiplicity=2 (prob≈25%): coefficient = `4·weight·re_phase(tmp.phase)·coeff_a`

Expected MulAssigns/entry: 0×0.25 + 2×0.5 + 2×0.25 = **1.5** (vs. current **2.5**).

---

## Task 17 — Packed `comm_parity` using byte-level bit operations

**Mathematical basis.** `comm_parity(a, b) = popcount((a.xbits & b.zbits) XOR (a.zbits & b.xbits)) mod 2`.
For `BitArray<[u8; N]>`, the raw byte array is accessible via the bitvec API. For NBYTES=3:
```rust
let mut p = 0u8;
for i in 0..N {
    p ^= ((a.xbits.data[i] & b.zbits.data[i]) ^ (a.zbits.data[i] & b.xbits.data[i]))
           .count_ones() as u8;
}
parity = p & 1
```
This is 3 iterations instead of n=20 per-qubit iterations.

---

## Task 18 — Rayon parallelism over Lindblad terms

Parallelize the outer `self.terms` loop in `LindbladOp::apply` using rayon with
thread-local accumulators via `fold`/`reduce`. No manual merge loop needed.

---

## Task 19 — `Budget` truncation strategy

Add a `Budget` truncation strategy that caps `|P|` at a target entry count. Coupling:
set `rtol = min_threshold` in `SolverConfig` to match ODE error to truncation error.
