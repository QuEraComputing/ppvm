# Brainstorm: Direct Ladder-Operator Kernel for Lindblad Apply

## Motivation

`S_± = X ± iY` is ubiquitous (spontaneous emission, superradiance). Currently, a `CollapseOp`
with `c = X_k + iY_k` (2 terms) expands into **4 `LindbladTerm`s** per operator pair. For
n=5 with a dense rate matrix that is 100 terms. The hot loop in `apply` does 2× n-qubit
`comm_parity` scans + ~2× chunk-wise `MulAssign` per (term, state word). For ladder operators
this is wasteful: the action depends only on the **2-bit Pauli at qubit k**, not all n qubits.

## Core insight

For a ladder operator on qubit k with weight γ, the on-site Lindblad action on a Pauli
word W depends only on the 2-bit Pauli at qubit k:

| W_k | Lower (`X+iY`) output    | Raise (`X−iY`) output   |
|-----|--------------------------|--------------------------|
| I   | 0 — skip                 | 0 — skip                 |
| X   | same word × −4γ          | same word × −4γ          |
| Y   | same word × −4γ          | same word × −4γ          |
| Z   | I at k × +8γ, W × −8γ   | I at k × −8γ, W × −8γ   |

The I/Z rows are identical for both directions; only X and Y differ (sign flip). No n-qubit
scan, no Pauli multiplication — just 2 bit reads and at most 2 hash updates.

For cross-site pairs (i≠j where both are Ladders), the action factorises: the left operator
acts on qi independently of the right acting on qj, so the output word is obtained by
`set_new_2(qi, out_qi, qj, out_qj)` with per-qubit lookups. Still O(1) per state word.

Estimated speedup for all-ladder systems (superradiance n=5): **8–15× in the RHS hot loop**.

---

## API design

### New types

```rust
pub enum LadderDirection { Raise, Lower }

pub struct LadderOp {
    pub qubit: usize,
    pub direction: LadderDirection,
}

// User-facing input to LindbladOp::new
pub enum JumpOp<T: Config> {
    Generic(CollapseOp<T>),
    Ladder(LadderOp),
}
```

`LadderOp` carries no generics and no stored Pauli words — just qubit index and direction.
`LindbladOp::new` changes its signature from `Vec<CollapseOp<T>>` to `Vec<JumpOp<T>>`.
The rate matrix still applies uniformly across all operators by position index.

**Also included**: `LadderOp::expand::<T>(n_qubits: usize) -> CollapseOp<T>` converts a
ladder to a two-term `CollapseOp` (X with phase 0 + Y with phase 1 or 3), useful for
feeding into a Hamiltonian or other context. Implementation is ~15 lines: create an
all-identity `T::PauliWordType::new(n_qubits)`, call `set_new(self.qubit, Pauli::X)` and
`set_new(self.qubit, Pauli::Y)`, attach phases, push into a fresh `CollapseOp`. The main
effort is the trait bounds on the generic parameter.

---

## Internal representation

`LindbladTerm` becomes an enum. `LindbladOp` retains a single `Vec<LindbladTerm<T>>`:

```rust
pub(crate) enum LindbladTerm<T: Config> {
    Generic {
        left:      PhasedPauliWord<T::...>,
        right:     PhasedPauliWord<T::...>,
        weight:    f64,
    },
    Ladder {
        qi:        usize,            // qubit of left operator (= ops[i]†)
        qj:        usize,            // qubit of right operator (= ops[j]); qi==qj for on-site
        left_dir:  LadderDirection,  // direction after conjugation: ops[i].direction.flip()
        right_dir: LadderDirection,  // ops[j].direction
        weight:    f64,
    },
}
```

The direction fields are needed because the X/Y entries in the lookup table differ by sign
between Raise and Lower (see core insight table above). `apply` dispatches via `match` on
the variant; `PAR_THRESHOLD` checks `self.terms.len()` as before.

---

## Routing in `LindbladOp::new`

For each (i, j) operator pair with rate `γ_ij`:

| `ops[i]` | `ops[j]` | Result |
|---|---|---|
| `Ladder(li)` | `Ladder(lj)` | One `LindbladTerm::Ladder { qi: li.qubit, qj: lj.qubit, left_dir: li.direction.flip(), right_dir: lj.direction, weight }` |
| `Generic` | `Generic` | M×N `LindbladTerm::Generic` entries (existing path) |
| `Ladder` | `Generic` | Expand Ladder → 2 Pauli terms, then 2×N `LindbladTerm::Generic` |
| `Generic` | `Ladder` | Expand Ladder → 2 Pauli terms, then M×2 `LindbladTerm::Generic` |

`LadderDirection::flip()` is a trivial one-liner (Raise↔Lower swap).

**Why no "half-ladder" optimisation for mixed pairs?** The left-side comm_parity reduces from
an n-bit scan to 1–2 bit reads, but `MulAssign` is still O(n) and dominates. The saving
(~one short loop per term) doesn't justify an extra `LindbladTerm` variant. In practice,
mixed systems are dominated by the Generic-side cost anyway.

---

## Deferred / out of scope

- **Dressed ladders** (`Z_1 ⊗ S_±_k`): users write them as `CollapseOp` manually.
- **LadderOp in Pauli strings**: not needed for the Lindblad use case.
