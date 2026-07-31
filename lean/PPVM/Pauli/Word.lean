/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase
import Mathlib.Algebra.BigOperators.Fin

/-!
# The n-qubit Pauli word: lifting the phase cocycle

`PPVM.PauliPhase` proved the single-qubit facts. The traits-2 design's key type
`PauliWord` is an n-qubit tensor product, and its packed kernel
(`phase/mul.rs`) accumulates a per-qubit phase and reduces the total mod 4:

```rust
sign_count += sign.count_ones();
imag_count += imag.count_ones();
// ...
self.add_phase(((2 * sign_count + imag_count) % 4) as u8);
```

So the n-qubit phase exponent is the **sum over qubits** of the single-qubit
`phaseExp`, and the n-qubit symplectic form is the sum of the single-qubit
`omega`. This file lifts the group cocycle and the commutation law to that sum,
which is exactly what makes the packed multi-qubit `MulAssign` well defined.

The product bits are pointwise (`x = a^c`, `z = b^d` per qubit).

Scope note: the single-qubit `phaseExp` is grounded against genuine `ℤ[i]`
matrices in `PPVM.PauliMatrix`. The n-qubit phase here is validated as the *sum*
of those single-qubit exponents (which is exactly what the Rust kernel computes),
not re-checked against `n`-fold tensor-product matrices.
-/

namespace PPVM.PauliWord

open PPVM.PauliPhase

variable {n : ℕ}

/-- An n-qubit Pauli word up to global phase: an `(x, z)` bit pair per qubit.
`i`-th slot is qubit `i`, matching the packed `xbits`/`zbits` planes. -/
abbrev Word (n : ℕ) := Fin n → Bool × Bool

/-- Pointwise product bits: `⊕` on each qubit's X and Z plane. -/
def mulWord (p q : Word n) : Word n :=
  fun i => mulBits (p i).1 (p i).2 (q i).1 (q i).2

/-- The n-qubit phase exponent: the sum over qubits of the per-qubit
`phaseExp`, reduced in `ℤ/4ℤ` — the spec of `(2*sign_count + imag_count) % 4`. -/
def phaseExpN (p q : Word n) : ZMod 4 :=
  ∑ i, phaseExp (p i).1 (p i).2 (q i).1 (q i).2

/-- The n-qubit symplectic form: sum over qubits of the single-qubit `omega`. -/
def omegaN (p q : Word n) : ZMod 4 :=
  ∑ i, omega (p i).1 (p i).2 (q i).1 (q i).2

/-- **The n-qubit 2-cocycle law.** The summed phase is still a 2-cocycle, so the
packed multi-qubit product is associative. Lifts `phaseExp_cocycle` termwise. -/
theorem phaseExpN_cocycle (p q r : Word n) :
    phaseExpN p q + phaseExpN (mulWord p q) r
      = phaseExpN q r + phaseExpN p (mulWord q r) := by
  simp only [phaseExpN, mulWord, mulBits]
  rw [← Finset.sum_add_distrib, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => phaseExp_cocycle _ _ _ _ _ _

/-- **The n-qubit commutation law.** `phaseExpN p q − phaseExpN q p = 2·ω(p,q)`,
i.e. `P·Q = (−1)^{ω(P,Q)} Q·P` for n-qubit words. Lifts `phaseExp_sub_comm`. -/
theorem phaseExpN_sub_comm (p q : Word n) :
    phaseExpN p q - phaseExpN q p = 2 * omegaN p q := by
  simp only [phaseExpN, omegaN, Finset.mul_sum]
  rw [← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun i _ => phaseExp_sub_comm _ _ _ _

/-- Product bits are involutive per qubit, so `P · P` has identity bits. -/
theorem mulWord_self (p : Word n) : mulWord p p = fun _ => (false, false) := by
  funext i
  simp only [mulWord, mulBits, Bool.xor_self]

/-- `P · P` is phase-free: the summed self-phase is `0`. Combined with
`mulWord_self`, `P² = +I`. -/
theorem phaseExpN_self (p : Word n) : phaseExpN p p = 0 := by
  simp only [phaseExpN]
  rw [Finset.sum_eq_zero]
  exact fun i _ => phaseExp_self _ _

end PPVM.PauliWord
