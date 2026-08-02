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
matrices in `PPVM.PauliMatrix`. The n-qubit phase `phaseExpN` here is defined as
the *sum* of those single-qubit exponents (which is exactly what the Rust kernel
computes); it is *also* tied back to genuine `n`-fold tensor-product matrices in
`PPVM.PauliMatrix.tensorPauli_mul`, which proves `phaseExpN` is the base-`i`
exponent of the real `2ⁿ×2ⁿ` operator product `g(p)·g(q)` (via the
phase-multiplicativity `prod_iuPow`, `∏ᵢ iᵏⁱ = i^{Σ kᵢ}`).
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

/-- **The product bits are the identity word exactly on the diagonal.**
`mulWord p q = I ↔ q = p` — the converse of `mulWord_self`. Each qubit's
`(x, z)` pair is `⊕`-cancelling, and `⊕` is injective, so the only word whose
product with `p` has trivial bits is `p` itself. This is what collapses the L4
outer product to its diagonal when one reads off the identity coefficient (see
`PPVM.Twisted.twistedConv_apply_id`). -/
theorem mulWord_eq_id_iff (p q : Word n) :
    mulWord p q = (fun _ => (false, false)) ↔ q = p := by
  have xor_cancel : ∀ a b : Bool, xor a b = false → b = a := by decide
  constructor
  · intro h
    funext i
    have hi := congrFun h i
    simp only [mulWord, mulBits, Prod.mk.injEq] at hi
    exact Prod.ext (xor_cancel _ _ hi.1) (xor_cancel _ _ hi.2)
  · rintro rfl
    exact mulWord_self _

/-- **The product bits are associative.** Per qubit both planes are `⊕`, and
`⊕` is associative, so `(P·Q)·R` and `P·(Q·R)` carry the same bits. Together
with `phaseExpN_cocycle` (the phase half) this is what makes the packed n-qubit
product associative; it is the `kmul`-associativity hypothesis of
`PPVM.Twisted.gtmul_assoc` for the Pauli key. -/
theorem mulWord_assoc (p q r : Word n) :
    mulWord (mulWord p q) r = mulWord p (mulWord q r) := by
  funext i
  simp only [mulWord, mulBits, Prod.mk.injEq]
  exact ⟨Bool.xor_assoc _ _ _, Bool.xor_assoc _ _ _⟩

/-- **Right multiplication by a fixed word is injective.**
`P₁ · Q = P₂ · Q → P₁ = P₂`, because each plane is `p ⊕ q` and `⊕ q` is a
bijection on `Bool`.

This is the machine-checked licence for the `RekeyBijective` fast path in the L4
`Multiply` component: `Sum::mul_word_assign` re-keys every term through
`k ↦ mulWord k q` and merges with a *plain* `insert` plus a `debug_assert!`,
because injectivity means two distinct source keys can never land on the same
destination key. (Injectivity is a real precondition there, not a hint: a
collision would silently drop a term in release rather than sum it.) It is the
word-product analogue of the Clifford re-key's `PPVM.Symplectic.*_bijective` /
`PPVM.Conjugation.conj*_injective`. -/
theorem mulWord_right_injective (q : Word n) :
    Function.Injective (fun p : Word n => mulWord p q) := by
  have xor_cancel : ∀ a b c : Bool, xor a c = xor b c → a = b := by decide
  intro p₁ p₂ h
  funext i
  have hi := congrFun h i
  simp only [mulWord, mulBits, Prod.mk.injEq] at hi
  exact Prod.ext (xor_cancel _ _ _ hi.1) (xor_cancel _ _ _ hi.2)

/-- `P · P` is phase-free: the summed self-phase is `0`. Combined with
`mulWord_self`, `P² = +I`. -/
theorem phaseExpN_self (p : Word n) : phaseExpN p p = 0 := by
  simp only [phaseExpN]
  rw [Finset.sum_eq_zero]
  exact fun i _ => phaseExp_self _ _

end PPVM.PauliWord
