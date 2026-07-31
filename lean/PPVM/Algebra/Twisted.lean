/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase

/-!
# The twisted product `key_mul` is associative

The design's L4 keys `PauliSum` on *mod-phase* Pauli words, with `KeyProduct::
key_mul(v, w) = (v ⊕ w, iᵏ)` folding the phase `iᵏ` onto the coefficient
(`traits-2-configuration-and-hashing.md`, "The map is a graded algebra"). That
makes `C[PauliWord]` a **twisted** group algebra: the product of two monomials
`(c, v) · (d, w) = (c·d·i^{phaseExp(v,w)}, v ⊕ w)`.

For this to be an associative algebra the phase cochain must be a 2-cocycle —
which `PPVM.PauliPhase.phaseExp_cocycle` proves. Here we close the loop: over any
commutative ring `C` with a designated fourth root of unity `i` (the design's
`ComplexCoefficient` bound — `i⁴ = 1` is all that is needed), the twisted product
is associative. This is `key_mul` realized directly on the mod-phase key, without
the phase living in a redundant group element.
-/

namespace PPVM.Twisted

open PPVM.PauliPhase

variable {C : Type*} [CommRing C] (i : C)

/-- `iᵏ` for `k ∈ ℤ/4ℤ` — the phase factor `key_mul` folds onto the coefficient. -/
def iPow (k : ZMod 4) : C := i ^ k.val

/-- **`iᵏ` is multiplicative** (`ℤ/4ℤ → C` is a monoid hom into the units), using
only `i⁴ = 1`. This is what turns the additive phase cocycle into a multiplicative
coefficient factor. -/
theorem iPow_add (hi : i ^ 4 = 1) (a b : ZMod 4) :
    iPow i (a + b) = iPow i a * iPow i b := by
  have key : ∀ m : ℕ, i ^ (m % 4) = i ^ m := by
    intro m
    conv_rhs => rw [← Nat.mod_add_div m 4]
    rw [pow_add, pow_mul, hi, one_pow, mul_one]
  simp only [iPow, ZMod.val_add]
  rw [key (a.val + b.val), pow_add]

/-- A mod-phase single-qubit Pauli monomial: a `C` coefficient and `(x,z)` bits.
(`C[PauliWord]` is the free module on these keys; a `Mono` is a single term.) -/
abbrev Mono (C : Type*) := C × Bool × Bool

/-- The twisted product — `KeyProduct::key_mul`: bits `⊕`, coefficients multiply,
and the phase `i^{phaseExp}` is folded onto the coefficient. -/
def tmul (a b : Mono C) : Mono C :=
  (a.1 * b.1 * iPow i (phaseExp a.2.1 a.2.2 b.2.1 b.2.2),
    xor a.2.1 b.2.1, xor a.2.2 b.2.2)

/-- **The twisted product is associative** — so `C[PauliWord]` with `key_mul` is
an associative algebra. Bit associativity is `Bool.xor_assoc`; coefficient
associativity is the phase 2-cocycle `phaseExp_cocycle` transported through
`iPow_add`. -/
theorem tmul_assoc (hi : i ^ 4 = 1) (a b c : Mono C) :
    tmul i (tmul i a b) c = tmul i a (tmul i b c) := by
  simp only [tmul]
  refine Prod.ext ?_ (Prod.ext ?_ ?_)
  · -- coefficient: fold both phase factors, then apply the cocycle
    have hcoc :
        iPow i (phaseExp a.2.1 a.2.2 b.2.1 b.2.2)
            * iPow i (phaseExp (xor a.2.1 b.2.1) (xor a.2.2 b.2.2) c.2.1 c.2.2)
          = iPow i (phaseExp b.2.1 b.2.2 c.2.1 c.2.2)
            * iPow i (phaseExp a.2.1 a.2.2 (xor b.2.1 c.2.1) (xor b.2.2 c.2.2)) := by
      rw [← iPow_add i hi, ← iPow_add i hi, phaseExp_cocycle]
    linear_combination (a.1 * b.1 * c.1) * hcoc
  · exact Bool.xor_assoc a.2.1 b.2.1 c.2.1
  · exact Bool.xor_assoc a.2.2 b.2.2 c.2.2

/-- The identity monomial `1 · I` is a two-sided unit for the twisted product. -/
theorem one_tmul (a : Mono C) : tmul i (1, false, false) a = a := by
  simp only [tmul, phaseExp_id_left, iPow, ZMod.val_zero, pow_zero, mul_one, one_mul,
    Bool.false_xor]

theorem tmul_one (a : Mono C) : tmul i a (1, false, false) = a := by
  simp only [tmul, phaseExp_id_right, iPow, ZMod.val_zero, pow_zero, mul_one,
    Bool.xor_false]

end PPVM.Twisted
