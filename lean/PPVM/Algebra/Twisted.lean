/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase
import PPVM.Pauli.Word
import PPVM.Algebra.GradedMap

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
`ImaginaryUnit` bound, weakened from the earlier `ComplexCoefficient` — `i⁴ = 1`
is all that associativity needs), the twisted product is associative. This is
`key_mul` realized directly on the mod-phase key, without the phase living in a
redundant group element.
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

/-! ### The obligation on an arbitrary `KeyProduct`

`tmul_assoc` above is about the *concrete* Pauli phase cochain. But the Rust
trait `KeyProduct` is deliberately key-agnostic (`key_mul(&self, &Self) ->
(Self, Phase)`) and its doc asserts, for every implementer, that "the phase
exponent is a genuine 2-cocycle (hence the twisted product is associative)". The
design anticipates further keys (an ordered fermionic word over `FermionSite`),
for which L4 `Multiply` would inherit that claim with nothing backing it.

So we state the claim where it belongs: for an arbitrary key set `K` with a
product `kmul` and a `ℤ/4`-valued 2-cochain `β : K → K → ZMod 4`, the twisted
product is associative **exactly when** `kmul` is associative and `β` is a
2-cocycle. Every `KeyProduct` impl owes these two hypotheses; the Pauli case
discharges them via `Bool.xor_assoc` and `phaseExp_cocycle` and is recovered
below as an instance (`tmul_assoc_of_gtmul`) rather than being the whole content. -/

variable {K : Type*}

/-- The twisted product on `C × K` for an abstract key product `kmul` and phase
2-cochain `β` — the general shape of `KeyProduct::key_mul` folded onto the
coefficient. -/
def gtmul (kmul : K → K → K) (β : K → K → ZMod 4) (a b : C × K) : C × K :=
  (a.1 * b.1 * iPow i (β a.2 b.2), kmul a.2 b.2)

/-- The 2-cocycle condition on an abstract phase cochain: the obligation a
`KeyProduct` implementer must discharge. -/
def IsCocycle (kmul : K → K → K) (β : K → K → ZMod 4) : Prop :=
  ∀ u v w, β u v + β (kmul u v) w = β v w + β u (kmul v w)

/-- **`KeyProduct`'s associativity obligation, in general.** Over any commutative
ring `C` with `i⁴ = 1`, any key product `kmul`, and any `ℤ/4`-valued phase
cochain `β`: if `kmul` is associative and `β` is a 2-cocycle, then the twisted
product `(a,u) ⊙ (b,v) = (a·b·i^{β u v}, kmul u v)` is associative. This is the
precise law every `KeyProduct` impl owes — the Pauli word is one instance. -/
theorem gtmul_assoc (hi : i ^ 4 = 1) (kmul : K → K → K)
    (hk : ∀ u v w, kmul (kmul u v) w = kmul u (kmul v w))
    (β : K → K → ZMod 4) (hβ : IsCocycle kmul β) (a b c : C × K) :
    gtmul i kmul β (gtmul i kmul β a b) c = gtmul i kmul β a (gtmul i kmul β b c) := by
  simp only [gtmul]
  refine Prod.ext ?_ (hk a.2 b.2 c.2)
  have hcoc :
      iPow i (β a.2 b.2) * iPow i (β (kmul a.2 b.2) c.2)
        = iPow i (β b.2 c.2) * iPow i (β a.2 (kmul b.2 c.2)) := by
    rw [← iPow_add i hi, ← iPow_add i hi, hβ]
  linear_combination (a.1 * b.1 * c.1) * hcoc

/-- The Pauli key is an *instance* of the abstract obligation: `tmul` is `gtmul`
at `kmul = componentwise xor` and `β = phaseExp`. -/
theorem tmul_eq_gtmul (a b : Mono C) :
    tmul i a b
      = gtmul i (fun u v => (xor u.1 v.1, xor u.2 v.2))
          (fun u v => phaseExp u.1 u.2 v.1 v.2) a b := rfl

/-- `phaseExp` discharges the abstract 2-cocycle obligation for the Pauli key. -/
theorem phaseExp_isCocycle :
    IsCocycle (K := Bool × Bool) (fun u v => (xor u.1 v.1, xor u.2 v.2))
      (fun u v => phaseExp u.1 u.2 v.1 v.2) :=
  fun u v w => phaseExp_cocycle u.1 u.2 v.1 v.2 w.1 w.2

/-- `tmul_assoc` re-derived from the general theorem, confirming the Pauli case
really is the instantiation and not a coincidence. -/
theorem tmul_assoc_of_gtmul (hi : i ^ 4 = 1) (a b c : Mono C) :
    tmul i (tmul i a b) c = tmul i a (tmul i b c) := by
  simp only [tmul_eq_gtmul]
  exact gtmul_assoc i hi _
    (fun u v w => Prod.ext (Bool.xor_assoc u.1 v.1 w.1) (Bool.xor_assoc u.2 v.2 w.2))
    _ phaseExp_isCocycle a b c

/-! ### L3 is the identity component of L4

`tmul`/`gtmul` above are the *monomial* product; the L4 `Multiply::multiply_into`
on a whole map is the outer product over both supports with collisions
accumulated (`twistedConv` below), i.e. the twisted convolution of
`C[PauliWord]`. L3 `Pair::overlap` is defined independently as the formal
bilinear form `∑ₖ aₖ bₖ` (`PPVM.GradedMap.overlap`), and the design labels it
"the symmetric bilinear Hilbert–Schmidt trace pairing `⟨A,B⟩ = Tr(A B)/2ⁿ`".

Nothing tied the two layers together: L3 and L4 were two unrelated definitions.
`twistedConv_apply_id` closes that — the pairing is exactly the **identity
coefficient of the L4 product**, `⟨A, B⟩ = (A · B)_I`, which is the
operator-algebra definition of the Hilbert–Schmidt pairing in an orthonormal
basis (and, on the Pauli basis, literally `Tr(AB)/2ⁿ`, since only the identity
word has nonzero trace — see `PPVM.PauliMatrix.trace_toOperator_mul` for that
statement against genuine `2ⁿ×2ⁿ` matrices). Two ingredients do the work:
`mulWord_eq_id_iff` collapses the outer product to its diagonal, and
`phaseExpN_self` kills the `iᵏ` twist there, leaving `∑ᵥ aᵥ bᵥ`. Note no
`i⁴ = 1` hypothesis is needed: the surviving phase is `i⁰`. -/

section Convolution

open PPVM.PauliWord PPVM.GradedMap

variable {n : ℕ}

/-- `Multiply::multiply_into` on `C[PauliWord]` — the twisted convolution: the
outer product over both supports, each pair contributing to key `p ⊕ q` with the
`i^{phaseExpN p q}` twist folded onto the coefficient. This is `tmul` (the
monomial case) extended to whole maps. -/
noncomputable def twistedConv (A B : CMap (Word n) C) : CMap (Word n) C :=
  A.sum fun p a => B.sum fun q b =>
    Finsupp.single (mulWord p q) (a * b * iPow i (phaseExpN p q))

/-- **The L3 pairing is the identity coefficient of the L4 product.**
`(A · B) I = ⟪A, B⟫` for the twisted convolution `A · B` and the bilinear
`overlap`. Since `mulWord p q = I ↔ q = p`, the outer product's contribution to
the identity key comes only from the diagonal, where `phaseExpN p p = 0` makes
the twist trivial — what is left is `∑ᵥ aᵥ bᵥ`, which *is* `overlap`. This is
the correctness spec relating `Pair::overlap` to `Multiply`, and it is what
licenses calling `overlap` the Hilbert–Schmidt trace pairing rather than merely a
convenient sum. -/
theorem twistedConv_apply_id (A B : CMap (Word n) C) :
    twistedConv i A B (fun _ => (false, false)) = overlap A B := by
  classical
  simp only [twistedConv, overlap, Finsupp.sum]
  rw [Finset.sum_apply']
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Finset.sum_apply']
  -- Only the diagonal `q = p` survives, and there the phase twist is `i⁰ = 1`.
  have hterm : ∀ q, (Finsupp.single (mulWord p q)
        (A p * B q * iPow i (phaseExpN p q)) : CMap (Word n) C) (fun _ => (false, false))
      = if q = p then A p * B q else 0 := by
    intro q
    rw [Finsupp.single_apply]
    by_cases hq : q = p
    · subst hq
      simp [mulWord_self, phaseExpN_self, iPow]
    · rw [if_neg (fun hc => hq ((mulWord_eq_id_iff p q).mp hc)), if_neg hq]
  rw [Finset.sum_congr rfl fun q _ => hterm q,
    Finset.sum_ite_eq' B.support p fun q => A p * B q]
  by_cases hp : p ∈ B.support
  · rw [if_pos hp]
  · rw [if_neg hp, Finsupp.notMem_support_iff.mp hp, mul_zero]

end Convolution

end PPVM.Twisted
