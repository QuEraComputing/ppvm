/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.DeriveFintype
import Mathlib.Tactic.LinearCombination
import Mathlib.Algebra.MonoidAlgebra.Basic

/-!
# The single-qubit Pauli phase cocycle

This file formalizes the **phase extension** of the single-qubit Pauli algebra —
the `ℤ/4ℤ`-valued cocycle that ppvm's packed representation computes with the
boolean `sign` / `imag` formulas in
`crates/ppvm-pauli-word/src/phase/mul.rs:42`. It is the headline Tier-1 target
of `lean/README.md`.

## The design being validated

The traits-2 design
(`docs/design/traits-2-configuration-and-hashing.md`, "Pauli algebra traits")
factors every Clifford operation into

* the **`Sp` part** — a point of `GF(2)^{2n}`, where multiplication is bit
  addition `⊕` (`PPVM.Pauli.mul`, proved phase-free in `PPVM.Pauli`); and
* the **phase extension** — role-dependent bookkeeping in `ℤ₄` for a phased
  word (`Y = iXZ` needs the `i`).

`KeyProduct::key_mul` (the L4 group-algebra product,
`traits-2-configuration-and-hashing.md#the-map-is-a-graded-algebra-over-ck`) is
specified there as `v · w = ± iᵏ (v ⊕ w)`: the symplectic bits `⊕`, and a phase
`iᵏ` that gets folded onto the coefficient. This file pins down that `k` and
proves it is a genuine group-defining 2-cocycle.

## Conventions

A bare symplectic pair `(x, z) ∈ 𝔽₂²` denotes the physical operator
`g(x,z) = iˣᶻ · Xˣ · Zᶻ` (so `g(1,1) = iXZ = Y`), matching `PPVM.Pauli`. Here we
carry the X/Z bits as `Bool`, because the Rust kernel is bitwise boolean logic
and every statement then closes by `decide` over the finite cases.

## What is proved

* `phaseExp_eq_ref` — the Rust `2·sign + imag` boolean formula equals the
  matrix-model reference `ab + cd + 2bc − (a⊕c)(b⊕d) (mod 4)` derived from
  `g(a,b)·g(c,d) = iˢⁱᵍⁿⁱᵈ g(a⊕c, b⊕d)`. This is the exact refinement claim for
  `phase/mul.rs`.
* `phaseExp_cocycle` — the 2-cocycle identity, i.e. associativity of the twisted
  product. This is what makes the phased Pauli product *well defined*.
* `PhasedPauli` group laws (`mul_assoc'`, `one_mul'`, `mul_one'`,
  `inv_mul_cancel'`) — the single-qubit Pauli group, with multiplication given
  verbatim by the Rust booleans, is a group.
* `phaseExp_sub_comm` — `phaseExp p q − phaseExp q p = 2 · ω(p,q)`: the phase
  asymmetry is exactly twice the symplectic form, i.e. `P·Q = (−1)^{ω} Q·P`.
-/

namespace PPVM.PauliPhase

/-- `𝔽₂` bit as an element of `ℤ/4ℤ` (`0 ↦ 0`, `1 ↦ 1`). Lets the boolean bit
formulas be read as integers mod 4. -/
def b4 (x : Bool) : ZMod 4 := if x then 1 else 0

/-! ### The Rust boolean formulas (`crates/ppvm-pauli-word/src/phase/mul.rs:42`)

With `(a,b)` the left operand's `(x,z)` bits and `(c,d)` the right operand's,
`sign` contributes a `−1 = i²` and `imag` contributes an `i`. -/

/-- `sign` bit of the per-qubit product phase. Mirrors
`(a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d)`. -/
def signBit (a b c d : Bool) : Bool :=
  (a && b && c && !d) || (a && !b && !c && d) || (!a && b && c && d)

/-- `imag` bit of the per-qubit product phase. Mirrors
`(a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d)`. -/
def imagBit (a b c d : Bool) : Bool :=
  (a && !b && d) || (a && !c && d) || (!a && b && c) || (b && c && !d)

/-- The base-`i` exponent (in `ℤ/4ℤ`) of the phase produced by multiplying the
Paulis `g(a,b)` and `g(c,d)`, i.e. `2·sign + imag`. This is the `k` in the
design's `v·w = ± iᵏ (v⊕w)`. -/
def phaseExp (a b c d : Bool) : ZMod 4 :=
  2 * b4 (signBit a b c d) + b4 (imagBit a b c d)

/-! ### The matrix-model reference

Using `X^a Z^b · X^c Z^d = (−1)^{bc} X^{a+c} Z^{b+d}` and the `iˣᶻ`
normalization `g(x,z) = iˣᶻ Xˣ Zᶻ`, one computes

  `g(a,b) · g(c,d) = i^{ab + cd + 2bc} X^{a⊕c} Z^{b⊕d}
                   = i^{ab + cd + 2bc − (a⊕c)(b⊕d)} · g(a⊕c, b⊕d)`,

so the reference exponent is `ab + cd + 2bc − (a⊕c)(b⊕d) (mod 4)`. This
`phaseRef` is an *analytic* formula (not itself a matrix); it is grounded in an
actual `ℤ[i]` matrix model in `PPVM.PauliMatrix` (`pauliMat_mul`), which proves
`phaseExp` is the exponent of the genuine 2×2 matrix product. So the chain is
`phase/mul.rs booleans = phaseExp = phaseRef = real matrix exponent`. -/

/-- Analytic reference exponent (an algebraic formula; `PPVM.PauliMatrix.pauliMat_mul`
grounds it in a genuine `ℤ[i]` matrix product). -/
def phaseRef (a b c d : Bool) : ZMod 4 :=
  b4 a * b4 b + b4 c * b4 d + 2 * (b4 b * b4 c) - b4 (xor a c) * b4 (xor b d)

/-- **Refinement of `phase/mul.rs`.** The packed boolean `2·sign + imag`
computes exactly the analytic reference exponent, for all 16 bit patterns; and
`PPVM.PauliMatrix.pauliMat_mul` checks that reference against real matrices. -/
theorem phaseExp_eq_ref : ∀ a b c d, phaseExp a b c d = phaseRef a b c d := by
  decide

/-! ### The product-bit rule (the `Sp` part) -/

/-- Result X/Z bits of the product: `(a ⊕ c, b ⊕ d)` (`x = a^c`, `z = b^d` in
`phase/mul.rs:45`). -/
def mulBits (a b c d : Bool) : Bool × Bool := (xor a c, xor b d)

/-! ### The 2-cocycle identity (associativity of the twisted product)

For a phased product `(φ, p) · (ψ, q) = (φ + ψ + phaseExp p q, p ⊕ q)` to be
associative, the phase cochain must satisfy the group-cohomology 2-cocycle
identity. Writing the operands as bit pairs `p = (a,b)`, `q = (c,d)`,
`r = (e,f)`: -/

/-- **The 2-cocycle law.** `phaseExp` is a 2-cocycle on `(𝔽₂², ⊕)`, so the phased
Pauli product is associative. -/
theorem phaseExp_cocycle :
    ∀ a b c d e f,
      phaseExp a b c d + phaseExp (xor a c) (xor b d) e f
        = phaseExp c d e f + phaseExp a b (xor c e) (xor d f) := by
  decide

/-- `𝔽₂` bit as an element of `ℤ/2ℤ` (the symplectic-coordinate view of a
`Bool`). -/
def b2 (x : Bool) : ZMod 2 := if x then 1 else 0

@[simp] theorem b2_xor (a b : Bool) : b2 (xor a b) = b2 a + b2 b := by
  cases a <;> cases b <;> decide

/-! ### Small phase facts used to build the group structurally

These are proved on the boolean formula and then feed the `Group` axioms, so the
axioms never `decide` through the bundled `*` (which reduces poorly). -/

/-- Squaring is phase-free: every Pauli is an involution up to sign, `P² = +I`. -/
theorem phaseExp_self (a b : Bool) : phaseExp a b a b = 0 := by
  cases a <;> cases b <;> decide

/-- `I · Q` picks up no phase. -/
theorem phaseExp_id_left (c d : Bool) : phaseExp false false c d = 0 := by
  cases c <;> cases d <;> decide

/-- `P · I` picks up no phase. -/
theorem phaseExp_id_right (a b : Bool) : phaseExp a b false false = 0 := by
  cases a <;> cases b <;> decide

/-! ### The single-qubit phased Pauli group `𝒫₁`

`PhasedPauli` bundles a `ℤ₄` phase with the symplectic X/Z bits. Its
multiplication folds the cocycle onto the phase — exactly what `Sum` does when a
`KeyProduct` "ships the phase to the coefficient." Rather than state loose laws,
we install a genuine Mathlib `Group` instance (axioms by `decide` over the
16-element type), so `𝒫₁` is a first-class group and all of Mathlib's group
theory applies. This is the design's **non-split central extension**
`1 → ℤ₄ → 𝒫₁ → 𝔽₂² → 1` made literal: `toSymplectic` below is the quotient onto
the symplectic bits, and `not_commutative` witnesses that it does not split (so
it is a genuine central extension, **not** a semidirect product `⋉`). The
symplectic *group* `Sp(2n,2)` acts one level up, on this quotient — see
`PPVM.Symplectic`. -/

/-- A single-qubit Pauli with explicit phase: `⟨phase, x, z⟩` denotes
`i^{phase} · g(x, z)`. -/
@[ext]
structure PhasedPauli where
  /-- Base-`i` phase exponent in `ℤ/4ℤ`. -/
  phase : ZMod 4
  /-- X bit. -/
  x : Bool
  /-- Z bit. -/
  z : Bool
deriving DecidableEq, Fintype

namespace PhasedPauli

/-- The twisted product: bits `⊕`, phases add, plus the cocycle. Mirrors the
`MulAssign` in `phase/mul.rs`. -/
def mul (p q : PhasedPauli) : PhasedPauli where
  phase := p.phase + q.phase + phaseExp p.x p.z q.x q.z
  x := xor p.x q.x
  z := xor p.z q.z

/-- The identity `+1 · I`. -/
def one : PhasedPauli := ⟨0, false, false⟩

/-- Inverse: same bits (each Pauli is an `⊕`-involution), negated phase — valid
because every Pauli squares to `+I` (`phaseExp p p = 0`). -/
def inv (p : PhasedPauli) : PhasedPauli := ⟨-p.phase, p.x, p.z⟩

-- The group axioms, proved field-wise on the raw `mul`/`one`/`inv` functions:
-- the symplectic bits by `Bool` `xor` laws, the phase by the 2-cocycle and
-- `ℤ/4ℤ` ring arithmetic. Feeding these to the `Group` instance keeps `decide`
-- (and its slow reduction through the bundled `*`) out of the instance entirely.
theorem mul_assoc' (p q r : PhasedPauli) : mul (mul p q) r = mul p (mul q r) := by
  have h := phaseExp_cocycle p.x p.z q.x q.z r.x r.z
  ext
  · simp only [mul]; linear_combination h
  · simp only [mul, Bool.xor_assoc]
  · simp only [mul, Bool.xor_assoc]

theorem one_mul' (p : PhasedPauli) : mul one p = p := by
  ext <;> simp only [mul, one, phaseExp_id_left, Bool.false_xor, zero_add, add_zero]

theorem mul_one' (p : PhasedPauli) : mul p one = p := by
  ext <;> simp only [mul, one, phaseExp_id_right, Bool.xor_false, add_zero]

theorem inv_mul_cancel' (p : PhasedPauli) : mul (inv p) p = one := by
  ext <;>
    simp only [mul, inv, one, phaseExp_self, Bool.xor_self, neg_add_cancel, add_zero]

/-- **The single-qubit Pauli group.** Its multiplication is exactly the Rust
`phase/mul.rs` booleans; associativity is the 2-cocycle `phaseExp_cocycle`. -/
instance : Group PhasedPauli where
  mul := mul
  one := one
  inv := inv
  mul_assoc := mul_assoc'
  one_mul := one_mul'
  mul_one := mul_one'
  inv_mul_cancel := inv_mul_cancel'

@[simp] theorem mul_def (p q : PhasedPauli) :
    p * q = ⟨p.phase + q.phase + phaseExp p.x p.z q.x q.z, xor p.x q.x, xor p.z q.z⟩ :=
  rfl

@[simp] theorem one_def : (1 : PhasedPauli) = ⟨0, false, false⟩ := rfl

/-- The Pauli group is genuinely noncommutative: `X · Z ≠ Z · X`. Since the
kernel `ℤ₄` and quotient `(ℤ/2)²` of `toSymplectic` are both abelian, a *split*
central extension would force `𝒫₁` abelian; this witness therefore proves the
extension below is **non-split** — a central extension, not a semidirect
product. -/
theorem not_commutative : ∃ p q : PhasedPauli, p * q ≠ q * p :=
  ⟨⟨0, true, false⟩, ⟨0, false, true⟩, by decide⟩

/-- **The non-split central extension `1 → ℤ₄ → 𝒫₁ → 𝔽₂² → 1`.** Forgetting the
phase is a group homomorphism onto the symplectic bits (the additive group
`(ℤ/2)²`, written multiplicatively): the kernel is the central phase group, and
by `not_commutative` the extension does not split. The symplectic *group*
`Sp(2n,2)` acts one level up, on this quotient (`Clifford / Pauli ≅ Sp`) — see
`PPVM.Symplectic`. -/
def toSymplectic : PhasedPauli →* Multiplicative (ZMod 2 × ZMod 2) where
  toFun p := Multiplicative.ofAdd (b2 p.x, b2 p.z)
  map_one' := rfl
  map_mul' p q := by
    simp only [mul_def]
    rw [← ofAdd_add, Prod.mk_add_mk, b2_xor, b2_xor]

/-! ### L4 for the Pauli key: the group algebra

Because `PhasedPauli` is a `Group`, the design's L4 "`KeyProduct` lifts `C[K]`
from a module to an algebra" is Mathlib's `MonoidAlgebra C 𝒫₁` — an associative,
unital `C`-algebra for free. -/

/-- The Pauli group algebra is a genuine associative unital `C`-algebra. -/
noncomputable example (C : Type*) [CommRing C] : Algebra C (MonoidAlgebra C PhasedPauli) :=
  inferInstance

/-- **`key_mul` is the group-algebra product.** Multiplying two basis monomials
multiplies coefficients and takes the Pauli group product of the keys — exactly
`KeyProduct::key_mul`, with the phase carried in the group element (the mod-phase
`PauliSum` keying folds that central phase into the coefficient instead). -/
theorem monoidAlgebra_single_mul {C : Type*} [Semiring C] (p q : PhasedPauli) (a b : C) :
    (MonoidAlgebra.single p a : MonoidAlgebra C PhasedPauli) * MonoidAlgebra.single q b
      = MonoidAlgebra.single (p * q) (a * b) := by
  rw [MonoidAlgebra.single_mul_single]

end PhasedPauli

/-! ### Commutation is the symplectic form

The canonical symplectic form is `PPVM.Symplectic.omega`, a Mathlib
`LinearMap.BilinForm (ℤ/2)` proven alternating there. The version here is its
single-qubit shadow valued in `ℤ/4ℤ` (not `ℤ/2`), because it appears *inside* the
`ℤ/4ℤ` phase equation below: the phase asymmetry of the product is exactly `2·ω`,
i.e. `P·Q = (−1)^{ω(P,Q)} Q·P`. -/

/-- The single-qubit symplectic form as a `ℤ/4ℤ` value (`0` = commute, `1` =
anticommute) — the `ℤ/4ℤ`-valued shadow of `PPVM.Symplectic.omega`, kept here so
`phaseExp_sub_comm` can state the phase law inside `ℤ/4ℤ`. -/
def omega (a b c d : Bool) : ZMod 4 := b4 a * b4 d + b4 b * b4 c

/-- **Commutation ↔ symplectic form.** `phaseExp p q − phaseExp q p = 2·ω(p,q)`;
equivalently `P·Q = (−1)^{ω(P,Q)} · Q·P`. -/
theorem phaseExp_sub_comm :
    ∀ a b c d, phaseExp a b c d - phaseExp c d a b = 2 * omega a b c d := by
  decide

end PPVM.PauliPhase
