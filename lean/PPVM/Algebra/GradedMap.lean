/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.MonoidAlgebra.Basic

/-!
# The coefficient map is the free module `C[K]`

The traits-2 design
(`docs/design/traits-2-configuration-and-hashing.md`, "The map is a graded
algebra over `C[K]`") states that the associative coefficient map — the thing
`HashMap`, `Vec`, and `ColumnStore` all implement — is *algebraically* the

> free `C`-module on a set of keys `K`: a finitely-supported function `K ⇀ C`,
> an element of `C[K]`.

In Mathlib that object is `K →₀ C` (`Finsupp`), and when `K` carries a product
it is `MonoidAlgebra C K` / `AddMonoidAlgebra C K`. This file identifies each
graded layer the design names (L0–L4) with the corresponding Mathlib structure,
so the module and algebra axioms the design relies on are inherited rather than
re-proved, and pins down the two places the design deliberately steps *outside*
the algebra:

* `reduce` (drop zero coefficients) is **structural** in `Finsupp` — the ideal
  model has no zeros in support, which is exactly why the design says
  canonicalization "runs only at finalize."
* `truncate` (drop small coefficients) is **not additive**, so it cannot be an
  algebra operation and must live on `Policy` / `Retain`.

## The layer ↔ Mathlib dictionary

| Design layer | Trait | Mathlib |
| --- | --- | --- |
| L0 finite partial function `K ⇀ C` | `Support` | `Finsupp` (`K →₀ C`) |
| L1 abelian-monoid + canonical support | `Accumulate` | `+` and the `Finsupp` support invariant |
| L2 `C`-module action | `Scale` | `Module C (K →₀ C)` (`•`) |
| L3 bilinear form | `Pair` | the `Finsupp` pairing |
| L4 group-algebra product | `Multiply` | `AddMonoidAlgebra`/`MonoidAlgebra` |
-/

namespace PPVM.GradedMap

open Finsupp

/-- `C[K]` — the free `C`-module on the key set `K`, i.e. finitely-supported
`K ⇀ C`. This is the mathematical object every `Sum` backend refines. -/
abbrev CMap (K C : Type*) [Zero C] := K →₀ C

variable {K C : Type*}

/-! ### L0 `Support` — the finite partial function -/

section Support
variable [Zero C]

/-- `Support::get` — read-out of a coordinate. -/
def get (f : CMap K C) (k : K) : C := f k

/-- `Support::len` — the size of the canonical (zero-free) support. -/
def len (f : CMap K C) : ℕ := f.support.card

/-- **`get` reads the coordinate function**: the map really is `K ⇀ C`. -/
theorem get_eq (f : CMap K C) (k : K) : get f k = f k := rfl

end Support

/-! ### L1 `Accumulate` — form linear combinations, then canonicalize

`accumulate_batch` merges produced terms into the map; algebraically that is
addition in `K →₀ C`, whose `AddCommMonoid` instance is Mathlib's. The design's
separate `reduce()` step drops keys whose coefficient is zero — but a `Finsupp`
*is* its reduced form by construction. -/

section Accumulate
variable [AddCommMonoid C]

/-- `Accumulate::accumulate_batch` (merge two contributions) is `+` on `C[K]`.
(`noncomputable` only because `Finsupp`'s bundled `+` is; the design's backends
implement it concretely.) -/
noncomputable def accumulate (f g : CMap K C) : CMap K C := f + g

@[simp] theorem accumulate_apply (f g : CMap K C) (k : K) :
    accumulate f g k = f k + g k := by
  simp [accumulate]

/-- Accumulation is commutative and associative — the abelian-monoid formation
the design's L1 requires, inherited from `Finsupp`'s `AddCommMonoid`. -/
theorem accumulate_comm (f g : CMap K C) : accumulate f g = accumulate g f :=
  add_comm f g

theorem accumulate_assoc (f g h : CMap K C) :
    accumulate (accumulate f g) h = accumulate f (accumulate g h) :=
  add_assoc f g h

/-- **`reduce()` is structural.** A key is in the canonical support *iff* its
coefficient is nonzero: the ideal `C[K]` model carries no zero coordinates to
drop, so canonicalization is meaningful only for a backend (`HashMap`) that
fails to maintain this invariant. This is the formal content of the design's
"`reduce` runs only at finalize." -/
theorem reduce_structural (f : CMap K C) (k : K) :
    f k ≠ 0 ↔ k ∈ f.support :=
  (mem_support_iff).symm

/-- Consequently a coordinate outside the support is exactly zero: there is
nothing for a correct `reduce` to remove from the ideal model. -/
theorem zero_off_support (f : CMap K C) (k : K) (hk : k ∉ f.support) : f k = 0 := by
  by_contra h
  exact hk (mem_support_iff.2 h)

end Accumulate

/-! ### L2 `Scale` — the `C`-module action -/

section Scale
variable [Semiring C]

/-- `Scale::scale` — multiply every coefficient by `s`; the `C`-module action,
Mathlib's `•` on `K →₀ C`. -/
noncomputable def scale (s : C) (f : CMap K C) : CMap K C := s • f

@[simp] theorem scale_apply (s : C) (f : CMap K C) (k : K) :
    scale s f k = s * f k := by
  simp [scale, Finsupp.smul_apply, smul_eq_mul]

/-- Scaling distributes over accumulation — the module distributive law, so L2
is compatible with L1 exactly as a module action must be. -/
theorem scale_accumulate (s : C) (f g : CMap K C) :
    scale s (f + g) = scale s f + scale s g :=
  smul_add s f g

/-- Scaling is associative with `C`-multiplication (`s·(t·f) = (s·t)·f`). -/
theorem scale_scale (s t : C) (f : CMap K C) :
    scale s (scale t f) = scale (s * t) f := by
  simp [scale, mul_smul]

end Scale

/-! ### L3 `Pair` — the bilinear form (observable read-out)

`Pair::overlap` is the `∑ₖ fₖ · gₖ` read side of the hash join, used for
expectation values. We prove it is **additive in each argument** (biadditive) —
the property the batched probe relies on to split work across produced terms.
(Full `C`-bilinearity additionally needs homogeneity; over a general semiring `C`
biadditivity is the load-bearing part.) -/

section Pair
variable [Semiring C]

/-- `Pair::overlap` — `⟪f, g⟫ = ∑ₖ fₖ gₖ`. Written as a `Finsupp.sum` over `f`'s
support, since terms with `fₖ = 0` contribute nothing. -/
def overlap (f g : CMap K C) : C := f.sum fun k a => a * g k

/-- **Left-additive** — `⟪f₁ + f₂, g⟫ = ⟪f₁, g⟫ + ⟪f₂, g⟫`. -/
theorem overlap_add_left (f₁ f₂ g : CMap K C) :
    overlap (f₁ + f₂) g = overlap f₁ g + overlap f₂ g :=
  Finsupp.sum_add_index' (fun a => zero_mul (g a)) (fun a b₁ b₂ => add_mul b₁ b₂ (g a))

/-- **Right-additive** — `⟪f, g₁ + g₂⟫ = ⟪f, g₁⟫ + ⟪f, g₂⟫`. Together with
`overlap_add_left` this is biadditivity of the pairing. -/
theorem overlap_add_right (f g₁ g₂ : CMap K C) :
    overlap f (g₁ + g₂) = overlap f g₁ + overlap f g₂ := by
  simp only [overlap, Finsupp.add_apply, mul_add]
  exact Finsupp.sum_add

end Pair

/-! ### L4 `Multiply` — the group-algebra product

The design's L4 needs a *`KeyProduct`* on the keys and computes an outer product
over the two operands' supports, accumulating collisions. When the key set is an
additive monoid (the Pauli symplectic bits under `⊕`), that is exactly
`AddMonoidAlgebra` convolution. -/

section Multiply
variable [Semiring C] [AddMonoid K]

/-- `Multiply::multiply_into` — the ring product on `C[K]`, i.e.
`AddMonoidAlgebra` convolution. -/
noncomputable def multiply (f g : AddMonoidAlgebra C K) : AddMonoidAlgebra C K := f * g

/-- **The product of two basis monomials is the design's `key_mul`.**
`(c·|v⟩)·(d·|w⟩) = (c·d)·|v ⊕ w⟩` — one key product `v ⊕ w`, coefficients
multiplied. This is the phase-free core of `KeyProduct::key_mul`; the Pauli case
multiplies the coefficient by the extra `i^{phaseExp}` factor from
`PPVM.PauliPhase` (a 2-cocycle *twist* of this untwisted convolution). -/
theorem multiply_single (v w : K) (c d : C) :
    multiply (Finsupp.single v c) (Finsupp.single w d)
      = Finsupp.single (v + w) (c * d) := by
  simp [multiply, AddMonoidAlgebra.single_mul_single]

end Multiply

/-! ### Truncation is *not* an algebra operation

The design keeps `truncate` off the algebra table and on `Policy`, because
dropping terms still in the support breaks module exactness. Here is the witness:
coefficient-threshold truncation, acting on a single coordinate, is not additive.
`truncMag t` keeps a coefficient only when `|c| ≥ t` — the `CoefficientThreshold`
policy's per-term decision. -/

/-- Coefficient-threshold truncation of one coordinate: keep `c` iff `t ≤ |c|`,
else drop to `0`. This is `CoefficientThreshold::truncate`'s pointwise action. -/
def truncMag (t : ℤ) (c : ℤ) : ℤ := if t ≤ |c| then c else 0

/-- **Truncation breaks additivity.** With threshold `2`, two coefficients that
each round to `0` sum to one that survives: `trunc 1 + trunc 1 = 0 ≠ 2 =
trunc (1 + 1)`. So `truncate` is not a module homomorphism and correctly lives
outside the graded algebra (on `Policy`/`Retain`), not inside it. -/
theorem truncMag_not_additive :
    truncMag 2 1 + truncMag 2 1 ≠ truncMag 2 (1 + 1) := by decide

end PPVM.GradedMap
