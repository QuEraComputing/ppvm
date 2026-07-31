/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.LinearAlgebra.BilinearForm.Properties
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.LinearCombination

/-!
# The symplectic space `GF(2)^{2n}` and the Clifford group as `Sp(2n,2)`

The traits-2 design
(`traits-2-configuration-and-hashing.md`, "Pauli algebra traits"):

> a Pauli operator modulo phase is a vector in the symplectic space
> `GF(2)^{2n}` (the X and Z bit planes), Pauli multiplication is vector addition
> (`⊕`), commutation is the symplectic form `ω(P,Q) = x_P·z_Q ⊕ z_P·x_Q`, and
> conjugation by a Clifford factors into a symplectic map on the bits together
> with a phase. Concretely, the Pauli group is a **non-split central extension**
> `1 → ℤ₄ → 𝒫ₙ → GF(2)^{2n} → 1`; the symplectic group `Sp(2n,2)` enters one level
> up, acting on the quotient (`𝒞ₙ/𝒫ₙ ≅ Sp(2n,2)`).

This file makes the clauses about the quotient literal with Mathlib's linear
algebra (the phase extension itself is `PPVM.PauliPhase`):

* the mod-phase Pauli space `Sp n` is the `ℤ/2`-module `(ℤ/2)² ^ n` (`inferInstance`);
* Pauli multiplication is the module addition (`pauliMul = (· + ·)`, `rfl`);
* `ω` is a genuine `LinearMap.BilinForm (ℤ/2) (Sp n)` that is **alternating**
  (`IsAlt`) — i.e. `(Sp n, ω)` is a symplectic space in Mathlib's sense; and
* the Clifford generators `H`, `S`, `CNOT` are **`ω`-isometries**, i.e. Clifford
  conjugation *lands in* `Sp(2n,2)` acting on the quotient (the `phases` part of
  the extension is `PPVM.PauliPhase.PhasedPauli.toSymplectic`). This proves the
  containment `Clifford conjugation ⊆ Sp`; the surjectivity that would give the
  full isomorphism `𝒞ₙ/𝒫ₙ ≅ Sp(2n,2)` is not formalized here.
-/

namespace PPVM.Symplectic

variable {n : ℕ}

/-- A Pauli word modulo phase: a point of `GF(2)^{2n}`, stored as an `(x, z)` bit
pair per qubit. This is a `ℤ/2`-module by the ambient `Pi`/`Prod` instances. -/
abbrev Sp (n : ℕ) := Fin n → ZMod 2 × ZMod 2

/-- The mod-phase Pauli space really is a `ℤ/2`-vector space. -/
example : Module (ZMod 2) (Sp n) := inferInstance

/-! ### Pauli multiplication is vector addition -/

/-- Pauli multiplication modulo phase. -/
def pauliMul (p q : Sp n) : Sp n := p + q

/-- **"Pauli multiplication is vector addition."** Definitionally the module `+`. -/
theorem pauliMul_eq_add (p q : Sp n) : pauliMul p q = p + q := rfl

/-! ### The symplectic form -/

/-- The bilinear form underlying `ω`: `∑ᵢ x_pᵢ z_qᵢ + z_pᵢ x_qᵢ`. -/
def omegaFun (p q : Sp n) : ZMod 2 :=
  ∑ i, ((p i).1 * (q i).2 + (p i).2 * (q i).1)

/-- **The symplectic form `ω` as a Mathlib bilinear form.** `ω(P,Q) = 0` iff `P`
and `Q` commute; `ω(P,Q) = 1` iff they anticommute. -/
def omega : LinearMap.BilinForm (ZMod 2) (Sp n) :=
  LinearMap.mk₂ (ZMod 2) omegaFun
    (fun p₁ p₂ q => by
      simp only [omegaFun, Pi.add_apply, Prod.fst_add, Prod.snd_add]
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun i _ => by ring)
    (fun c p q => by
      simp only [omegaFun, Pi.smul_apply, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ => by ring)
    (fun p q₁ q₂ => by
      simp only [omegaFun, Pi.add_apply, Prod.fst_add, Prod.snd_add]
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun i _ => by ring)
    (fun c p q => by
      simp only [omegaFun, Pi.smul_apply, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun i _ => by ring)

@[simp] theorem omega_apply (p q : Sp n) : omega p q = omegaFun p q := rfl

/-- **`(Sp n, ω)` is a symplectic space.** The form is alternating: `ω(P,P) = 0`,
because each qubit contributes `x·z + z·x = 2xz = 0` over `𝔽₂`. -/
theorem omega_isAlt : (omega : LinearMap.BilinForm (ZMod 2) (Sp n)).IsAlt := by
  intro v
  simp only [omega_apply, omegaFun]
  refine Finset.sum_eq_zero fun i _ => ?_
  have h : ∀ x y : ZMod 2, x * y + y * x = 0 := by decide
  exact h _ _

/-- The form is symmetric (over `𝔽₂`, alternating ⇒ symmetric). -/
theorem omega_comm (p q : Sp n) : omega p q = omega q p := by
  simp only [omega_apply, omegaFun]
  exact Finset.sum_congr rfl fun i _ => by ring

/-! ### The Clifford generators are `ω`-isometries (the `Sp` part)

Each gate acts on the bit planes; we check it preserves `ω`, i.e. lies in the
symplectic group `Sp(2n,2)`. The gate rules are the design's `SymplecticColumns`
primitives (`swap_xz`, `z ⊕= x`, `x_t ⊕= x_c`, `z_c ⊕= z_t`). -/

/-- `H` on qubit `q`: swap its X and Z bits (`swap_xz`). -/
def hAct (q : Fin n) (v : Sp n) : Sp n := Function.update v q (v q).swap

/-- `S` on qubit `q`: `z ⊕= x`. -/
def sAct (q : Fin n) (v : Sp n) : Sp n :=
  Function.update v q ((v q).1, (v q).1 + (v q).2)

theorem two_zmod2 : (2 : ZMod 2) = 0 := by decide

/-- **`H` is a symplectic isometry.** `ω(H v, H w) = ω(v, w)`, so `H ∈ Sp(2n,2)`.
`H` swaps each qubit's `(x,z)`, and the per-qubit form `xz' + zx'` is invariant
under swapping both operands, so the isometry holds term-by-term. -/
theorem hAct_isometry (q : Fin n) (v w : Sp n) :
    omega (hAct q v) (hAct q w) = omega v w := by
  simp only [omega_apply, omegaFun]
  refine Finset.sum_congr rfl fun i _ => ?_
  rcases eq_or_ne i q with h | h
  · subst h; simp only [hAct, Function.update_self, Prod.fst_swap, Prod.snd_swap]; ring
  · simp only [hAct, Function.update_of_ne h]

/-- **`S` is a symplectic isometry.** Only qubit `q` changes, and its per-qubit
contribution shifts by `2·x·x = 0` over `𝔽₂`, so `S ∈ Sp(2n,2)`. -/
theorem sAct_isometry (q : Fin n) (v w : Sp n) :
    omega (sAct q v) (sAct q w) = omega v w := by
  simp only [omega_apply, omegaFun]
  refine Finset.sum_congr rfl fun i _ => ?_
  rcases eq_or_ne i q with h | h
  · rw [h]
    simp only [sAct, Function.update_self]
    show (v q).1 * ((w q).1 + (w q).2) + ((v q).1 + (v q).2) * (w q).1
        = (v q).1 * (w q).2 + (v q).2 * (w q).1
    linear_combination ((v q).1 * (w q).1) * two_zmod2
  · simp only [sAct, Function.update_of_ne h]

/-! ### Two-qubit gates: a cross-index isometry argument

Unlike `H`/`S`, `CNOT`/`CZ` change two qubits, and the `ω`-contributions of the
control and target qubits each shift; only their *sum* vanishes. So the isometry
is proved by isolating the `{c, t}` terms of the `ω` sum. -/

/-- `CNOT` with control `c`, target `t` (`c ≠ t`): `x_t ⊕= x_c`, `z_c ⊕= z_t`. -/
def cnotAct (c t : Fin n) (v : Sp n) : Sp n :=
  Function.update (Function.update v t ((v c).1 + (v t).1, (v t).2))
    c ((v c).1, (v c).2 + (v t).2)

/-- `CZ` on `c, t` (`c ≠ t`): `z_c ⊕= x_t`, `z_t ⊕= x_c` (symmetric). -/
def czAct (c t : Fin n) (v : Sp n) : Sp n :=
  Function.update (Function.update v t ((v t).1, (v t).2 + (v c).1))
    c ((v c).1, (v c).2 + (v t).1)

variable {c t : Fin n}

@[simp] theorem cnotAct_c (_h : c ≠ t) (v : Sp n) :
    cnotAct c t v c = ((v c).1, (v c).2 + (v t).2) := by
  simp [cnotAct, Function.update_self]

@[simp] theorem cnotAct_t (h : c ≠ t) (v : Sp n) :
    cnotAct c t v t = ((v c).1 + (v t).1, (v t).2) := by
  simp [cnotAct, Function.update_of_ne (Ne.symm h), Function.update_self]

theorem cnotAct_other (v : Sp n) {i} (hc : i ≠ c) (ht : i ≠ t) : cnotAct c t v i = v i := by
  simp [cnotAct, Function.update_of_ne hc, Function.update_of_ne ht]

@[simp] theorem czAct_c (_h : c ≠ t) (v : Sp n) :
    czAct c t v c = ((v c).1, (v c).2 + (v t).1) := by
  simp [czAct, Function.update_self]

@[simp] theorem czAct_t (h : c ≠ t) (v : Sp n) :
    czAct c t v t = ((v t).1, (v t).2 + (v c).1) := by
  simp [czAct, Function.update_of_ne (Ne.symm h), Function.update_self]

theorem czAct_other (v : Sp n) {i} (hc : i ≠ c) (ht : i ≠ t) : czAct c t v i = v i := by
  simp [czAct, Function.update_of_ne hc, Function.update_of_ne ht]

/-- Reduce a preserved `ω` sum to its `{c, t}` terms: if `φ` fixes every qubit
outside `{c, t}` (for both operands), the `ω`-difference is `(diff at c) +
(diff at t)`. -/
private theorem omega_cross (h : c ≠ t) (φ : Sp n → Sp n)
    (hother : ∀ (v : Sp n) i, i ≠ c → i ≠ t → φ v i = v i) (v w : Sp n)
    (hc : (φ v c).1 * (φ w c).2 + (φ v c).2 * (φ w c).1
            + ((φ v t).1 * (φ w t).2 + (φ v t).2 * (φ w t).1)
          = (v c).1 * (w c).2 + (v c).2 * (w c).1
            + ((v t).1 * (w t).2 + (v t).2 * (w t).1)) :
    omega (φ v) (φ w) = omega v w := by
  simp only [omega_apply, omegaFun]
  rw [← sub_eq_zero, ← Finset.sum_sub_distrib,
    ← Finset.sum_subset (Finset.subset_univ ({c, t} : Finset (Fin n)))]
  · rw [Finset.sum_pair h, sub_add_sub_comm, sub_eq_zero]
    exact hc
  · intro i _ hi
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hi
    rw [hother v i hi.1 hi.2, hother w i hi.1 hi.2, sub_self]

/-- **`CNOT` is a symplectic isometry** (`∈ Sp(2n,2)`). -/
theorem cnotAct_isometry (h : c ≠ t) (v w : Sp n) :
    omega (cnotAct c t v) (cnotAct c t w) = omega v w := by
  refine omega_cross h _ (fun u i hc ht => cnotAct_other u hc ht) v w ?_
  simp only [cnotAct_c h, cnotAct_t h]
  linear_combination ((v c).1 * (w t).2 + (v t).2 * (w c).1) * two_zmod2

/-- **`CZ` is a symplectic isometry** (`∈ Sp(2n,2)`). -/
theorem czAct_isometry (h : c ≠ t) (v w : Sp n) :
    omega (czAct c t v) (czAct c t w) = omega v w := by
  refine omega_cross h _ (fun u i hc ht => czAct_other u hc ht) v w ?_
  simp only [czAct_c h, czAct_t h]
  linear_combination ((v c).1 * (w t).1 + (v t).1 * (w c).1) * two_zmod2

end PPVM.Symplectic
