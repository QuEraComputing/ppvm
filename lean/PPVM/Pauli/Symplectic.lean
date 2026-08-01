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

/-! ### The phase-stripped bit maps are bijections (the no-collision invariant)

`ppvm-pauli-sum-2` re-keys every term of a `PauliSum` by the phase-stripped
symplectic bit map and asserts — in `src/producer.rs` and `src/clifford.rs` —
that "A Clifford re-key is a bijection, so colliding re-keyed terms never occur;
`reduce` is a no-op on the support size." That no-collision guarantee is exactly
**injectivity of the word-level bit map** `φ_G : 𝔽₂^{2n} → 𝔽₂^{2n}`: two distinct
keys never land on one, so no coefficient combination/cancellation happens and
the post-`reduce` support size is preserved.

The `*_isometry` theorems above give only `ω`-preservation, which does not by
itself name a bijection; and `conjH_injective`/`conjCNOT_injective`
(`Conjugation.lean`) are injectivity on the *phased* groups `𝒫₁`/`𝒫₂`, not on the
phase-stripped `Sp n` word space the re-key actually keys on. Here we prove each
generator's `Sp n` bit map is an **involution** (`G² = I` on the bits) — note `S`
too: the phase makes `conjS` order 4, but the phase-stripped transvection
`z ⊕= x` squares to the identity over `𝔽₂` — hence a bijection. -/

/-- **`H`'s bit map is an involution**, hence a bijection: swapping `(x,z)` twice
is the identity. -/
theorem hAct_involutive (q : Fin n) : Function.Involutive (hAct q) := by
  intro v; funext i
  by_cases h : i = q
  · subst h; simp [hAct, Function.update_self]
  · simp [hAct, Function.update_of_ne h]

theorem hAct_bijective (q : Fin n) : Function.Bijective (hAct q) :=
  (hAct_involutive q).bijective

/-- **`S`'s bit map is an involution**, hence a bijection: the transvection
`z ⊕= x` applied twice shifts `z` by `2x = 0` over `𝔽₂`. (The *phase* makes
`conjS` order 4; the phase-stripped bit map has order 2.) -/
theorem sAct_involutive (q : Fin n) : Function.Involutive (sAct q) := by
  intro v; funext i
  by_cases h : i = q
  · rw [h]
    refine Prod.ext_iff.mpr ⟨by simp [sAct, Function.update_self], ?_⟩
    simp only [sAct, Function.update_self]
    linear_combination (v q).1 * two_zmod2
  · simp [sAct, Function.update_of_ne h]

theorem sAct_bijective (q : Fin n) : Function.Bijective (sAct q) :=
  (sAct_involutive q).bijective

/-- **`CNOT`'s bit map is an involution** (`c ≠ t`), hence a bijection. -/
theorem cnotAct_involutive (h : c ≠ t) : Function.Involutive (cnotAct c t) := by
  intro v; funext i
  by_cases hic : i = c
  · rw [hic]
    refine Prod.ext_iff.mpr ⟨by simp [cnotAct_c h, cnotAct_t h], ?_⟩
    simp only [cnotAct_c h, cnotAct_t h]
    linear_combination (v t).2 * two_zmod2
  · by_cases hit : i = t
    · rw [hit]
      refine Prod.ext_iff.mpr ⟨?_, by simp [cnotAct_t h, cnotAct_c h]⟩
      simp only [cnotAct_t h, cnotAct_c h]
      linear_combination (v c).1 * two_zmod2
    · simp [cnotAct_other _ hic hit]

theorem cnotAct_bijective (h : c ≠ t) : Function.Bijective (cnotAct c t) :=
  (cnotAct_involutive h).bijective

/-- **`CZ`'s bit map is an involution** (`c ≠ t`), hence a bijection. -/
theorem czAct_involutive (h : c ≠ t) : Function.Involutive (czAct c t) := by
  intro v; funext i
  by_cases hic : i = c
  · rw [hic]
    refine Prod.ext_iff.mpr ⟨by simp [czAct_c h, czAct_t h], ?_⟩
    simp only [czAct_c h, czAct_t h]
    linear_combination (v t).1 * two_zmod2
  · by_cases hit : i = t
    · rw [hit]
      refine Prod.ext_iff.mpr ⟨by simp [czAct_t h, czAct_c h], ?_⟩
      simp only [czAct_t h, czAct_c h]
      linear_combination (v c).1 * two_zmod2
    · simp [czAct_other _ hic hit]

theorem czAct_bijective (h : c ≠ t) : Function.Bijective (czAct c t) :=
  (czAct_involutive h).bijective

/-! ### The loss-guarded action preserves the canonical loss invariant

The lossy word (`crates/ppvm-lossy-pauli-word-2/src/clifford.rs`) differs from
the bare `Sp(2n,2)` action in exactly one way: **a gate touching a lost qubit is
a no-op**. Physically a lost qubit no longer participates; canonically it is held
at X/Z identity (`word-data-structures.md` §"Canonical loss invariant":
`lost[q] ⇒ x[q] = 0 ∧ z[q] = 0`). We model loss by a decidable predicate on
qubits, guard each generator with it, and prove the two properties the crate
relies on but that the phaseless `Sp` action above cannot express:

* **(a)** the guarded action preserves the loss invariant — in particular the
  critical `CNOT` case *present control, lost target* leaves the target at
  `(0,0)` (`clifford.rs` `cnot_present_control_lost_target_preserves_invariant`);
* **(b)** on the present-qubit sub-block (no operand lost) it coincides with the
  already-proven `Sp(2n,2)` isometry, so no symplectic structure is lost where it
  applies.

This is the crate's sole behavioral difference from the bare `PauliWord`, so it
is the one semantic property of the loss guard with a machine-checked counterpart. -/

variable (lost : Fin n → Prop) [DecidablePred lost]

/-- The **canonical loss invariant**: every lost qubit is X/Z identity `(0,0)`. -/
def LossInv (v : Sp n) : Prop := ∀ q, lost q → v q = 0

/-- Loss-guarded `H`: a no-op on a lost qubit, else the bare `hAct`. -/
def hActL (q : Fin n) (v : Sp n) : Sp n := if lost q then v else hAct q v

/-- Loss-guarded `S`: a no-op on a lost qubit, else the bare `sAct`. -/
def sActL (q : Fin n) (v : Sp n) : Sp n := if lost q then v else sAct q v

/-- Loss-guarded `CNOT`: a no-op if **either** operand is lost, else `cnotAct`.
This is the whole-gate skip of `clifford.rs` (`is_lost(ctrl) || is_lost(tgt)`). -/
def cnotActL (c t : Fin n) (v : Sp n) : Sp n :=
  if lost c ∨ lost t then v else cnotAct c t v

/-- Loss-guarded `CZ`: a no-op if either operand is lost, else `czAct`. -/
def czActL (c t : Fin n) (v : Sp n) : Sp n :=
  if lost c ∨ lost t then v else czAct c t v

/-- **(a) `H` preserves the loss invariant.** -/
theorem hActL_preserves_loss (q : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (hActL lost q v) := by
  intro p hp
  by_cases hq : lost q
  · simpa [hActL, hq] using hv p hp
  · have hpq : p ≠ q := fun h => hq (h ▸ hp)
    simpa [hActL, hq, hAct, Function.update_of_ne hpq] using hv p hp

/-- **(a) `S` preserves the loss invariant.** -/
theorem sActL_preserves_loss (q : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (sActL lost q v) := by
  intro p hp
  by_cases hq : lost q
  · simpa [sActL, hq] using hv p hp
  · have hpq : p ≠ q := fun h => hq (h ▸ hp)
    simpa [sActL, hq, sAct, Function.update_of_ne hpq] using hv p hp

/-- **(a) `CNOT` preserves the loss invariant.** The guard skips the whole gate
whenever an operand is lost, so a lost qubit's `(x, z)` never changes. -/
theorem cnotActL_preserves_loss (c t : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (cnotActL lost c t v) := by
  intro p hp
  by_cases hg : lost c ∨ lost t
  · simpa [cnotActL, hg] using hv p hp
  · have hg' := not_or.mp hg
    have hpc : p ≠ c := fun h => hg'.1 (h ▸ hp)
    have hpt : p ≠ t := fun h => hg'.2 (h ▸ hp)
    simpa [cnotActL, hg, cnotAct_other v hpc hpt] using hv p hp

/-- **(a) `CZ` preserves the loss invariant.** -/
theorem czActL_preserves_loss (c t : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (czActL lost c t v) := by
  intro p hp
  by_cases hg : lost c ∨ lost t
  · simpa [czActL, hg] using hv p hp
  · have hg' := not_or.mp hg
    have hpc : p ≠ c := fun h => hg'.1 (h ▸ hp)
    have hpt : p ≠ t := fun h => hg'.2 (h ▸ hp)
    simpa [czActL, hg, czAct_other v hpc hpt] using hv p hp

/-- **The critical `CNOT` case: present control, lost target.** With the target
lost, the guarded gate leaves the target's `(x, z)` at `(0,0)` regardless of the
control — the invariant that a present control must not write onto a lost target
(`clifford.rs` `cnot_present_control_lost_target_preserves_invariant`). -/
theorem cnotActL_lost_target_stays_identity {c t : Fin n} {v : Sp n}
    (hv : LossInv lost v) (ht : lost t) : cnotActL lost c t v t = 0 := by
  simp only [cnotActL, if_pos (Or.inr ht)]
  exact hv t ht

/-- **(b) On present qubits, guarded `H` is the `Sp(2n,2)` isometry.** -/
theorem hActL_present_isometry (q : Fin n) (hq : ¬ lost q) (v w : Sp n) :
    omega (hActL lost q v) (hActL lost q w) = omega v w := by
  simp only [hActL, if_neg hq]; exact hAct_isometry q v w

/-- **(b) On present qubits, guarded `S` is the `Sp(2n,2)` isometry.** -/
theorem sActL_present_isometry (q : Fin n) (hq : ¬ lost q) (v w : Sp n) :
    omega (sActL lost q v) (sActL lost q w) = omega v w := by
  simp only [sActL, if_neg hq]; exact sAct_isometry q v w

/-- **(b) With both operands present, guarded `CNOT` is the `Sp(2n,2)`
isometry.** -/
theorem cnotActL_present_isometry (h : c ≠ t) (hc : ¬ lost c) (ht : ¬ lost t)
    (v w : Sp n) :
    omega (cnotActL lost c t v) (cnotActL lost c t w) = omega v w := by
  simp only [cnotActL, if_neg (not_or.mpr ⟨hc, ht⟩)]; exact cnotAct_isometry h v w

/-- **(b) With both operands present, guarded `CZ` is the `Sp(2n,2)`
isometry.** -/
theorem czActL_present_isometry (h : c ≠ t) (hc : ¬ lost c) (ht : ¬ lost t)
    (v w : Sp n) :
    omega (czActL lost c t v) (czActL lost c t w) = omega v w := by
  simp only [czActL, if_neg (not_or.mpr ⟨hc, ht⟩)]; exact czAct_isometry h v w

/-! ### `CNOT` as two independently loss-guarded column primitives

The blanket `Clifford` in `ppvm-traits-2` does not skip a whole `CNOT` atomically;
it composes two `SymplecticColumns` primitives, each guarded on the *same*
predicate `lost c ∨ lost t` (`crates/ppvm-lossy-pauli-word-2/src/clifford.rs`,
`ppvm-traits-2/src/pauli.rs` `cnot = xor_x_col(c,t); xor_z_col(t,c)`):

* `xor_x_col(c, t)` writes `x_t ⊕= x_c` (`xorXCol`), and
* `xor_z_col(t, c)` writes `z_c ⊕= z_t` (`xorZCol`).

The crate's module doc claims this per-primitive guarding "reproduces the old
whole-gate skip". We machine-check that claim at primitive granularity:
**(i)** each guarded primitive preserves `LossInv` on its own, and **(ii)** their
composition equals the atomic `cnotActL`. Because both primitives test the
identical predicate and neither reads or writes the loss plane, weakening either
guard in isolation would break (ii) — the link the atomic `cnotActL_*` theorems
cannot see. (`CZ` needs no such lemma: the crate emits the single primitive
`cz_bits(a,b)`, which *is* `czAct`, so `czActL` already models the guarded
primitive directly.) -/

/-- `CNOT` column primitive one (`xor_x_col`): `x_t ⊕= x_c`. -/
def xorXCol (c t : Fin n) (v : Sp n) : Sp n :=
  Function.update v t ((v c).1 + (v t).1, (v t).2)

/-- `CNOT` column primitive two (`xor_z_col`): `z_c ⊕= z_t`. -/
def xorZCol (c t : Fin n) (v : Sp n) : Sp n :=
  Function.update v c ((v c).1, (v c).2 + (v t).2)

/-- Loss-guarded `xor_x_col`: no-op if either operand is lost. -/
def xorXColL (c t : Fin n) (v : Sp n) : Sp n :=
  if lost c ∨ lost t then v else xorXCol c t v

/-- Loss-guarded `xor_z_col`: no-op if either operand is lost (same guard). -/
def xorZColL (c t : Fin n) (v : Sp n) : Sp n :=
  if lost c ∨ lost t then v else xorZCol c t v

/-- **(i) `xor_x_col` preserves the loss invariant.** It touches only qubit `t`,
and the guard skips whenever `t` (or `c`) is lost, so a lost qubit is untouched. -/
theorem xorXColL_preserves_loss (c t : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (xorXColL lost c t v) := by
  intro p hp
  by_cases hg : lost c ∨ lost t
  · simpa [xorXColL, hg] using hv p hp
  · have hpt : p ≠ t := fun h => (not_or.mp hg).2 (h ▸ hp)
    simpa [xorXColL, hg, xorXCol, Function.update_of_ne hpt] using hv p hp

/-- **(i) `xor_z_col` preserves the loss invariant.** It touches only qubit `c`,
guarded on the same predicate. -/
theorem xorZColL_preserves_loss (c t : Fin n) {v : Sp n} (hv : LossInv lost v) :
    LossInv lost (xorZColL lost c t v) := by
  intro p hp
  by_cases hg : lost c ∨ lost t
  · simpa [xorZColL, hg] using hv p hp
  · have hpc : p ≠ c := fun h => (not_or.mp hg).1 (h ▸ hp)
    simpa [xorZColL, hg, xorZCol, Function.update_of_ne hpc] using hv p hp

/-- The unguarded columns compose to the atomic `CNOT` bit map (`c ≠ t`): first
`x_t ⊕= x_c`, then `z_c ⊕= z_t` reading the (unchanged) `z_t`. -/
theorem xorZCol_xorXCol_eq_cnotAct (h : c ≠ t) (v : Sp n) :
    xorZCol c t (xorXCol c t v) = cnotAct c t v := by
  simp only [xorZCol, xorXCol, cnotAct, Function.update_of_ne h, Function.update_self]

/-- **(ii) The two guarded columns compose to the atomic whole-gate skip.** When
`c ≠ t`, `xor_z_col ∘ xor_x_col` under the shared guard equals `cnotActL`; a lost
operand short-circuits both primitives to the identity, a present pair runs both.
This is the machine-checked form of the crate's "reproduces the old whole-gate
skip" claim. -/
theorem xorZColL_xorXColL_eq_cnotActL (h : c ≠ t) (v : Sp n) :
    xorZColL lost c t (xorXColL lost c t v) = cnotActL lost c t v := by
  by_cases hg : lost c ∨ lost t
  · simp [xorZColL, xorXColL, cnotActL, hg]
  · simp only [xorZColL, xorXColL, cnotActL, if_neg hg]
    exact xorZCol_xorXCol_eq_cnotAct h v

end PPVM.Symplectic
