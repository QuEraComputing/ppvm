/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-!
# Non-Clifford rotation as a branch on `C[K]`

The design's non-Clifford producer
(`traits-2-configuration-and-hashing.md`, "Every gate is a producer"): a
rotation `extend[s]-by-linearity` a single term into a small fan-out,

  `e^{-iθG/2} · (c · P) · e^{iθG/2} = c·cos θ · P + c·sin θ · (iGP)`

when `G` and `P` anticommute (`{G,P}=0`), and leaves `c·P` unchanged when they
commute. This file validates two things:

* **The new key `P' = iGP` is a distinct Pauli** whose symplectic bits are
  `G ⊕ P` — the same `⊕` key product the L4 algebra uses (`PPVM.PauliPhase`) —
  so the branch produces exactly one extra term, and it collides with nothing
  already indexed by `P`.
* **The branch is a norm-preserving 2-D rotation** on the `(P, P')` coefficient
  plane. The single fact `sin²θ + cos²θ = 1` gives reversibility
  (`R_{-θ} ∘ R_θ = id`), `ℓ²`-norm preservation (the physical *unitarity* of the
  gate), and the angle-addition law — the mathematical guarantees the
  branch-free `eps` / `comm_2` / `SIGN_NEG` kernel (`sum/rot1.rs`, `sum/rot2.rs`)
  must inherit.
-/

namespace PPVM.Rotation

open PPVM.PauliPhase

/-! ### The new key `iGP` is a distinct Pauli with bits `G ⊕ P` -/

/-- **The anticommuting branch produces a genuinely new key.** If `G` and `P`
anticommute (`ω = 1`), then `iGP`, whose symplectic bits are `G ⊕ P`
(`mulBits`), differs from `P`. So a rotation adds exactly one fresh term per
input; it never silently merges back into `P`. -/
theorem anticommute_new_key :
    ∀ a b c d, omega a b c d = 1 → mulBits a b c d ≠ (c, d) := by decide

/-- **The commuting case is inert.** If `G` and `P` commute (`ω = 0`) and `G ≠ I`,
the product bits `G ⊕ P` still differ from `P`; but physically the `sin` branch
carries coefficient `0`, so `c·P` is unchanged. (The bit fact is recorded; the
zeroing is the `sin`-coefficient of `rot 0` below.) -/
theorem commute_bits (a b c d : Bool) :
    mulBits a b c d = (xor a c, xor b d) := rfl

/-! ### The branch is a norm-preserving 2-D rotation

Model the `(coefficient of P, coefficient of P')` pair as a point of `ℝ²`. The
rotation sends `P ↦ cos θ · P + sin θ · P'` and `P' ↦ −sin θ · P + cos θ · P'`,
i.e. the standard rotation matrix. -/

/-- The rotation acting on the `(P, P')` coefficient plane. Starting from a pure
`P` (`v = (1,0)`) gives `(cos θ, sin θ)` — the design's `cos·w + sin·w'`. -/
noncomputable def rot (θ : ℝ) (v : ℝ × ℝ) : ℝ × ℝ :=
  (Real.cos θ * v.1 - Real.sin θ * v.2, Real.sin θ * v.1 + Real.cos θ * v.2)

@[simp] theorem rot_fst (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).1 = Real.cos θ * v.1 - Real.sin θ * v.2 := rfl

@[simp] theorem rot_snd (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).2 = Real.sin θ * v.1 + Real.cos θ * v.2 := rfl

/-- A pure `P` rotates to `cos θ · P + sin θ · P'` — the exact branch shape of the
design's rotation producer. -/
theorem rot_basis (θ : ℝ) : rot θ (1, 0) = (Real.cos θ, Real.sin θ) := by
  simp [rot]

/-- `θ = 0` is the identity (the commuting / no-op branch). -/
theorem rot_zero (v : ℝ × ℝ) : rot 0 v = v := by
  simp [rot]

/-- **Reversibility.** `R_{-θ} ∘ R_θ = id`: rotating by `θ` and then `−θ` returns
the original coefficients. This is the algebraic core of the gate being
invertible, and it uses exactly `sin²θ + cos²θ = 1`. -/
theorem rot_neg_rot (θ : ℝ) (v : ℝ × ℝ) : rot (-θ) (rot θ v) = v := by
  have h := Real.sin_sq_add_cos_sq θ
  refine Prod.ext ?_ ?_ <;>
    simp only [rot_fst, rot_snd, Real.cos_neg, Real.sin_neg]
  · linear_combination v.1 * h
  · linear_combination v.2 * h

/-- **Unitarity.** The rotation preserves the `ℓ²` norm of the coefficient pair,
`‖cos·P + sin·P'‖² = ‖P‖²` — the physical statement that the non-Clifford gate is
unitary, again from `sin²θ + cos²θ = 1`. -/
theorem rot_norm_sq (θ : ℝ) (v : ℝ × ℝ) :
    (rot θ v).1 ^ 2 + (rot θ v).2 ^ 2 = v.1 ^ 2 + v.2 ^ 2 := by
  have h := Real.sin_sq_add_cos_sq θ
  simp only [rot_fst, rot_snd]
  linear_combination (v.1 ^ 2 + v.2 ^ 2) * h

/-- **Angle addition.** Composing two rotations adds the angles, so a Trotter
step of many small rotations about the same axis is one rotation — the identity
a rotation-merging optimization relies on. -/
theorem rot_rot (θ φ : ℝ) (v : ℝ × ℝ) : rot θ (rot φ v) = rot (θ + φ) v := by
  refine Prod.ext ?_ ?_ <;>
    simp only [rot_fst, rot_snd, Real.cos_add, Real.sin_add] <;> ring

end PPVM.Rotation
