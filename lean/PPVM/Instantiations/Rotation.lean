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
commute. This file validates two things about that branch — **at the level of
the bits and the coefficient plane**, not by deriving the operator conjugation:

* **The new key `P' = iGP` is a distinct Pauli** whose symplectic bits are
  `G ⊕ P` — the same `⊕` key product the L4 algebra uses (`PPVM.PauliPhase`) —
  so the branch produces exactly one extra term, distinct from both `G` and `P`.
  This part *is* derived (from `mulBits` and the anticommutation `ω = 1`).
* **The coefficient update `(c_P, c_{P'})` is a 2-D rotation.** We formalize the
  standard rotation on `ℝ²` and prove it is norm-preserving, reversible, and
  angle-additive (all from `sin²θ + cos²θ = 1`). This is the abstract shape the
  branch applies to the coefficient pair; the operator identity
  `e^{-iθG/2}(cP)e^{iθG/2} = …` is *modeled* by this 2-D rotation, **not** derived
  from operator algebra here. It gives the algebraic guarantees the branch-free
  `eps` / `comm_2` / `SIGN_NEG` kernel (`sum/rot1.rs`, `sum/rot2.rs`) inherits.
-/

namespace PPVM.Rotation

open PPVM.PauliPhase

/-! ### The new key `iGP` is a distinct Pauli with bits `G ⊕ P` -/

/-- **The anticommuting branch produces a genuinely new key, distinct from both
operands.** If `G = (a,b)` and `P = (c,d)` anticommute (`ω = 1`), then `iGP`,
whose symplectic bits are `G ⊕ P` (`mulBits`), differs from `P` *and* from `G`.
So a rotation adds exactly one fresh term per input; it never silently merges
back into `P` or into `G`. -/
theorem anticommute_new_key :
    ∀ a b c d, omega a b c d = 1 →
      mulBits a b c d ≠ (c, d) ∧ mulBits a b c d ≠ (a, b) := by decide

/-- The product bits are `(a⊕c, b⊕d)` — a definitional restatement of `mulBits`.
(In the *commuting* case the rotation is inert not because of these bits but
because the `sin` branch carries coefficient `0`; `rot 0` below is the identity.) -/
theorem commute_bits (a b c d : Bool) :
    mulBits a b c d = (xor a c, xor b d) := rfl

/-! ### The branch sign `ε` is the real phase of `iGP`, derived from `phaseExp`

`Rotation.lean`'s 2-D `rot` *models* the coefficient plane but does not derive the
per-axis `±1` sign the branch coefficient carries. That sign is a genuine Pauli
phase fact: the stored branch key is the **real** word `g(G⊕P)`, and the physical
branch term is `sinθ · (iGP)`, with

  `iGP = i · (G·P) = i · i^{phaseExp(G,P)} · g(G⊕P) = i^{1 + phaseExp(G,P)} · g(G⊕P)`.

So the coefficient sign is `ε = i^{1 + phaseExp(G,P)}` as a base-`i` exponent
(`branchExp` below). When `G` and `P` anticommute (`ω = 1`) the single-qubit
product `G·P` carries an odd power of `i`, so `1 + phaseExp` is *even* and `ε` is a
real `±1` — this is why the leading `i` "cancels" and the branch stays real. The
three axis theorems then check that this derived `±1` equals the hand-ported table
in `producer.rs:141-143` (`RotationProducer::produce`) case-by-case. -/

/-- The base-`i` exponent of the prefactor sitting on the *real* branch word
`g(G⊕P)`: `iGP = i^{1 + phaseExp(G,P)} · g(G⊕P)`. -/
def branchExp (gx gz x z : Bool) : ZMod 4 := 1 + phaseExp gx gz x z

/-- **The `i` in `iGP` cancels the anticommuting product's `i`.** When `G = (gx,gz)`
and `P = (x,z)` anticommute (`ω = 1`), `branchExp` is even (`∈ {0,2}`), i.e. the
branch prefactor is a real `±1`, never `±i` — so `producer.rs` may store a real
coefficient sign `ε` on the real word `g(G⊕P)`. -/
theorem branchExp_isRealPhase (gx gz x z : Bool) (h : omega gx gz x z = 1) :
    IsRealPhase (branchExp gx gz x z) := by
  change branchExp gx gz x z = 0 ∨ branchExp gx gz x z = 2
  revert h; revert gx gz x z; decide

/-- **`rx` (`G = X = (1,0)`): the table `ε = −1 iff x`.** For an anticommuting `P`
(`ω = 1`, i.e. `z = 1`), the phase-derived branch sign `i^{1+phaseExp(X,P)}` equals
`−1` (exponent `2`) exactly when `x`, matching `producer.rs:141`
(`RotAxis::X => … eps = if x { -1 } else { 1 }`). -/
theorem rx_eps_from_product :
    ∀ x z, omega true false x z = 1 →
      branchExp true false x z = (if x then 2 else 0) := by decide

/-- **`ry` (`G = Y = (1,1)`): the table `ε = −1 iff z`.** For an anticommuting `P`
(`ω = 1`, i.e. `x ≠ z`), `i^{1+phaseExp(Y,P)}` equals `−1` exactly when `z`,
matching `producer.rs:142` (`RotAxis::Y => … eps = if z { -1 } else { 1 }`). -/
theorem ry_eps_from_product :
    ∀ x z, omega true true x z = 1 →
      branchExp true true x z = (if z then 2 else 0) := by decide

/-- **`rz` (`G = Z = (0,1)`): the table `ε = +1 iff z`.** For an anticommuting `P`
(`ω = 1`, i.e. `x = 1`), `i^{1+phaseExp(Z,P)}` equals `+1` (exponent `0`) exactly
when `z`, matching `producer.rs:143`
(`RotAxis::Z => … eps = if z { 1 } else { -1 }`). -/
theorem rz_eps_from_product :
    ∀ x z, omega false true x z = 1 →
      branchExp false true x z = (if z then 0 else 2) := by decide

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

/-- **Norm preservation.** The 2-D rotation preserves the `ℓ²` norm of the
coefficient pair, `‖cos·P + sin·P'‖² = ‖P‖²`, from `sin²θ + cos²θ = 1`. This is
the algebraic reason the branch is norm-preserving; the operator-level unitarity
of `e^{-iθG/2}` is *modeled* by this, not derived from operator algebra. -/
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
