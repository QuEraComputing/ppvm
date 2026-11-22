/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Kernels.Model

/-!
The gate rules of `ppvm-pauli-kernels` (`crates/ppvm-pauli-kernels/src/heisenberg.rs`),
transcribed one-to-one as Boolean functions, and proved equal to the
Heisenberg conjugation `G† P G` computed from the gate matrices in `Model.lean`.

A rule maps the pre-gate bits of one Pauli to its post-gate bits and a sign
flip `f`, meaning `G† P G = (−1)^f · P'`. Inputs are listed in the Rust order:
`(x, z)` for one qubit, `(x_a, z_a, x_b, z_b)` for two.
-/

namespace PPVM.Kernels

/-- A one-qubit rule: `(x, z) ↦ ((x', z'), flip)`. -/
abbrev Rule1 := Bool → Bool → (Bool × Bool) × Bool

/-- A two-qubit rule: `(x_a, z_a, x_b, z_b) ↦ ((x_a', z_a', x_b', z_b'), flip)`. -/
abbrev Rule2 := Bool → Bool → Bool → Bool → (Bool × Bool × Bool × Bool) × Bool

/-- `rule` is the Heisenberg action of `g` on every one-qubit Pauli. -/
def Correct1 (g : Gate 2) (rule : Rule1) : Prop :=
  ∀ x z : Bool,
    Mat.Eq (g.heis (pauli x z))
      (Mat.smul (GI.ofInt g.scale * sgn (rule x z).2) (pauli (rule x z).1.1 (rule x z).1.2))

/-- `rule` is the Heisenberg action of `g` on every two-qubit Pauli. -/
def Correct2 (g : Gate 4) (rule : Rule2) : Prop :=
  ∀ xa za xb zb : Bool,
    let r := rule xa za xb zb
    Mat.Eq (g.heis (pauli2 xa za xb zb))
      (Mat.smul (GI.ofInt g.scale * sgn r.2) (pauli2 r.1.1 r.1.2.1 r.1.2.2.1 r.1.2.2.2))

instance (g : Gate 2) (rule : Rule1) : Decidable (Correct1 g rule) := by
  unfold Correct1; infer_instance

instance (g : Gate 4) (rule : Rule2) : Decidable (Correct2 g rule) := by
  unfold Correct2; infer_instance

/-! ## The rules, as in `heisenberg.rs` -/

def ruleX : Rule1 := fun x z => ((x, z), z)
def ruleY : Rule1 := fun x z => ((x, z), x ^^ z)
def ruleZ : Rule1 := fun x z => ((x, z), x)
def ruleH : Rule1 := fun x z => ((z, x), x && z)
def ruleS : Rule1 := fun x z => ((x, z ^^ x), x && !z)
def ruleSdag : Rule1 := fun x z => ((x, z ^^ x), x && z)
def ruleSqrtX : Rule1 := fun x z => ((x ^^ z, z), x && z)
def ruleSqrtXdag : Rule1 := fun x z => ((x ^^ z, z), !x && z)
def ruleSqrtY : Rule1 := fun x z => ((z, x), !x && z)
def ruleSqrtYdag : Rule1 := fun x z => ((z, x), x && !z)

def ruleCnot : Rule2 := fun xc zc xt zt =>
  ((xc, zc ^^ zt, xt ^^ xc, zt), xc && zt && !(xt ^^ zc))
def ruleCz : Rule2 := fun xa za xb zb =>
  ((xa, za ^^ xb, xb, zb ^^ xa), xa && xb && (za ^^ zb))
def ruleCy : Rule2 := fun xc zc xt zt =>
  ((xc, zc ^^ xt ^^ zt, xt ^^ xc, zt ^^ xc), xc && (xt ^^ zt) && !(zc ^^ zt))

/-! ## The gates are unitary after scaling -/

theorem gX_unitary : gX.Unitary := by decide
theorem gY_unitary : gY.Unitary := by decide
theorem gZ_unitary : gZ.Unitary := by decide
theorem gH_unitary : gH.Unitary := by decide
theorem gS_unitary : gS.Unitary := by decide
theorem gSdag_unitary : gSdag.Unitary := by decide
theorem gSqrtX_unitary : gSqrtX.Unitary := by decide
theorem gSqrtXdag_unitary : gSqrtXdag.Unitary := by decide
theorem gSqrtY_unitary : gSqrtY.Unitary := by decide
theorem gSqrtYdag_unitary : gSqrtYdag.Unitary := by decide
theorem gCX_unitary : gCX.Unitary := by decide
theorem gCZ_unitary : gCZ.Unitary := by decide
theorem gCY_unitary : gCY.Unitary := by decide

/-! ## Each rule is its gate's conjugation -/

theorem ruleX_correct : Correct1 gX ruleX := by decide
theorem ruleY_correct : Correct1 gY ruleY := by decide
theorem ruleZ_correct : Correct1 gZ ruleZ := by decide
theorem ruleH_correct : Correct1 gH ruleH := by decide
theorem ruleS_correct : Correct1 gS ruleS := by decide
theorem ruleSdag_correct : Correct1 gSdag ruleSdag := by decide
theorem ruleSqrtX_correct : Correct1 gSqrtX ruleSqrtX := by decide
theorem ruleSqrtXdag_correct : Correct1 gSqrtXdag ruleSqrtXdag := by decide
theorem ruleSqrtY_correct : Correct1 gSqrtY ruleSqrtY := by decide
theorem ruleSqrtYdag_correct : Correct1 gSqrtYdag ruleSqrtYdag := by decide
theorem ruleCnot_correct : Correct2 gCX ruleCnot := by decide
theorem ruleCz_correct : Correct2 gCZ ruleCz := by decide
theorem ruleCy_correct : Correct2 gCY ruleCy := by decide

/-! ## Sanity: the gate matrices are the intended gates

The squares pin `√X` and `√Y` (and their daggers) to the right square roots,
and `S² = Z`, so a wrong matrix cannot slip through as "some other Clifford". -/

theorem gS_sq : Mat.Eq (Mat.mul gS.mat gS.mat) gZ.mat := by decide
theorem gSdag_mul_gS : Mat.Eq (Mat.mul gSdag.mat gS.mat) (Mat.id 2) := by decide
/-- `(2·√X)² = 4·X`, i.e. `√X² = X`. -/
theorem gSqrtX_sq : Mat.Eq (Mat.mul gSqrtX.mat gSqrtX.mat) (Mat.smul (GI.ofInt 4) gX.mat) := by decide
theorem gSqrtXdag_mul_gSqrtX :
    Mat.Eq (Mat.mul gSqrtXdag.mat gSqrtX.mat) (Mat.smul (GI.ofInt 4) (Mat.id 2)) := by decide
/-- `(2·√Y)² = 4·Y`, i.e. `√Y² = Y`. -/
theorem gSqrtY_sq : Mat.Eq (Mat.mul gSqrtY.mat gSqrtY.mat) (Mat.smul (GI.ofInt 4) gY.mat) := by decide
theorem gSqrtYdag_mul_gSqrtY :
    Mat.Eq (Mat.mul gSqrtYdag.mat gSqrtY.mat) (Mat.smul (GI.ofInt 4) (Mat.id 2)) := by decide
/-- `(√2·H)² = 2·I`, i.e. `H² = I`. -/
theorem gH_sq : Mat.Eq (Mat.mul gH.mat gH.mat) (Mat.smul (GI.ofInt 2) (Mat.id 2)) := by decide
/-- `CY = (I ⊗ S) · CX · (I ⊗ S†)`. -/
theorem gCY_eq : Mat.Eq gCY.mat
    (Mat.mul (Mat.kron (Mat.id 2) gS.mat) (Mat.mul gCX.mat (Mat.kron (Mat.id 2) gSdag.mat))) := by
  decide

end PPVM.Kernels
