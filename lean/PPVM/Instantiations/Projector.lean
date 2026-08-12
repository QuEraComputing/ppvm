/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Matrix

/-!
# `Projection::p0` / `p1` on `C[K]`: linearity, idempotence, and the old kernel

`crates/ppvm-pauli-sum-2/src/proj.rs` is the one shipped public gate on `Sum`
whose semantics had **no** oracle: `grep -rl proj lean/PPVM/` found only
`Tableau/Projection.lean`, which is about the *generalized tableau's* amplitude
vector (`Amp n = Bitstring n → ℂ`, the case-a/case-b measurement arms) and says
nothing about the Heisenberg action of a computational-basis projector on a
Pauli-keyed `C[K]`. No old test or benchmark exercises `p0`/`p1` either, so old's
kernel was never checked against anything.

## What the crate computes

Old (`ppvm-pauli-sum/src/sum/proj.rs:16-30`, ported bit-for-bit into
`ppvm-pauli-sum-2/src/proj.rs:57-81`) reads

```text
let half = v.half();          -- `half == v/2`, a *value*, not the ring's ½
match k.get(pos) {
  Pauli::I => { *v *= half; Some((k.set_new(pos, Z), v.clone())) }
  Pauli::Z => { *v *= half; Some((k.set_new(pos, I), v.clone())) }
  _ => None,                  -- X / Y: coefficient untouched, no branch
}
```

so the survivor **and** the branch both come out at `c²/2`, not `c/2`
(`oldProj`). This file is the oracle that adjudicates it.

## What is proved

* `projLin_add`, `projLin_smul` — the intended map (halving by the ring constant
  `½`) is **linear** on `C[K]`.
* `projLin_idem` — and **idempotent**, `P_b ∘ P_b = P_b`, for either branch sign.
* `oldProj_ne_projLin`, `oldProj_not_additive`, `oldProj_not_idem` — old's
  `c ↦ c²/2` is neither linear nor idempotent, and
  `oldStep_eq_half_iff` pins exactly when it coincides with the correct map:
  `c²/2 = c/2 ↔ c ∈ {0, 1}`. Unit-coefficient stabilizer sums — old's only usage —
  are precisely the blind spot. **Old is wrong; the Lean-correct value is `c/2`.**
* `twoProj_conj_I` … `twoProj_conj_Y` — the genuine operator-level Heisenberg
  action `A ↦ Π A Π` with `Π = (I+Z)/2`, over honest `ℤ[i]` matrices. It agrees
  with `projLin` on the `I`/`Z` block and **kills** `X`/`Y` — whereas the crate's
  `_ => None` leaves them untouched. `projLin_p0_add_p1` exhibits the observable
  consequence: `p0 + p1` is the identity on `I`/`Z` but *doubles* `X`/`Y`, where
  completeness `Π₀ + Π₁ = 1` forces the dephasing channel (`X, Y ↦ 0`). That is a
  second, independent defect in old's kernel, recorded here rather than silently
  patched.
-/

namespace PPVM.Projector

open PPVM.PauliMatrix

/-! ### The site alphabet

One qubit's Pauli as its `(x, z)` bits, and the coefficient vector of a `C[K]`
restricted to that site — the four coordinates `I = (0,0)`, `X = (1,0)`,
`Z = (0,1)`, `Y = (1,1)` the gate reads. -/

/-- A single-qubit Pauli's `(x, z)` bits, `g(x,z) = iˣᶻ Xˣ Zᶻ` as everywhere else. -/
abbrev Bits := Bool × Bool

/-- The coefficient vector of a `C[K]` at one site (the four coordinates the gate
touches). Addition and scaling are pointwise, i.e. the free-module structure of
`PPVM.GradedMap.CMap` restricted to the site. -/
abbrev Coeffs := Bits → ℝ

/-- The branch key: `p0`/`p1` toggle the `z` bit at the qubit (`I ↔ Z`), leaving
`x` alone — `k.toggled_bits(qubit, false, true)`. -/
def branchKey (p : Bits) : Bits := (p.1, !p.2)

@[simp] theorem branchKey_involutive (p : Bits) : branchKey (branchKey p) = p := by
  simp only [branchKey, Bool.not_not]

/-! ### The intended (linear) gate -/

/-- **The projector gate as the crate intends it**: halve by the *ring constant*
`½` and add the `Z`-toggled partner with sign `ε` (`ε = +1` for `p0`, `ε = −1`
for `p1`); an `X`/`Y` key (`x = true`) is left alone, which is old's `_ => None`.

Read as a map on coefficient vectors: the coefficient landing at key `p` is the
halved survivor `A p` plus the halved branch arriving from `branchKey p`. -/
noncomputable def projLin (ε : ℝ) (A : Coeffs) : Coeffs :=
  fun p => if p.1 then A p else (A p + ε * A (branchKey p)) / 2

/-- **`p0`/`p1` are additive.** -/
theorem projLin_add (ε : ℝ) (A B : Coeffs) :
    projLin ε (fun p => A p + B p) = fun p => projLin ε A p + projLin ε B p := by
  funext p
  cases hp : p.1
  · simp only [projLin, hp, Bool.false_eq_true, if_false]; ring
  · simp only [projLin, hp, if_true]

/-- **`p0`/`p1` are homogeneous** — the halving factor is the ring's `½` and does
not depend on the coefficient, so scaling the input scales the output. This is
exactly the law old's `let half = v.half()` breaks. -/
theorem projLin_smul (ε c : ℝ) (A : Coeffs) :
    projLin ε (fun p => c * A p) = fun p => c * projLin ε A p := by
  funext p
  cases hp : p.1
  · simp only [projLin, hp, Bool.false_eq_true, if_false]; ring
  · simp only [projLin, hp, if_true]

/-- **`p0`/`p1` are idempotent**, `P_b ∘ P_b = P_b`, for either branch sign
(`ε² = 1`). Together with `projLin_add`/`projLin_smul` this is the defining
property of a projector: the gate is a genuine idempotent linear operator on
`C[K]`. -/
theorem projLin_idem (ε : ℝ) (hε : ε * ε = 1) (A : Coeffs) :
    projLin ε (projLin ε A) = projLin ε A := by
  funext p
  by_cases hp : p.1 = true
  · simp only [projLin, hp, if_true]
  · have hb : (branchKey p).1 = p.1 := rfl
    simp only [projLin, hp, Bool.false_eq_true, if_false, hb, branchKey_involutive]
    linear_combination (A p / 4) * hε

/-! ### Old's kernel: `c ↦ c²/2` -/

/-- Old's per-term coefficient update, verbatim: `let half = v.half(); *v *= half`
leaves `v = c·(c/2) = c²/2`, and the branch is pushed from the *already mutated*
`v`. -/
noncomputable def oldStep (c : ℝ) : ℝ := c * (c / 2)

/-- **Old's halving coincides with the correct one exactly on `{0, 1}`.**
`c²/2 = c/2 ↔ c = 0 ∨ c = 1` — so old's only usage, unit-coefficient stabilizer
sums, is precisely the blind spot that hid the defect. -/
theorem oldStep_eq_half_iff (c : ℝ) : oldStep c = c / 2 ↔ c = 0 ∨ c = 1 := by
  constructor
  · intro h
    have h0 : c * (c - 1) = 0 := by
      simp only [oldStep] at h; linarith [h]
    rcases mul_eq_zero.mp h0 with h1 | h1
    · exact Or.inl h1
    · exact Or.inr (by linarith)
  · rintro (h | h) <;> subst h <;> simp [oldStep]

/-- **Old's gate is not homogeneous.** Halving by the term's own value is
quadratic, so doubling the input quadruples — not doubles — the output. -/
theorem oldStep_not_homogeneous : oldStep (2 * 1) ≠ 2 * oldStep 1 := by
  simp only [oldStep]; norm_num

/-- **Old's gate is not additive.** `(1+1)²/2 = 2 ≠ 1 = 1²/2 + 1²/2`. -/
theorem oldStep_not_additive : oldStep (1 + 1) ≠ oldStep 1 + oldStep 1 := by
  simp only [oldStep]; norm_num

/-- Old's `p0` as a whole-map operator: survivor and branch both at `c²/2`. -/
noncomputable def oldProj (ε : ℝ) (A : Coeffs) : Coeffs :=
  fun p => if p.1 then A p else oldStep (A p) + ε * oldStep (A (branchKey p))

/-- **Old's gate is not idempotent.** On the single term `2·I` the correct
projector gives `1·I + 1·Z` and is stable under a second application, while old
gives `2·I + 2·Z` and then `4·I + 4·Z`: the "projector" *grows* the state. -/
theorem oldProj_not_idem :
    oldProj 1 (oldProj 1 (fun p => if p = (false, false) then (2 : ℝ) else 0)) (false, false)
      ≠ oldProj 1 (fun p => if p = (false, false) then (2 : ℝ) else 0) (false, false) := by
  norm_num [oldProj, oldStep, branchKey]

/-- **Old differs from the linear projector already at `c = 2`.** The correct
value on the survivor is `c/2 = 1`; old computes `c²/2 = 2`. -/
theorem oldProj_ne_projLin :
    oldProj 1 (fun p => if p = (false, false) then (2 : ℝ) else 0) (false, false)
      ≠ projLin 1 (fun p => if p = (false, false) then (2 : ℝ) else 0) (false, false) := by
  norm_num [oldProj, projLin, oldStep, branchKey]

/-! ### The operator-level projector, and old's second defect (`X`/`Y`)

`projLin` above is the *crate's intended* map. Independently of the halving bug,
is that map the Heisenberg action of `|0⟩⟨0|`? Over honest `ℤ[i]` matrices, with
`2Π = I + Z` (so the statements below are `4·Π P Π`):

* `Π I Π = Π Z Π = ½(I + Z)` — matching `projLin`; but
* `Π X Π = Π Y Π = 0` — where the crate's `_ => None` leaves `X`/`Y` **untouched**.
-/

/-- `2Π = I + Z` — the unnormalized computational-basis projector, kept in `ℤ[i]`
(which has no `½`) so every statement below closes by `decide`. -/
def twoProj : Matrix (Fin 2) (Fin 2) GaussianInt :=
  pauliMat false false + pauliMat false true

/-- `(2Π) I (2Π) = 2·(2Π)`, i.e. `Π I Π = ½(I + Z)` — the `I` column of `projLin`. -/
theorem twoProj_conj_I :
    twoProj * pauliMat false false * twoProj = twoProj + twoProj := by decide

/-- `(2Π) Z (2Π) = 2·(2Π)`, i.e. `Π Z Π = ½(Z + I)` — the `Z` column of `projLin`. -/
theorem twoProj_conj_Z :
    twoProj * pauliMat false true * twoProj = twoProj + twoProj := by decide

/-- **`Π X Π = 0`** — the projector *annihilates* `X`, it does not fix it. -/
theorem twoProj_conj_X : twoProj * pauliMat true false * twoProj = 0 := by decide

/-- **`Π Y Π = 0`** — likewise for `Y`. -/
theorem twoProj_conj_Y : twoProj * pauliMat true true * twoProj = 0 := by decide

/-- **The observable consequence of old's `_ => None`.** Completeness
`Π₀ + Π₁ = 1` makes `A ↦ Π₀ A Π₀ + Π₁ A Π₁` the dephasing channel, which sends
`X, Y ↦ 0` (`twoProj_conj_X`/`_Y` plus the mirrored `(I−Z)/2` arm). The crate's
`p0 + p1` instead returns `A` on the `I`/`Z` block — correct — and `2A` on
`X`/`Y`. So even after the halving is fixed, the `X`/`Y` arm is a separate
divergence from the projector semantics. -/
theorem projLin_p0_add_p1 (A : Coeffs) (p : Bits) :
    projLin 1 A p + projLin (-1) A p = if p.1 then 2 * A p else A p := by
  cases hp : p.1
  · simp only [projLin, hp, Bool.false_eq_true, if_false]; ring
  · simp only [projLin, hp, if_true]; ring

end PPVM.Projector
