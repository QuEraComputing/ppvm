/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Symplectic

/-!
# The stabilizer tableau as a symplectic basis

The design's second key type
(`traits-2-configuration-and-hashing.md`, "One engine, two key types" and
`tableau-data-structure.md`): a `Tableau` keys the `TableauMixture = C[Tableau]`
and, internally, is a `2n`-generator stabilizer/destabilizer frame. Its defining
invariant (Aaronson–Gottesman; `tableau-data-structure.md` logical model,
`data.rs:45`) is that

> the `2n` rows form a symplectic basis of the Pauli group,

and *every Clifford preserves it*. This file makes that precise on the symplectic
space `PPVM.Symplectic.Sp` built in `Pauli/Symplectic.lean`:

* `IsSymplecticFrame` — the `ω`-orthonormality relations of a symplectic basis;
* `identityFrame` (destabilizers `Xᵢ`, stabilizers `Zᵢ`) **is** one; and
* **any `ω`-isometry preserves the frame** (`IsSymplecticFrame.map`), so — since
  `H`/`S` are isometries (`Symplectic.hAct_isometry` / `sAct_isometry`) — Clifford
  gates map a symplectic tableau to a symplectic tableau.

The mixture side needs nothing new: a tableau is just a key `K`, so
`TableauMixture = K →₀ C` is the *same* graded algebra `C[K]` of
`PPVM.GradedMap`, keyed on a `Frame` instead of a `PauliWord` — the design's
"one engine, two key types."
-/

namespace PPVM.Tableau

open PPVM.Symplectic

variable {n : ℕ}

/-- A stabilizer/destabilizer tableau (frame) on `n` qubits: `n` destabilizer
generators and `n` stabilizer generators, each a vector of `Sp n = GF(2)^{2n}`. -/
structure Frame (n : ℕ) where
  /-- Destabilizer generators. -/
  destab : Fin n → Sp n
  /-- Stabilizer generators. -/
  stab : Fin n → Sp n

/-- **The symplectic-basis invariant.** Destabilizers pairwise commute,
stabilizers pairwise commute, and destabilizer `i` anticommutes with stabilizer
`j` exactly when `i = j` (`ω(dᵢ, sⱼ) = δᵢⱼ`). These are the defining relations of
a symplectic basis of `(Sp n, ω)`. -/
def IsSymplecticFrame (T : Frame n) : Prop :=
  (∀ i j, omega (T.destab i) (T.destab j) = 0) ∧
  (∀ i j, omega (T.stab i) (T.stab j) = 0) ∧
  (∀ i j, omega (T.destab i) (T.stab j) = if i = j then 1 else 0)

/-- Apply a Pauli-space map `φ` to every generator (Clifford conjugation acts on
the frame this way). -/
def Frame.map (φ : Sp n → Sp n) (T : Frame n) : Frame n where
  destab i := φ (T.destab i)
  stab i := φ (T.stab i)

/-- An `ω`-isometry: preserves the symplectic form (a member of `Sp(2n,2)`). -/
def IsIsometry (φ : Sp n → Sp n) : Prop := ∀ u v, omega (φ u) (φ v) = omega u v

/-- **Every Clifford preserves the symplectic frame.** If `φ` preserves `ω` then
mapping the whole tableau by `φ` keeps all the symplectic-basis relations — the
tableau invariant is a `Clifford`-invariant. -/
theorem IsSymplecticFrame.map {φ : Sp n → Sp n} (hφ : IsIsometry φ) {T : Frame n}
    (hT : IsSymplecticFrame T) : IsSymplecticFrame (T.map φ) := by
  obtain ⟨hdd, hss, hds⟩ := hT
  refine ⟨fun i j => ?_, fun i j => ?_, fun i j => ?_⟩
  · exact (hφ _ _).trans (hdd i j)
  · exact (hφ _ _).trans (hss i j)
  · exact (hφ _ _).trans (hds i j)

/-- `H` on qubit `q` preserves the symplectic frame. -/
theorem isSymplecticFrame_hAct (q : Fin n) {T : Frame n} (hT : IsSymplecticFrame T) :
    IsSymplecticFrame (T.map (hAct q)) :=
  IsSymplecticFrame.map (fun u v => hAct_isometry q u v) hT

/-- `S` on qubit `q` preserves the symplectic frame. -/
theorem isSymplecticFrame_sAct (q : Fin n) {T : Frame n} (hT : IsSymplecticFrame T) :
    IsSymplecticFrame (T.map (sAct q)) :=
  IsSymplecticFrame.map (fun u v => sAct_isometry q u v) hT

/-- **`CNOT` preserves the symplectic frame** — the two-qubit frame update. -/
theorem isSymplecticFrame_cnotAct (c t : Fin n) (h : c ≠ t) {T : Frame n}
    (hT : IsSymplecticFrame T) : IsSymplecticFrame (T.map (cnotAct c t)) :=
  IsSymplecticFrame.map (fun u v => cnotAct_isometry h u v) hT

/-- `CZ` preserves the symplectic frame. -/
theorem isSymplecticFrame_czAct (c t : Fin n) (h : c ≠ t) {T : Frame n}
    (hT : IsSymplecticFrame T) : IsSymplecticFrame (T.map (czAct c t)) :=
  IsSymplecticFrame.map (fun u v => czAct_isometry h u v) hT

/-! ### The initial tableau is a symplectic basis

Destabilizer `i` is `Xᵢ` (bit `x = 1` on qubit `i`), stabilizer `i` is `Zᵢ`
(bit `z = 1`). -/

/-- The initial computational-basis tableau: `dᵢ = Xᵢ`, `sᵢ = Zᵢ`. -/
def identityFrame (n : ℕ) : Frame n where
  destab i := fun k => (if k = i then 1 else 0, 0)
  stab i := fun k => (0, if k = i then 1 else 0)

/-- **The initial tableau is a symplectic basis.** So a stabilizer state's
tableau starts as a symplectic basis, and by `IsSymplecticFrame.map` every
Clifford circuit keeps it one. -/
theorem isSymplecticFrame_identity : IsSymplecticFrame (identityFrame n) := by
  refine ⟨fun i j => ?_, fun i j => ?_, fun i j => ?_⟩ <;>
    simp only [identityFrame, omega_apply, omegaFun]
  · -- destabilizers commute: both Z-planes vanish
    refine Finset.sum_eq_zero fun k _ => ?_
    simp
  · -- stabilizers commute: both X-planes vanish
    refine Finset.sum_eq_zero fun k _ => ?_
    simp
  · -- ω(Xᵢ, Zⱼ) = δᵢⱼ
    simp only [mul_zero, zero_mul, add_zero, ite_mul, one_mul, zero_mul,
      Finset.sum_ite_eq', Finset.mem_univ, if_true]

/-- **The `2n` generators of a symplectic frame are linearly independent** — the
crux of being a symplectic basis of the `2n`-dimensional space `Sp n` (spanning
then follows by a dimension count, which we do not formalize here). If a linear
combination of destabilizers and stabilizers vanishes, every coefficient is `0`:
pairing with `ω(·, sⱼ)` extracts the destabilizer coefficient `aⱼ` (everything
else is `ω`-orthogonal), and pairing with `ω(·, dⱼ)` extracts `bⱼ`. -/
theorem frame_linearIndependent (T : Frame n) (hT : IsSymplecticFrame T)
    (a b : Fin n → ZMod 2)
    (h : ∑ i, a i • T.destab i + ∑ i, b i • T.stab i = 0) :
    (∀ j, a j = 0) ∧ (∀ j, b j = 0) := by
  obtain ⟨hdd, hss, hds⟩ := hT
  have hsd : ∀ i j, omega (T.stab i) (T.destab j) = if j = i then 1 else 0 := by
    intro i j; rw [omega_comm]; exact hds j i
  have expand : ∀ w : Sp n,
      omega (∑ i, a i • T.destab i + ∑ i, b i • T.stab i) w
        = (∑ i, a i • omega (T.destab i) w) + ∑ i, b i • omega (T.stab i) w := by
    intro w
    simp only [map_add, LinearMap.add_apply, map_sum, LinearMap.sum_apply,
      map_smul, LinearMap.smul_apply]
  refine ⟨fun j => ?_, fun j => ?_⟩
  · have hj := expand (T.stab j)
    rw [h, map_zero, LinearMap.zero_apply] at hj
    simp only [hds, hss, smul_eq_mul, mul_ite, mul_one, mul_zero,
      Finset.sum_const_zero, add_zero, Finset.sum_ite_eq', Finset.mem_univ, if_true] at hj
    exact hj.symm
  · have hj := expand (T.destab j)
    rw [h, map_zero, LinearMap.zero_apply] at hj
    simp only [hdd, hsd, smul_eq_mul, mul_ite, mul_one, mul_zero,
      Finset.sum_const_zero, zero_add, Finset.sum_ite_eq, Finset.mem_univ, if_true] at hj
    exact hj.symm

/-! ### Measurement: reading coordinates via `ω`, and the outcome dichotomy

Aaronson–Gottesman measurement of a Pauli `M` inspects how `M` (anti)commutes
with the generators. On the **initial** frame those `ω`-values are literally
`M`'s bits (the coordinate read-outs below); on a general symplectic frame they
are `M`'s coordinates *in the frame basis*, which is a genuine basis by
`frame_linearIndependent`. Whether any stabilizer anticommutes with `M` decides
the deterministic-vs-random branch of the `measure` case split. -/

/-- On the initial frame, `ω(M, Zᵢ)` reads off `M`'s X-bit at qubit `i`. -/
theorem omega_stab_identity (M : Sp n) (i : Fin n) :
    omega M ((identityFrame n).stab i) = (M i).1 := by
  simp only [identityFrame, omega_apply, omegaFun, mul_zero, add_zero, mul_ite, mul_one]
  rw [Finset.sum_ite_eq' Finset.univ i fun k => (M k).1]
  simp

/-- On the initial frame, `ω(M, Xᵢ)` reads off `M`'s Z-bit at qubit `i`. -/
theorem omega_destab_identity (M : Sp n) (i : Fin n) :
    omega M ((identityFrame n).destab i) = (M i).2 := by
  simp only [identityFrame, omega_apply, omegaFun, mul_zero, zero_add, mul_ite, mul_one]
  rw [Finset.sum_ite_eq' Finset.univ i fun k => (M k).2]
  simp

/-- The pivot-search dichotomy (the shape of the `measure` case split): either
`M` commutes with every stabilizer — the pivot search fails and the outcome is
deterministic — or some stabilizer anticommutes with it, giving the pivot for the
random branch. (Structurally this is just the dichotomy "a `𝔽₂`-valued function
is either identically `0` or hits `1` somewhere"; the *content* of each branch is
in the two theorems below.) -/
theorem measurement_dichotomy (T : Frame n) (M : Sp n) :
    (∀ i, omega M (T.stab i) = 0) ∨ ∃ i, omega M (T.stab i) = 1 := by
  by_cases h : ∀ i, omega M (T.stab i) = 0
  · exact Or.inl h
  · rw [not_forall] at h
    obtain ⟨i, hi⟩ := h
    refine Or.inr ⟨i, ?_⟩
    rcases (by decide : ∀ x : ZMod 2, x = 0 ∨ x = 1) (omega M (T.stab i)) with h0 | h1
    · exact absurd h0 hi
    · exact h1

/-- **Deterministic ⇔ `X`-free.** On the initial frame, `M` commutes with every
stabilizer (the measurement is deterministic) *iff* `M` has no `X` component —
i.e. `M` is a `Z`-type (diagonal) Pauli, exactly the operators measured with
certainty on `|0…0⟩`. This gives the deterministic branch real content, via the
coordinate read-out. -/
theorem measure_deterministic_iff_xfree (M : Sp n) :
    (∀ i, omega M ((identityFrame n).stab i) = 0) ↔ ∀ i, (M i).1 = 0 := by
  simp only [omega_stab_identity]

/-- Measuring `Zq` on the initial state is **deterministic**: it commutes with
every stabilizer, so the pivot search fails. -/
theorem measure_Z_deterministic (q i : Fin n) :
    omega ((identityFrame n).stab q) ((identityFrame n).stab i) = 0 :=
  isSymplecticFrame_identity.2.1 q i

/-- Measuring `Xq` is **random**: it anticommutes with stabilizer `q` (the pivot),
so the frame update fires. -/
theorem measure_X_pivot (q : Fin n) :
    omega ((identityFrame n).destab q) ((identityFrame n).stab q) = 1 := by
  simpa using isSymplecticFrame_identity.2.2 q q

end PPVM.Tableau
