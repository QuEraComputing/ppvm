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
follows by the dimension count `frame_span`, and the resulting coordinate formula
is `frame_coordinate_expansion`). If a linear
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

/-! ### The measurement projection preserves the symplectic frame

`IsSymplecticFrame.map` covers only the **unitary** generators: it requires an
`ω`-isometry, which the Aaronson–Gottesman *measurement projection* is not (it is
not even injective on `Sp n` — it overwrites two rows). Yet the projection is the
only non-unitary frame mutation in the crate
(`crates/ppvm-tableau-2/src/data.rs`, `update_tableau_according_to_outcome`), and
`StabilizerFrame::canonicalize` (`crates/ppvm-tableau-2/src/clifford.rs`) is a
deliberate **no-op** justified precisely by the claim that the projection
"restores the destabilizer/stabilizer pairing in place". Every downstream
`compute_decomposition` (Yoder-2012 Lemma 5) assumes a valid symplectic basis, so
that claim carries the whole component; it is proved here.

With measured qubit `q` and pivot `p` (a stabilizer with `ω(sₚ, Z_q) = 1`, i.e.
`sₚ` has its `x`-bit set at `q` — the `stabilizers[i].xbits[addr0]` test) the
projection is

* every **other** row (destabilizer or stabilizer) whose `x`-bit at `q` is set
  gets `sₚ` multiplied in;
* `dₚ := sₚ` (the old stabilizer becomes the new destabilizer);
* `sₚ := (−1)^b Z_q`.

The outcome bit `b` is a *phase*, and `Sp n` is the phase-stripped space, so `b`
does not appear below — the frame relations are insensitive to it. -/

/-- `Z_q` as a point of `Sp n`: the row the projection installs as the new
stabilizer `p` (the `(−1)^b` sign is a phase, invisible here). -/
def zVec (q : Fin n) : Sp n := fun k => (0, if k = q then 1 else 0)

/-- `ω(v, Z_q)` reads off `v`'s X-bit at `q` — the crate's
`stab.anticommutes_at(addr0, (false, true))` / `row.xbits[addr0]` test. -/
theorem omega_zVec (v : Sp n) (q : Fin n) : omega v (zVec q) = (v q).1 :=
  omega_stab_identity v q

private theorem omega_addSmul_left (g v w : Sp n) (a : ZMod 2) :
    omega (v + a • g) w = omega v w + a * omega g w := by
  simp only [map_add, LinearMap.add_apply, map_smul, LinearMap.smul_apply, smul_eq_mul]

private theorem omega_addSmul_right (g v w : Sp n) (b : ZMod 2) :
    omega v (w + b • g) = omega v w + b * omega v g := by
  simp only [map_add, map_smul, smul_eq_mul]

private theorem omega_addSmul (g v w : Sp n) (a b : ZMod 2) :
    omega (v + a • g) (w + b • g) = omega v w + b * omega v g + a * omega g w := by
  rw [omega_addSmul_left, omega_addSmul_right, omega_addSmul_right, omega_isAlt g,
    mul_zero, add_zero]

/-- The conditional row multiply the projection performs on every non-pivot row:
multiply the pivot generator `g` in **iff** the row's `x`-bit at the measured
qubit is set. Written with the `𝔽₂` scalar `(v q).1` so the case split is
algebraic. -/
def rowUpdate (g : Sp n) (q : Fin n) (v : Sp n) : Sp n := v + (v q).1 • g

/-- `rowUpdate` really is the crate's `if row.xbits[addr0] { row.mul_assign(&g_q) }`. -/
theorem rowUpdate_eq_ite (g : Sp n) (q : Fin n) (v : Sp n) :
    rowUpdate g q v = if (v q).1 = 1 then v + g else v := by
  simp only [rowUpdate]
  rcases (by decide : ∀ x : ZMod 2, x = 0 ∨ x = 1) ((v q).1) with h | h <;> rw [h]
  · rw [zero_smul, add_zero, if_neg (by decide)]
  · rw [one_smul, if_pos rfl]

/-- **The Aaronson–Gottesman measurement projection on the frame**, for measured
qubit `q` and pivot `p` (`update_tableau_according_to_outcome`). -/
def projectFrame (T : Frame n) (q p : Fin n) : Frame n where
  destab i := if i = p then T.stab p else rowUpdate (T.stab p) q (T.destab i)
  stab i := if i = p then zVec q else rowUpdate (T.stab p) q (T.stab i)

/-- **The measurement projection preserves the symplectic-frame invariant.**

The one non-unitary frame mutation in the crate maps a symplectic basis to a
symplectic basis, so `canonicalize` genuinely has nothing to do and every
subsequent `compute_decomposition` runs on a valid basis. The pivot hypothesis
`ω(sₚ, Z_q) = 1` is exactly what `find_z_anticommuting_stabilizer` returns.

Each of the nine cases is a two-term bilinear expansion; the two interesting ones
are the new stabilizer row `Z_q` against an updated row (where the `𝔽₂` identity
`x + x = 0` kills the cross term, using `ω(sₚ, Z_q) = 1`) and the new
destabilizer/stabilizer pair `(sₚ, Z_q)`, whose pairing is the pivot hypothesis
itself. -/
theorem isSymplecticFrame_projectFrame {T : Frame n} (hT : IsSymplecticFrame T)
    (q p : Fin n) (hpiv : omega (T.stab p) (zVec q) = 1) :
    IsSymplecticFrame (projectFrame T q p) := by
  obtain ⟨hdd, hss, hds⟩ := hT
  have hself : ∀ x : ZMod 2, x + x = 0 := by decide
  have hgz : omega (T.stab p) (zVec q) = 1 := hpiv
  have hzg : omega (zVec q) (T.stab p) = 1 := by rw [omega_comm]; exact hpiv
  refine ⟨fun i j => ?_, fun i j => ?_, fun i j => ?_⟩
  · -- destabilizers pairwise commute
    by_cases hi : i = p <;> by_cases hj : j = p <;>
      simp only [projectFrame, rowUpdate, if_pos, hi, hj, if_false]
    · exact hss p p
    · rw [omega_addSmul_right, hss p p, mul_zero, add_zero, omega_comm]
      simpa [hj] using hds j p
    · rw [omega_addSmul_left, omega_isAlt (T.stab p), mul_zero, add_zero]
      simpa [hi] using hds i p
    · rw [omega_addSmul, hdd i j,
        show omega (T.destab i) (T.stab p) = 0 by simpa [hi] using hds i p,
        show omega (T.stab p) (T.destab j) = 0 by
          rw [omega_comm]; simpa [hj] using hds j p]
      simp
  · -- stabilizers pairwise commute
    by_cases hi : i = p <;> by_cases hj : j = p <;>
      simp only [projectFrame, rowUpdate, if_pos, hi, hj, if_false]
    · exact omega_isAlt (zVec q)
    · rw [omega_addSmul_right, hzg, mul_one, omega_comm, omega_zVec]
      exact hself _
    · rw [omega_addSmul_left, hgz, mul_one, omega_zVec]
      exact hself _
    · rw [omega_addSmul, hss i j, hss i p, hss p j]
      simp
  · -- ω(dᵢ, sⱼ) = δᵢⱼ
    by_cases hi : i = p <;> by_cases hj : j = p <;>
      simp only [projectFrame, rowUpdate, if_pos, hi, hj, if_false]
    · exact hgz
    · rw [omega_addSmul_right, hss p j, omega_isAlt (T.stab p), mul_zero, add_zero,
        if_neg (fun h : p = j => hj h.symm)]
    · rw [omega_addSmul_left, hgz, mul_one, omega_zVec]
      exact hself _
    · rw [omega_addSmul, hss p j, mul_zero, add_zero,
        show omega (T.destab i) (T.stab p) = 0 by simpa [hi] using hds i p,
        mul_zero, add_zero]
      exact hds i j

/-! ### The frame *spans*, and `ω` reads off the coordinates (Yoder-2012 Lemma 5)

`frame_linearIndependent` is only half of "the `2n` rows are a symplectic basis":
it says the generators are independent, not that every Pauli is a product of
them. The other half is what `compute_decomposition`
(`crates/ppvm-tableau-2/src/data.rs`, Yoder-2012 Lemma 5) actually runs on. That
routine reads `destab_anticomm_bits` / `stab_anticomm_bits` straight off the
anticommutation tests, multiplies the corresponding generators into a running
word, and returns `p_word.phase` **without ever checking that the residual word
has collapsed to the identity**. If the anticommutation bits were not `v`'s exact
coordinates in the frame basis, the residual would be a non-identity Pauli and
the returned phase would be meaningless, with no assertion to catch it.

`frame_coordinate_expansion` is exactly that missing guarantee, in the
phase-stripped space: the `ω`-pairings against the *dual* generators are the
coordinates, so the residual really is the identity. Everything downstream rests
on it — the branch relabel `idx ^ stab_anticomm_bits` (every `t`/rotation), the
per-coefficient phase from `destab_anticomm_bits`, `get_deterministic_outcome`'s
claim that `∏_{i : (dᵢ)_x[q]=1} sᵢ = ±Z_q`, and `expectation` /
`compute_decomposition_word`.

The dimension count the old prose deferred is done here by *counting*: the
coordinate map `frameCombine` from `𝔽₂^{2n}` is injective (that is
`frame_linearIndependent`), and its domain and codomain are finite sets of the
same size `4ⁿ`, so it is surjective (`frame_surjective`). No `finrank`/`Basis`
machinery is needed. -/

/-- The coordinate map of a frame: `(a, b) ↦ Σᵢ aᵢ·dᵢ + Σᵢ bᵢ·sᵢ`. Its
surjectivity is "the `2n` generators span". -/
def frameCombine (T : Frame n) (ab : (Fin n → ZMod 2) × (Fin n → ZMod 2)) : Sp n :=
  (∑ i, ab.1 i • T.destab i) + ∑ i, ab.2 i • T.stab i

/-- Pairing the coordinate map against `ω` reads a coordinate back out. -/
private theorem omega_frameCombine (T : Frame n) (ab) (u : Sp n) :
    omega (frameCombine T ab) u
      = (∑ i, ab.1 i * omega (T.destab i) u) + ∑ i, ab.2 i * omega (T.stab i) u := by
  simp only [frameCombine, map_add, LinearMap.add_apply, map_sum, LinearMap.sum_apply,
    map_smul, LinearMap.smul_apply, smul_eq_mul]

/-- **The coordinate map is injective** — `frame_linearIndependent` restated for
the map. -/
theorem frameCombine_injective (T : Frame n) (hT : IsSymplecticFrame T) :
    Function.Injective (frameCombine T) := by
  rintro ⟨a₁, b₁⟩ ⟨a₂, b₂⟩ h
  have h0 : (∑ i, (a₁ i - a₂ i) • T.destab i) + ∑ i, (b₁ i - b₂ i) • T.stab i = 0 := by
    simp only [sub_smul, Finset.sum_sub_distrib]
    simp only [frameCombine] at h
    rw [show ((∑ i, a₁ i • T.destab i) - ∑ i, a₂ i • T.destab i)
          + ((∑ i, b₁ i • T.stab i) - ∑ i, b₂ i • T.stab i)
        = ((∑ i, a₁ i • T.destab i) + ∑ i, b₁ i • T.stab i)
          - ((∑ i, a₂ i • T.destab i) + ∑ i, b₂ i • T.stab i) by ring, h, sub_self]
  obtain ⟨ha, hb⟩ := frame_linearIndependent T hT _ _ h0
  exact Prod.ext (funext fun i => sub_eq_zero.mp (ha i)) (funext fun i => sub_eq_zero.mp (hb i))

/-- **The `2n` generators span** — the half of "symplectic basis" that
`frame_linearIndependent` left open. An injective map between finite sets of
equal cardinality (`|𝔽₂^n × 𝔽₂^n| = 4ⁿ = |Sp n|`) is surjective. -/
theorem frame_surjective (T : Frame n) (hT : IsSymplecticFrame T) :
    Function.Surjective (frameCombine T) := by
  have hcard : Fintype.card ((Fin n → ZMod 2) × (Fin n → ZMod 2)) = Fintype.card (Sp n) := by
    simp only [Fintype.card_prod, Fintype.card_fun, Fintype.card_fin, ZMod.card,
      Fintype.card_prod]
    rw [← mul_pow]
  exact ((Fintype.bijective_iff_injective_and_card _).mpr
    ⟨frameCombine_injective T hT, hcard⟩).2

/-- **The coordinate expansion — `compute_decomposition` is correct.**

For a symplectic frame `T` and *any* Pauli `v`,

  `v = Σᵢ ω(v, sᵢ) · dᵢ + Σᵢ ω(v, dᵢ) · sᵢ`,

i.e. the `2n` generators span *and* a vector's coordinates in that basis are
exactly its `ω`-pairings against the dual generators. `ω(v, sᵢ)` is the crate's
"does `v` anticommute with stabilizer `i`" test, accumulated into
`stab_anticomm_bits`; `ω(v, dᵢ)` is the destabilizer test, accumulated into
`destab_anticomm_bits`. So the two bitmasks `compute_decomposition` builds are
literally the two coordinate vectors above, and multiplying those generators into
`p_word` therefore cancels `v` down to the identity — which is why returning
`p_word.phase` is meaningful and why the unchecked residual can never be a
non-identity Pauli. -/
theorem frame_coordinate_expansion (T : Frame n) (hT : IsSymplecticFrame T) (v : Sp n) :
    (∑ i, omega v (T.stab i) • T.destab i) + (∑ i, omega v (T.destab i) • T.stab i) = v := by
  obtain ⟨hdd, hss, hds⟩ := hT
  have hsd : ∀ i j, omega (T.stab i) (T.destab j) = if j = i then 1 else 0 := by
    intro i j; rw [omega_comm]; exact hds j i
  obtain ⟨ab, hab⟩ := frame_surjective T ⟨hdd, hss, hds⟩ v
  have hs : ∀ j, omega v (T.stab j) = ab.1 j := by
    intro j
    rw [← hab, omega_frameCombine]
    simp only [hds, hss, mul_ite, mul_one, mul_zero, Finset.sum_const_zero, add_zero,
      Finset.sum_ite_eq', Finset.mem_univ, if_true]
  have hd : ∀ j, omega v (T.destab j) = ab.2 j := by
    intro j
    rw [← hab, omega_frameCombine]
    simp only [hdd, hsd, mul_ite, mul_one, mul_zero, Finset.sum_const_zero, zero_add,
      Finset.sum_ite_eq, Finset.mem_univ, if_true]
  simp only [hs, hd]
  exact hab

/-- **`ω` read in frame coordinates.** Pairing two Paulis is the crossed sum of
their `compute_decomposition` masks:

  `ω(v, u) = Σᵢ ω(v, sᵢ)·ω(u, dᵢ) + Σᵢ ω(v, dᵢ)·ω(u, sᵢ)`,

i.e. `⟨stab_anticomm_v, destab_anticomm_u⟩ + ⟨destab_anticomm_v,
stab_anticomm_u⟩`. Immediate from `frame_coordinate_expansion` and bilinearity,
and it is what turns "the two Paulis commute" into the `𝔽₂` hypothesis that
`Tableau/BranchPhase.lean`'s `rot2_order_irrelevant` needs: the `b`-before-`a`
ordering of `rotate_2`'s two single-site relabels shifts the accumulated `ℤ/4`
phase by exactly this quantity. -/
theorem omega_eq_frame_coords (T : Frame n) (hT : IsSymplecticFrame T) (v u : Sp n) :
    omega v u = (∑ i, omega v (T.stab i) * omega u (T.destab i))
      + ∑ i, omega v (T.destab i) * omega u (T.stab i) := by
  conv_lhs => rw [← frame_coordinate_expansion T hT v]
  simp only [map_add, LinearMap.add_apply, map_sum, LinearMap.sum_apply, map_smul,
    LinearMap.smul_apply, smul_eq_mul]
  congr 1 <;>
    exact Finset.sum_congr rfl fun i _ => by congr 1; exact omega_comm _ _

/-- **Paulis with disjoint support commute.** The two single-site Paulis
`rotate_2` conjugates live on *distinct* qubits, so `ω = 0` and — by
`omega_eq_frame_coords` — the order in which their frame-conjugated relabels are
applied cannot change the accumulated phase. -/
theorem omega_disjoint_support (v u : Sp n) (h : ∀ i, v i = 0 ∨ u i = 0) :
    omega v u = 0 := by
  rw [omega_apply]
  refine Finset.sum_eq_zero fun i _ => ?_
  rcases h i with hi | hi <;> rw [hi] <;> simp

end PPVM.Tableau
