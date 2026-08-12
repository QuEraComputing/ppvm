/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Tableau.Frame
import PPVM.Tableau.Projection

/-!
# The *concrete* branch phase: discharging `SelfInverse`, and the two-site product

`Tableau/Projection.lean` states **every** case-a and case-b theorem under the
hypothesis `SelfInverse s φ` (`shiftOp_involutive`, `shiftOp_selfAdjoint`,
`proj_idem`, `probOne_eq`, `selfInverse_zero_phase_even`, `proj_zero_apply`,
`proj_zero_eq_caseB_retain`) and models `φ` *abstractly*. The crate does not have
an abstract `φ`: it computes one particular function
(`crates/ppvm-tableau-2/src/data.rs`, `compute_phase_with_mask_static`, folded
with `phase_decomp` from `compute_decomposition`),

  `φ j = phase_decomp + 2·⟨destab_anticomm, j⟩ + 2·popcount(j ∧ stab_anticomm ∧ oddPhaseMask)`.

Until that function is shown to be `SelfInverse`, the whole projection oracle is
vacuously true of the implementation. This file closes that link, and then reuses
the same closed forms for the two-site rotation.

## The model

The amplitude basis is `|j⟩ = ∏_l d_l^{j_l} |ψ₀⟩` — a bitstring `j` selects a
subset of destabilizers to apply to the stabilizer state `|ψ₀⟩`. Three facts
about a *valid* frame turn the generators into closed-form operators on
`Amp n = Bitstring n → ℂ`:

* `s_k |ψ₀⟩ = |ψ₀⟩` and `ω(s_k, d_l) = δ_{kl}`, so `s_k |j⟩ = (−1)^{j_k} |j⟩` and
  hence `S_G |j⟩ = (−1)^{⟨G, j⟩} |j⟩` — `stabAction`;
* a row is `i^phase ·` (a tensor of Hermitian Paulis), so `d_l² = (−1)^{phase_l}`;
  destabilizers commute, so `D_L |j⟩ = (−1)^{popcount(j ∧ L ∧ m)} |j ⊕ L⟩` with
  `m` the odd-phase mask — `destabAction`, and this is *exactly* where the
  crate's `odd_phase_destabilizer_mask` term comes from;
* `compute_decomposition` visits **all stabilizers first, then all
  destabilizers**, multiplying each generator's inverse into the running word, so
  what it returns is the `φ` of `P = i^φ · D_L · S_G` — destabilizers on the
  *left*.

`frameOp` is that composite, and `frameOp_eq_shiftOp` proves it **is** the
crate's formula: the phase function is *derived* here, not assumed.

## What is proved

* `frameOp_eq_shiftOp` — `i^{pd} · D_L · S_G` acting on amplitudes is
  `shiftOp L (branchPhase pd (crateWeight G L m))`, i.e. literally
  `compute_phase_with_mask_static` plus `phase_decomp`. This also pins the visit
  order: `stab_destab_commute_sign` shows the *opposite* order shifts `φ` by
  `2·⟨G, L⟩`, so "all stabilizers first" is a genuine convention and the crate's
  use of the **original** index `idx` (not `idx ⊕ L`) in the `⟨destab, ·⟩` term
  is the one that matches it.
* `frameOp_sq` / `frameOp_involutive_iff` — `M² = i^{2·pd + 2(⟨G,L⟩ + ⟨L,m⟩)}·I`,
  so `M² = I` **iff** the frame identity
  `pd + ⟨G, L⟩ + popcount(L ∧ m) ≡ 0 (mod 2)` holds.
* `selfInverse_branchPhase` / `selfInverse_branchPhase_iff` — that same identity
  is *equivalent* to `SelfInverse L (branchPhase …)`. Since the conjugated Pauli
  is an involution (`Z_q² = I`, and conjugation preserves it), the hypothesis is
  discharged and every theorem in `Projection.lean` now applies to the function
  the crate actually computes: `shiftOp_involutive_crate`, `proj_idem_crate`,
  `probOne_eq_crate`.
* `frameInvolution_zero_iff` — in case b (`L = 0`) the frame identity collapses
  to `pd ∈ {0, 2}`. That is the crate's
  `debug_assert!(phase_decomp == 0 || phase_decomp == 2)` ("Measurement result
  cannot be imaginary!", `measure.rs`) — a runtime guard replaced by a theorem.
* `destabAction_sq` — `d_L² = (−1)^{popcount(L ∧ m)}`, i.e. a generator is its own
  inverse *up to its phase*. That is the algebraic content of
  `compute_decomposition`'s `p_word.add_phase(8 − 2·g.phase)` step
  (`add_phase_eight_sub` records `8 − 2·ph = 2·ph` in `ℤ/4`).
* `shiftOp_comp` — composing two of these operators gives another one, with shift
  `L₁ ⊕ L₂`, weight `w₁ + w₂` and phase `pd₁ + pd₂ + 2⟨w₁, L₂⟩`. This is what
  `rotate_2` does (`gates.rs`: `compute_coefficients_after_pauli_apply` at `b`,
  then at `a`): the sequential relabel **is** the frame-conjugated two-site Pauli,
  phase included.
* `shiftOp_comp_order_iff` / `rot2_order_irrelevant` — the two application orders
  differ by exactly `2·(⟨G₁, L₂⟩ + ⟨G₂, L₁⟩)`, which `omega_eq_frame_coords`
  (`Frame.lean`) identifies as `ω(P₁, P₂)` read in frame coordinates. Two Paulis
  on **disjoint sites commute** (`omega_disjoint_support`), so for `rotate_2` the
  `b`-before-`a` order is provably *irrelevant to the `ℤ/4` phase*. (The crate
  keeps the old order anyway: the prime directive pins the float summation order,
  which is a different concern.)

## Conformance

`crates/ppvm-conformance-2/tests/tableau_lean.rs`
(`decomposition_satisfies_the_lean_frame_identity`,
`rot2_application_order_is_phase_neutral`) checks the frame identity and the
`rot2_order_irrelevant` hypothesis on every decomposition a random Clifford+`T`
sweep produces (thousands of them), so the discharge is tied to the code, not
only to the model.

That sweep also records that `odd_phase_mask` is **always zero** in a valid
frame: rows are `i^phase ·` (Hermitian Pauli tensor), and Hermiticity forces
`phase ∈ {0, 2}`, so `dot L m = 0` and the mask term of `crateWeight` cannot
contribute. The crate keeps it because the old crate computes it (the
behaviour-preservation directive), and none of the theorems below need it to
vanish — they are stated for arbitrary `m`.

## Residual

The one input this file does not derive is the model itself — that a valid frame
satisfies the three bullet points above. That needs a Hilbert-space model of the
frame, the same gap `Projection.lean`'s scope note records. What it *does* remove
is the ∀-quantified `SelfInverse` unknown: everything now hangs on the single
scalar fact `P² = I`.
-/

namespace PPVM.Tableau.BranchPhase

open PPVM.GenTableau PPVM.Tableau.Projection PPVM.Symplectic

variable {n : ℕ}

/-! ### `𝔽₂` bit arithmetic: popcount parities and masks -/

/-- `⟨w, j⟩` — the parity of `popcount(w ∧ j)`, the crate's `symplectic_inner`
reduced mod 2. -/
def dot (w j : Bitstring n) : ZMod 2 := ∑ i, w i * j i

/-- Bitwise `AND` of two bitstrings (pointwise `𝔽₂` product). -/
def band (a b : Bitstring n) : Bitstring n := fun i => a i * b i

theorem dot_comm (w j : Bitstring n) : dot w j = dot j w :=
  Finset.sum_congr rfl fun _ _ => mul_comm _ _

theorem dot_add_right (w j k : Bitstring n) : dot w (j + k) = dot w j + dot w k := by
  simp only [dot, Pi.add_apply, mul_add]
  exact Finset.sum_add_distrib

theorem dot_add_left (w v j : Bitstring n) : dot (w + v) j = dot w j + dot v j := by
  rw [dot_comm, dot_add_right, dot_comm w, dot_comm v]

@[simp] theorem dot_zero_right (w : Bitstring n) : dot w 0 = 0 := by
  simp [dot]

@[simp] theorem dot_zero_left (w : Bitstring n) : dot 0 w = 0 := by
  simp [dot]

/-- `popcount(L ∧ m ∧ L) = popcount(L ∧ m)` — idempotence of the mask. -/
theorem dot_band_self (L m : Bitstring n) : dot (band L m) L = dot L m := by
  refine Finset.sum_congr rfl fun i _ => ?_
  have h : ∀ x : ZMod 2, x * x = x := by decide
  calc L i * m i * L i = L i * L i * m i := by ring
    _ = L i * m i := by rw [h]

/-- `popcount(L₁ ∧ m ∧ L₂)` is symmetric in `L₁`, `L₂` — the mask term of the
composite phase carries **no** ordering information. -/
theorem dot_band_symm (L₁ L₂ m : Bitstring n) :
    dot (band L₁ m) L₂ = dot (band L₂ m) L₁ :=
  Finset.sum_congr rfl fun i _ => by simp only [band]; ring

/-! ### `ℤ/2 → ℤ/4` doubling -/

/-- `x ↦ 2x`, the `ℤ/2 → ℤ/4` map every sign in the crate's phase is expressed
through. -/
def two2 (x : ZMod 2) : ZMod 4 := if x = 1 then 2 else 0

/-- Reduction `ℤ/4 → ℤ/2` (a phase is "imaginary" iff its reduction is `1`). -/
def parity4 (a : ZMod 4) : ZMod 2 := (a.val : ZMod 2)

@[simp] theorem two2_zero : two2 (0 : ZMod 2) = 0 := rfl

theorem two2_add (a b : ZMod 2) : two2 (a + b) = two2 a + two2 b := by decide +revert

theorem two2_self_add (a : ZMod 2) : two2 a + two2 a = 0 := by decide +revert

theorem two2_eq_zero_iff (a : ZMod 2) : two2 a = 0 ↔ a = 0 := by decide +revert

/-- `2·a` in `ℤ/4` only sees `a`'s parity. -/
theorem two2_parity4 (a : ZMod 4) : two2 (parity4 a) = 2 * a := by decide +revert

theorem parity4_eq_zero_iff (a : ZMod 4) : parity4 a = 0 ↔ (a = 0 ∨ a = 2) := by
  decide +revert

/-- `i^k = 1` only for `k = 0` — the injectivity used to read a scalar back out
of an operator identity. -/
theorem iPow_eq_one_iff (a : ZMod 4) : iPow a = 1 ↔ a = 0 := by
  constructor
  · intro h
    rcases zmod4_cases a with ha | ha | ha | ha <;> subst ha
    · rfl
    · rw [iPow_one] at h; exact absurd h (by simp [Complex.ext_iff])
    · rw [iPow_two] at h; exact absurd h (by norm_num)
    · rw [iPow_three] at h; exact absurd h (by simp [Complex.ext_iff])
  · rintro rfl; exact iPow_zero

/-! ### The crate's phase function -/

/-- The crate's per-coefficient phase, in closed form:
`φ j = phase_decomp + 2·⟨w, j⟩`. Both of `compute_phase_with_mask_static`'s
popcount terms are `𝔽₂`-linear in `j`, so they fuse into a single weight vector
`w` (`crateWeight`). -/
def branchPhase (pd : ZMod 4) (w : Bitstring n) : Bitstring n → ZMod 4 :=
  fun j => pd + two2 (dot w j)

/-- The weight the crate builds: `destab_anticomm + (stab_anticomm ∧
odd_phase_mask)`. -/
def crateWeight (G L m : Bitstring n) : Bitstring n := G + band L m

/-- **`branchPhase` on `crateWeight` is `compute_phase_with_mask_static` verbatim**:
`phase_decomp + 2·⟨destab_anticomm, j⟩ + 2·popcount(j ∧ stab_anticomm ∧ mask)`. -/
theorem branchPhase_crateWeight (pd : ZMod 4) (G L m j : Bitstring n) :
    branchPhase pd (crateWeight G L m) j
      = pd + two2 (dot G j) + two2 (dot (band L m) j) := by
  simp only [branchPhase, crateWeight, dot_add_left, two2_add, add_assoc]

/-- The phase difference across the branch pair is **independent of `j`**: the
`j`-linear part cancels, leaving `2·pd + 2⟨w, s⟩`. -/
theorem branchPhase_pair (pd : ZMod 4) (w s j : Bitstring n) :
    branchPhase pd w j + branchPhase pd w (j + s) = 2 * pd + two2 (dot w s) := by
  simp only [branchPhase, dot_add_right, two2_add]
  rw [show pd + two2 (dot w j) + (pd + (two2 (dot w j) + two2 (dot w s)))
      = 2 * pd + (two2 (dot w j) + two2 (dot w j)) + two2 (dot w s) by ring,
    two2_self_add, add_zero]

/-! ### The frame identity, and `SelfInverse` for the concrete phase -/

/-- **The frame identity.** `phase_decomp + ⟨destab_anticomm, stab_anticomm⟩ +
popcount(stab_anticomm ∧ oddPhaseMask) ≡ 0 (mod 2)` — equivalently `M² = I`
(`frameOp_involutive_iff`), i.e. the conjugated Pauli is an involution. -/
def FrameInvolution (pd : ZMod 4) (G L m : Bitstring n) : Prop :=
  parity4 pd + dot G L + dot L m = 0

/-- The single scalar that `SelfInverse` reduces to, in crate terms. -/
theorem branchPhase_pair_crateWeight (pd : ZMod 4) (G L m j : Bitstring n) :
    branchPhase pd (crateWeight G L m) j + branchPhase pd (crateWeight G L m) (j + L)
      = two2 (parity4 pd + dot G L + dot L m) := by
  rw [branchPhase_pair, crateWeight, dot_add_left, dot_band_self]
  simp only [two2_add, two2_parity4]
  ring

/-- **`SelfInverse` holds exactly on the frame identity.** The `∀ j` hypothesis
carried by every theorem in `Projection.lean` is equivalent to one `ℤ/2`
equation about the three masks `compute_decomposition` returns. -/
theorem selfInverse_branchPhase_iff (pd : ZMod 4) (G L m : Bitstring n) :
    SelfInverse L (branchPhase pd (crateWeight G L m)) ↔ FrameInvolution pd G L m := by
  constructor
  · intro h
    have h0 := h 0
    rw [branchPhase_pair_crateWeight] at h0
    exact (two2_eq_zero_iff _).mp h0
  · intro h j
    rw [branchPhase_pair_crateWeight, h]
    rfl

/-- **The crate's branch phase is `SelfInverse`** (given the frame identity). -/
theorem selfInverse_branchPhase {pd : ZMod 4} {G L m : Bitstring n}
    (h : FrameInvolution pd G L m) :
    SelfInverse L (branchPhase pd (crateWeight G L m)) :=
  (selfInverse_branchPhase_iff pd G L m).mpr h

/-! ### The generator actions, and that the crate's formula *is* their composite -/

/-- `S_G |j⟩ = (−1)^{⟨G, j⟩} |j⟩` on amplitudes — the diagonal action of the
stabilizer product selected by `destab_anticomm_bits`. -/
noncomputable def stabAction (G : Bitstring n) (c : Amp n) : Amp n :=
  fun j => iPow (two2 (dot G j)) * c j

/-- `D_L |j⟩ = (−1)^{popcount(j ∧ L ∧ m)} |j ⊕ L⟩` on amplitudes — the destabilizer
product selected by `stab_anticomm_bits`, with the sign coming from `d_l² =
(−1)^{phase_l}` on the doubly-selected generators (`m` = `odd_phase_mask`). -/
noncomputable def destabAction (L m : Bitstring n) (c : Amp n) : Amp n :=
  fun j => iPow (two2 (dot (band L m) (j + L))) * c (j + L)

/-- The frame-conjugated Pauli `P = i^{pd} · D_L · S_G` acting on amplitudes. -/
noncomputable def frameOp (pd : ZMod 4) (G L m : Bitstring n) (c : Amp n) : Amp n :=
  fun j => iPow pd * destabAction L m (stabAction G c) j

/-- **`d_L` is its own inverse up to its phase**: `D_L² = (−1)^{popcount(L ∧ m)}`.
This is what licenses `compute_decomposition`'s
`p_word.mul_assign(g); p_word.add_phase(8 − 2·g.phase)` — multiply the generator
in and divide its phase squared out, instead of inverting it. -/
theorem destabAction_sq (L m : Bitstring n) (c : Amp n) (j : Bitstring n) :
    destabAction L m (destabAction L m c) j = iPow (two2 (dot L m)) * c j := by
  have hjj : j + L + L = j := add_shift_shift L j
  have key : dot (band L m) (j + L) + dot (band L m) j = dot L m := by
    rw [← dot_add_right, show j + L + j = L by rw [add_right_comm, add_self, zero_add]]
    exact dot_band_self L m
  simp only [destabAction, hjj]
  rw [← mul_assoc, ← iPow_add, ← two2_add, key]

/-- `8 − 2·ph = 2·ph` in `ℤ/4`: the crate's `add_phase(8 - 2 * g.phase)` really is
"divide out the generator's phase squared". -/
theorem add_phase_eight_sub (ph : ZMod 4) : 8 - 2 * ph = 2 * ph := by decide +revert

/-- **The stabilizer/destabilizer visit order is a real convention**: swapping
`S_G` and `D_L` shifts the phase by `2·⟨G, L⟩`. `compute_decomposition`'s two
loops fix the order `D_L · S_G` (destabilizers outermost), which is why the
`⟨destab_anticomm, ·⟩` term is evaluated at the **original** index `idx` rather
than at the relabelled `idx ⊕ stab_anticomm`. -/
theorem stab_destab_commute_sign (G L m : Bitstring n) (c : Amp n) (j : Bitstring n) :
    destabAction L m (stabAction G c) j
      = iPow (two2 (dot G L)) * stabAction G (destabAction L m c) j := by
  simp only [destabAction, stabAction]
  rw [dot_add_right G j L, two2_add, iPow_add]
  ring

/-- **The crate's phase formula is the composite of the generator actions.**
`i^{pd} · D_L · S_G` on amplitudes is exactly
`shiftOp stab_anticomm (branchPhase phase_decomp (crateWeight …))`, i.e.
`compute_phase_with_mask_static` folded with `phase_decomp`. The formula is
*derived* here rather than assumed. -/
theorem frameOp_eq_shiftOp (pd : ZMod 4) (G L m : Bitstring n) (c : Amp n) :
    frameOp pd G L m c = shiftOp L (branchPhase pd (crateWeight G L m)) c := by
  funext j
  simp only [frameOp, destabAction, stabAction, shiftOp, branchPhase, crateWeight,
    dot_add_left, two2_add]
  rw [iPow_add, iPow_add]
  ring

/-! ### `M² = I` ⇔ the frame identity -/

/-- Squaring a phased XOR shift collapses to a scalar. -/
theorem shiftOp_sq (s : Bitstring n) (φ : Bitstring n → ZMod 4) (c : Amp n)
    (j : Bitstring n) :
    shiftOp s φ (shiftOp s φ c) j = iPow (φ j + φ (j + s)) * c j := by
  simp only [shiftOp, add_shift_shift]
  rw [← mul_assoc, ← iPow_add, add_comm (φ (j + s)) (φ j)]

/-- **`M² = i^{2·pd + 2(⟨G,L⟩ + ⟨L,m⟩)} · I`.** -/
theorem frameOp_sq (pd : ZMod 4) (G L m : Bitstring n) (c : Amp n) (j : Bitstring n) :
    frameOp pd G L m (frameOp pd G L m c) j
      = iPow (two2 (parity4 pd + dot G L + dot L m)) * c j := by
  rw [frameOp_eq_shiftOp, frameOp_eq_shiftOp, shiftOp_sq, branchPhase_pair_crateWeight]

/-- **`M² = I` iff the frame identity holds.** The conjugated Pauli is an
involution (`Z_q² = I` and conjugation preserves that), so the right-hand side is
a fact about the frame, not an assumption about `φ`. -/
theorem frameOp_involutive_iff (pd : ZMod 4) (G L m : Bitstring n) :
    (∀ (c : Amp n) (j : Bitstring n), frameOp pd G L m (frameOp pd G L m c) j = c j)
      ↔ FrameInvolution pd G L m := by
  constructor
  · intro h
    have h0 := h (fun _ => 1) 0
    rw [frameOp_sq, mul_one] at h0
    exact (two2_eq_zero_iff _).mp ((iPow_eq_one_iff _).mp h0)
  · intro h c j
    rw [frameOp_sq, h]
    simp

/-! ### The `Projection.lean` oracle, now applied to the crate's function -/

/-- **`M² = I` for the phase the crate computes** (`shiftOp_involutive`
discharged). -/
theorem shiftOp_involutive_crate {pd : ZMod 4} {G L m : Bitstring n}
    (h : FrameInvolution pd G L m) :
    Function.Involutive (shiftOp L (branchPhase pd (crateWeight G L m))) :=
  shiftOp_involutive (selfInverse_branchPhase h)

/-- **`P_b² = P_b` for the phase the crate computes** (`proj_idem` discharged). -/
theorem proj_idem_crate {pd : ZMod 4} {G L m : Bitstring n}
    (h : FrameInvolution pd G L m) (b : Bool) (c : Amp n) :
    proj b L (branchPhase pd (crateWeight G L m))
        (proj b L (branchPhase pd (crateWeight G L m)) c)
      = proj b L (branchPhase pd (crateWeight G L m)) c :=
  proj_idem (selfInverse_branchPhase h) b c

/-- **The Born rule for the phase the crate computes** (`probOne_eq` discharged):
`prob_1 = 0.5 − 0.5·z_overlap_re` is `⟨c, P₁ c⟩` for the *concrete* `φ`. -/
theorem probOne_eq_crate {pd : ZMod 4} {G L m : Bitstring n}
    (h : FrameInvolution pd G L m) (c : Amp n)
    (hnorm : ∑ j, (starRingEnd ℂ) (c j) * c j = 1) :
    (∑ j, (starRingEnd ℂ) (c j) *
        proj true L (branchPhase pd (crateWeight G L m)) c j).re
      = 1 / 2 - 1 / 2 * overlapRe L (branchPhase pd (crateWeight G L m)) c :=
  probOne_eq (selfInverse_branchPhase h) c hnorm

/-- **Case b: "measurement result cannot be imaginary" is a theorem.** With
`stab_anticomm = 0` the frame identity collapses to `phase_decomp ∈ {0, 2}` —
exactly the crate's `debug_assert!(phase_decomp == 0 || phase_decomp == 2)`
(`measure.rs`), and the premise `selfInverse_zero_phase_even` runs on. -/
theorem frameInvolution_zero_iff (pd : ZMod 4) (G m : Bitstring n) :
    FrameInvolution pd G 0 m ↔ (pd = 0 ∨ pd = 2) := by
  simp only [FrameInvolution, dot_zero_right, dot_zero_left, add_zero]
  exact parity4_eq_zero_iff pd

/-! ### `rotate_2`: the sequential relabel is the two-site product

`gates.rs`'s `rotate_2` applies `compute_coefficients_after_pauli_apply` at `b`
and then at `a` — two independent single-site relabels on the amplitude vector,
never a two-site decomposition. `shiftOp_comp` says the composite is again an
operator of exactly the same shape, so it *is* the frame-conjugated
`P_a ⊗ P_b`, phase included; `rot2_order_irrelevant` says the `b`-before-`a`
order does not affect that phase. -/

/-- The composite phase, at the level of `branchPhase` alone. -/
theorem branchPhase_comp (pd₁ pd₂ : ZMod 4) (w₁ w₂ s₂ k : Bitstring n) :
    branchPhase pd₁ w₁ (k + s₂) + branchPhase pd₂ w₂ k
      = branchPhase (pd₁ + pd₂ + two2 (dot w₁ s₂)) (w₁ + w₂) k := by
  simp only [branchPhase, dot_add_right, dot_add_left, two2_add]
  ring

/-- **Composition law.** Two phased XOR shifts compose to a phased XOR shift with
shift `s₁ ⊕ s₂`, weight `w₁ + w₂` and phase `pd₁ + pd₂ + 2⟨w₁, s₂⟩`. -/
theorem shiftOp_comp (pd₁ pd₂ : ZMod 4) (w₁ w₂ s₁ s₂ : Bitstring n) (c : Amp n)
    (j : Bitstring n) :
    shiftOp s₁ (branchPhase pd₁ w₁) (shiftOp s₂ (branchPhase pd₂ w₂) c) j
      = shiftOp (s₁ + s₂)
          (branchPhase (pd₁ + pd₂ + two2 (dot w₁ s₂)) (w₁ + w₂)) c j := by
  have hk : j + s₁ + s₂ = j + (s₁ + s₂) := add_assoc _ _ _
  have hj : j + (s₁ + s₂) + s₂ = j + s₁ := by
    rw [add_assoc, add_assoc, add_self, add_zero]
  simp only [shiftOp]
  rw [hk, ← hj, ← mul_assoc, ← iPow_add, branchPhase_comp]

/-- **The order dependence is exactly `⟨G₁, L₂⟩ + ⟨G₂, L₁⟩`.** The mask term is
symmetric (`dot_band_symm`), so *all* of it sits in the destabilizer/stabilizer
cross pairing — which is `ω(P₁, P₂)`. -/
theorem dot_crateWeight_order (G₁ L₁ G₂ L₂ m : Bitstring n) :
    dot (crateWeight G₁ L₁ m) L₂ + dot (crateWeight G₂ L₂ m) L₁
      = dot G₁ L₂ + dot G₂ L₁ := by
  have hcc : dot (band L₂ m) L₁ + dot (band L₂ m) L₁ = 0 := by
    have h : ∀ x : ZMod 2, x + x = 0 := by decide
    exact h _
  simp only [crateWeight, dot_add_left, dot_band_symm L₁ L₂ m]
  rw [show dot G₁ L₂ + dot (band L₂ m) L₁ + (dot G₂ L₁ + dot (band L₂ m) L₁)
      = dot G₁ L₂ + dot G₂ L₁ + (dot (band L₂ m) L₁ + dot (band L₂ m) L₁) from by ring,
    hcc, add_zero]

theorem two2_inj_iff (a b : ZMod 2) : two2 a = two2 b ↔ a + b = 0 := by decide +revert

/-- **The two application orders carry the same `ℤ/4` phase iff the two Paulis
commute** in frame coordinates. -/
theorem shiftOp_comp_order_iff (G₁ L₁ G₂ L₂ m : Bitstring n) :
    two2 (dot (crateWeight G₁ L₁ m) L₂) = two2 (dot (crateWeight G₂ L₂ m) L₁)
      ↔ dot G₁ L₂ + dot G₂ L₁ = 0 := by
  rw [two2_inj_iff, dot_crateWeight_order]

/-- **The `b`-before-`a` order in `rotate_2` is irrelevant to the `ℤ/4` phase**
whenever the two Paulis commute in frame coordinates. -/
theorem rot2_order_irrelevant (pd₁ pd₂ : ZMod 4) (G₁ L₁ G₂ L₂ m : Bitstring n)
    (hcomm : dot G₁ L₂ + dot G₂ L₁ = 0) (c : Amp n) (j : Bitstring n) :
    shiftOp L₁ (branchPhase pd₁ (crateWeight G₁ L₁ m))
        (shiftOp L₂ (branchPhase pd₂ (crateWeight G₂ L₂ m)) c) j
      = shiftOp L₂ (branchPhase pd₂ (crateWeight G₂ L₂ m))
          (shiftOp L₁ (branchPhase pd₁ (crateWeight G₁ L₁ m)) c) j := by
  rw [shiftOp_comp, shiftOp_comp,
    (shiftOp_comp_order_iff G₁ L₁ G₂ L₂ m).mpr hcomm,
    add_comm L₂ L₁, add_comm pd₂ pd₁,
    add_comm (crateWeight G₂ L₂ m) (crateWeight G₁ L₁ m)]

/-! ### Tying the order hypothesis to the symplectic form

`omega_eq_frame_coords` (`Frame.lean`) says `ω(v, u)` read in frame coordinates is
`⟨S_v, D_u⟩ + ⟨D_v, S_u⟩`, where `S` is the `stab_anticomm` mask and `D` the
`destab_anticomm` mask `compute_decomposition` returns. So the hypothesis of
`rot2_order_irrelevant` is literally "the two Paulis commute", and single-site
Paulis on distinct qubits do. -/

/-- The crate's `stab_anticomm_bits` as a bitstring: bit `i` is `ω(v, sᵢ)`. -/
def stabMask (T : Frame n) (v : Sp n) : Bitstring n := fun i => omega v (T.stab i)

/-- The crate's `destab_anticomm_bits` as a bitstring: bit `i` is `ω(v, dᵢ)`. -/
def destabMask (T : Frame n) (v : Sp n) : Bitstring n := fun i => omega v (T.destab i)

/-- **`ω` in the crate's own masks**: `ω(v, u) = ⟨stab_v, destab_u⟩ +
⟨destab_v, stab_u⟩`. -/
theorem omega_eq_masks (T : Frame n) (hT : IsSymplecticFrame T) (v u : Sp n) :
    omega v u = dot (stabMask T v) (destabMask T u) + dot (destabMask T v) (stabMask T u) :=
  omega_eq_frame_coords T hT v u

/-- **Commuting Paulis ⇒ the `rotate_2` order does not matter.** -/
theorem rot2_order_irrelevant_of_commuting (T : Frame n) (hT : IsSymplecticFrame T)
    (v u : Sp n) (hcomm : omega v u = 0) (pd₁ pd₂ : ZMod 4) (m : Bitstring n)
    (c : Amp n) (j : Bitstring n) :
    shiftOp (stabMask T v)
        (branchPhase pd₁ (crateWeight (destabMask T v) (stabMask T v) m))
        (shiftOp (stabMask T u)
          (branchPhase pd₂ (crateWeight (destabMask T u) (stabMask T u) m)) c) j
      = shiftOp (stabMask T u)
          (branchPhase pd₂ (crateWeight (destabMask T u) (stabMask T u) m))
          (shiftOp (stabMask T v)
            (branchPhase pd₁ (crateWeight (destabMask T v) (stabMask T v) m)) c) j := by
  refine rot2_order_irrelevant _ _ _ _ _ _ _ ?_ c j
  have h := omega_eq_masks T hT v u
  rw [hcomm, dot_comm (stabMask T v) (destabMask T u)] at h
  rw [add_comm]
  exact h.symm

end PPVM.Tableau.BranchPhase
