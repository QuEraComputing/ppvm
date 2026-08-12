/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Data.Complex.BigOperators
import Mathlib.Tactic.LinearCombination
import PPVM.Instantiations.Bitstring

/-!
# Case-a measurement: the `ℤ/4` overlap sign table and the Born projector

`crates/ppvm-tableau-2/src/measure.rs` (`measure_with_scratch`, case a) and
`project_case_a` carry the most load-bearing arithmetic in the generalized
tableau. `Tableau/Frame.lean` proves only *which* branch is taken
(`measurement_dichotomy`, `measure_deterministic_iff_xfree`) and
`Instantiations/Bitstring.lean` proves only that the XOR relabel is a bijection —
neither says anything about what the case-a branch *computes*. The `ℤ/4`
bookkeeping is exactly where a sign error yields a plausible-but-wrong
probability, and the Rust side only checks it against itself.

## The model

Conjugating `Z_q` through the frame yields, in the amplitude basis, an operator
of the shape

  `M |k⟩ = i^{φ k} |k ⊕ s⟩`,   `s = stab_anticomm_bits ≠ 0`,

with `φ : Bitstring n → ℤ/4` the per-index phase the crate computes. On
amplitudes that is `(M c)_j = i^{φ (j ⊕ s)} · c_{j ⊕ s}` (`shiftOp`). `M` is a
unitary conjugate of a Hermitian involution, so `φ k + φ (k ⊕ s) = 0` — the one
hypothesis (`SelfInverse`) everything below runs on.

## What is proved

* `rustTerm_eq` — the merge loop's four-way `ℤ/4` dispatch
  (`0 => +re_w, 1 => +im_w, 2 => -re_w, 3 => -im_w` on
  `re_w = aᵣbᵣ + aᵢbᵢ`, `im_w = aᵣbᵢ − aᵢbᵣ`) computes exactly
  `Re(conj(i^φ · a) · b)`. The convention is *not* `Re(i^φ · conj a · b)`: the odd
  branches carry the conjugated phase, and that is precisely where a sign slip
  would live.
* `shiftOp_involutive`, `shiftOp_selfAdjoint` — `M² = I` and `M† = M`.
* `overlap_eq_inner` — the crate's `z_overlap_re` is `Re ⟨c, M c⟩`.
* `proj_add`, `proj_idem` — `P₀ + P₁ = I` and `P_b² = P_b` for
  `P_b = (I + (−1)^b M)/2`: the case-a outcome really is a projective measurement.
* `probOne_eq` — `prob_1 = 0.5 − 0.5 · z_overlap_re` **is** `⟨c, P₁ c⟩`, the Born
  probability of outcome `1`, for a normalized `c`.
* `projectRaw_eq_two_proj` — the keep-`A`/transform-`B`/merge map of the case-a
  arm is `2 · P_b` on the surviving half, i.e. the projector up to the
  normalization that the subsequent unconditional `normalize()` supplies.

## Scope note

The residual link this file does **not** close: that the `ψ` the projection
actually uses (`alpha + 2·⟨idx, destab_anticomm⟩`, i.e. *without* the
odd-phase-destabilizer term that the *overlap* folds in through
`compute_phase_with_mask_static`) satisfies the `hψ` hypothesis of
`projectRaw_eq_two_proj` against the overlap's `φ`. The two phase functions are
read in different bases — the overlap in the pre-measurement frame, the projected
amplitudes in the post-projection frame installed by
`update_tableau_according_to_outcome` — so relating them needs a Hilbert-space
model of the frame that this development does not have. The old and the `-2`
crate agree here verbatim, so it is a *specification* gap, not a port divergence.
-/

namespace PPVM.Tableau.Projection

open PPVM.GenTableau

variable {n : ℕ}

/-! ### `i^k` for `k : ℤ/4` -/

/-- `i^k` for a `ℤ/4` exponent — the crate's `COMPLEX_PHASE_CONVERSION` table
`[1, i, −1, −i]`. -/
noncomputable def iPow (k : ZMod 4) : ℂ := Complex.I ^ k.val

theorem zmod4_cases : ∀ k : ZMod 4, k = 0 ∨ k = 1 ∨ k = 2 ∨ k = 3 := by decide

@[simp] theorem iPow_zero : iPow 0 = 1 := by
  simp only [iPow, show (0 : ZMod 4).val = 0 from rfl, pow_zero]

@[simp] theorem iPow_one : iPow 1 = Complex.I := by
  simp only [iPow, show (1 : ZMod 4).val = 1 from rfl, pow_one]

@[simp] theorem iPow_two : iPow 2 = -1 := by
  simp only [iPow, show (2 : ZMod 4).val = 2 from rfl, Complex.I_sq]

@[simp] theorem iPow_three : iPow 3 = -Complex.I := by
  simp only [iPow, show (3 : ZMod 4).val = 3 from rfl, Complex.I_pow_three]

/-- `i^{a+b} = i^a · i^b`: the `ℤ/4` phase exponent really is an exponent. -/
theorem iPow_add (a b : ZMod 4) : iPow (a + b) = iPow a * iPow b := by
  simp only [iPow, ← pow_add, ZMod.val_add]
  exact (Complex.I_pow_eq_pow_mod _).symm

/-- Conjugating a phase negates its `ℤ/4` exponent. -/
theorem iPow_conj (a : ZMod 4) : (starRingEnd ℂ) (iPow a) = iPow (-a) := by
  rcases zmod4_cases a with h | h | h | h <;> subst h
  · simp
  · rw [show -(1 : ZMod 4) = 3 by decide]; simp
  · rw [show -(2 : ZMod 4) = 2 by decide]; simp
  · rw [show -(3 : ZMod 4) = 1 by decide]; simp

/-! ### The overlap merge's `ℤ/4` sign table -/

/-- The case-a merge's per-key contribution, verbatim: `re_w = aᵣbᵣ + aᵢbᵢ`,
`im_w = aᵣbᵢ − aᵢbᵣ`, dispatched on the `ℤ/4` phase. -/
noncomputable def rustTerm (φ : ZMod 4) (a b : ℂ) : ℝ :=
  if φ = 0 then a.re * b.re + a.im * b.im
  else if φ = 1 then a.re * b.im - a.im * b.re
  else if φ = 2 then -(a.re * b.re + a.im * b.im)
  else -(a.re * b.im - a.im * b.re)

/-- **The merge's four-way `ℤ/4` dispatch is `Re(conj(i^φ · a) · b)`.** -/
theorem rustTerm_eq (φ : ZMod 4) (a b : ℂ) :
    rustTerm φ a b = ((starRingEnd ℂ) (iPow φ * a) * b).re := by
  have h10 : ¬((1 : ZMod 4) = 0) := by decide
  have h20 : ¬((2 : ZMod 4) = 0) := by decide
  have h21 : ¬((2 : ZMod 4) = 1) := by decide
  have h30 : ¬((3 : ZMod 4) = 0) := by decide
  have h31 : ¬((3 : ZMod 4) = 1) := by decide
  have h32 : ¬((3 : ZMod 4) = 2) := by decide
  rcases zmod4_cases φ with h | h | h | h <;> subst h <;>
    simp only [rustTerm, h10, h20, h21, h30, h31, h32, if_true, if_false, iPow_zero,
      iPow_one, iPow_two, iPow_three, map_mul,
      Complex.mul_re, Complex.mul_im, Complex.conj_re, Complex.conj_im,
      Complex.one_re, Complex.one_im, Complex.I_re, Complex.I_im, Complex.neg_re,
      Complex.neg_im] <;> ring

/-! ### The conjugated `Z_q` as an amplitude operator -/

/-- Amplitude vectors as plain functions over the (finite) index type; the
support-sparse `Finsupp` view is `PPVM.GenTableau.Amplitudes`. -/
abbrev Amp (n : ℕ) := Bitstring n → ℂ

/-- `M |k⟩ = i^{φ k} |k ⊕ s⟩` acting on amplitudes:
`(M c)_j = i^{φ (j ⊕ s)} c_{j ⊕ s}`. -/
noncomputable def shiftOp (s : Bitstring n) (φ : Bitstring n → ZMod 4) (c : Amp n) :
    Amp n := fun j => iPow (φ (j + s)) * c (j + s)

/-- `M` is a Hermitian involution: `φ k + φ (k ⊕ s) = 0`. Holds because `M` is a
unitary conjugate of `Z_q`, which is a Hermitian involution. -/
def SelfInverse (s : Bitstring n) (φ : Bitstring n → ZMod 4) : Prop :=
  ∀ k, φ k + φ (k + s) = 0

theorem add_shift_shift (s j : Bitstring n) : j + s + s = j := by
  rw [add_assoc, add_self, add_zero]

/-- **`M² = I`.** -/
theorem shiftOp_involutive {s : Bitstring n} {φ : Bitstring n → ZMod 4}
    (hφ : SelfInverse s φ) : Function.Involutive (shiftOp s φ) := by
  intro c
  funext j
  simp only [shiftOp, add_shift_shift]
  rw [← mul_assoc, ← iPow_add, add_comm (φ (j + s)) (φ j), hφ j, iPow_zero, one_mul]

/-- **`M† = M`**, stated as `⟨c, M d⟩ = ⟨M c, d⟩` over the full index type. -/
theorem shiftOp_selfAdjoint {s : Bitstring n} {φ : Bitstring n → ZMod 4}
    (hφ : SelfInverse s φ) (c d : Amp n) :
    ∑ j, (starRingEnd ℂ) (c j) * shiftOp s φ d j
      = ∑ j, (starRingEnd ℂ) (shiftOp s φ c j) * d j := by
  rw [← Equiv.sum_comp (xorRelabel s)
    (fun j => (starRingEnd ℂ) (shiftOp s φ c j) * d j)]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [shiftOp, xorRelabel_apply, add_shift_shift, map_mul, iPow_conj,
    neg_eq_of_add_eq_zero_right (hφ j)]
  ring

/-! ### The overlap and the Born probability -/

/-- The crate's `z_overlap_re`. -/
noncomputable def overlapRe (s : Bitstring n) (φ : Bitstring n → ZMod 4) (c : Amp n) :
    ℝ := ∑ k, rustTerm (φ k) (c k) (c (k + s))

/-- **The crate's overlap is `Re ⟨c, M c⟩`.** -/
theorem overlap_eq_inner {s : Bitstring n} {φ : Bitstring n → ZMod 4}
    (hφ : SelfInverse s φ) (c : Amp n) :
    overlapRe s φ c = (∑ j, (starRingEnd ℂ) (c j) * shiftOp s φ c j).re := by
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [rustTerm_eq]
  congr 1
  simp only [shiftOp, map_mul, iPow_conj, neg_eq_of_add_eq_zero_right (hφ k)]
  ring

/-- The outcome sign `(−1)^b`. -/
noncomputable def sgn (b : Bool) : ℂ := if b then -1 else 1

@[simp] theorem sgn_false : sgn false = 1 := rfl
@[simp] theorem sgn_true : sgn true = -1 := rfl

theorem sgn_sq (b : Bool) : sgn b * sgn b = 1 := by cases b <;> norm_num

/-- The measurement projector `P_b = (I + (−1)^b M)/2`. -/
noncomputable def proj (b : Bool) (s : Bitstring n) (φ : Bitstring n → ZMod 4)
    (c : Amp n) : Amp n := fun j => (c j + sgn b * shiftOp s φ c j) / 2

theorem proj_apply (b : Bool) (s : Bitstring n) (φ : Bitstring n → ZMod 4)
    (c : Amp n) (j : Bitstring n) :
    proj b s φ c j = (c j + sgn b * shiftOp s φ c j) / 2 := rfl

/-- **`P₀ + P₁ = I`** — the two outcomes are a complete measurement. -/
theorem proj_add (s : Bitstring n) (φ : Bitstring n → ZMod 4) (c : Amp n) (j) :
    proj false s φ c j + proj true s φ c j = c j := by
  simp only [proj_apply, sgn_false, sgn_true]
  ring

/-- **`P_b² = P_b`** — the case-a projection is idempotent (uses `M² = I`). -/
theorem proj_idem {s : Bitstring n} {φ : Bitstring n → ZMod 4}
    (hφ : SelfInverse s φ) (b : Bool) (c : Amp n) :
    proj b s φ (proj b s φ c) = proj b s φ c := by
  funext j
  have hph : iPow (φ (j + s)) * iPow (φ j) = 1 := by
    rw [← iPow_add, add_comm, hφ j, iPow_zero]
  have hlin : shiftOp s φ (proj b s φ c) j
      = (shiftOp s φ c j + sgn b * c j) / 2 := by
    simp only [shiftOp, proj, add_shift_shift]
    linear_combination (sgn b * c j / 2) * hph
  rw [proj_apply, hlin, proj_apply]
  linear_combination (c j / 4) * sgn_sq b

/-- **The Born rule.** For a normalized amplitude vector, the crate's
`prob_1 = 0.5 − 0.5 · z_overlap_re` is exactly `⟨c, P₁ c⟩`. -/
theorem probOne_eq {s : Bitstring n} {φ : Bitstring n → ZMod 4}
    (hφ : SelfInverse s φ) (c : Amp n)
    (hnorm : ∑ j, (starRingEnd ℂ) (c j) * c j = 1) :
    (∑ j, (starRingEnd ℂ) (c j) * proj true s φ c j).re
      = 1 / 2 - 1 / 2 * overlapRe s φ c := by
  have hsplit : ∑ j, (starRingEnd ℂ) (c j) * proj true s φ c j
      = (∑ j, (starRingEnd ℂ) (c j) * c j) * (1 / 2)
        - (∑ j, (starRingEnd ℂ) (c j) * shiftOp s φ c j) * (1 / 2) := by
    rw [Finset.sum_mul, Finset.sum_mul, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun j _ => ?_
    simp only [proj_apply, sgn_true]
    ring
  rw [hsplit, hnorm, overlap_eq_inner hφ]
  simp only [Complex.sub_re, Complex.mul_re, Complex.one_re, Complex.one_im]
  norm_num
  try ring

/-! ### The case-a projection map is the projector

The Rust keeps the `k`-bit-0 half `A` untouched and sends each `k`-bit-1 entry
`(idx, c)` to `(idx ⊕ s, i^{ψ idx} · c)`, summing the two streams. The pivot bit
is set in `s`, so `idx ↦ idx ⊕ s` maps the `k`-bit-1 half onto the `k`-bit-0
half and the merged output lives entirely on `A`. -/

/-- The raw case-a output: the `k`-bit-1 half is emptied, and each surviving
index `j` receives `c_j` plus the transformed partner `i^{ψ (j ⊕ s)} c_{j ⊕ s}`. -/
noncomputable def projectRaw (s : Bitstring n) (ψ : Bitstring n → ZMod 4)
    (kbit : Bitstring n → Bool) (c : Amp n) : Amp n :=
  fun j => if kbit j then 0 else c j + iPow (ψ (j + s)) * c (j + s)

/-- **The keep-`A`/transform-`B`/merge map is `2 · P_b`.** Given that the
projection's phase `ψ` is the overlap's `φ` twisted by the outcome sign
(`alpha = phase_decomp + 2·outcome`, i.e. `i^ψ = (−1)^b i^φ`), the case-a arm
computes exactly twice the Born projector on the surviving half — the factor `2`
being what the subsequent unconditional `normalize()` removes. -/
theorem projectRaw_eq_two_proj (b : Bool) (s : Bitstring n)
    (φ ψ : Bitstring n → ZMod 4) (hψ : ∀ k, iPow (ψ k) = sgn b * iPow (φ k))
    (kbit : Bitstring n → Bool) (c : Amp n) (j : Bitstring n) (hj : kbit j = false) :
    projectRaw s ψ kbit c j = 2 * proj b s φ c j := by
  simp only [projectRaw, hj, Bool.false_eq_true, reduceIte, proj_apply, shiftOp, hψ]
  ring

/-- The `k`-bit partition really is exchanged by the relabel: with the pivot bit
`q` set in `s`, `j ⊕ s` flips `j`'s bit at `q`, so exactly one of `{j, j ⊕ s}`
survives into the merged output. -/
theorem xor_flips_pivot (s : Bitstring n) (q : Fin n) (hs : s q = 1) (j : Bitstring n) :
    (j + s) q = j q + 1 := by
  simp only [Pi.add_apply, hs]

/-! ### Case b (`s = 0`): the projector with factor **1**

Everything above is case a. The crate's *other* measurement arm —
`measure_with_scratch`'s `stab_anticomm_bits == I::zero()` branch
(`crates/ppvm-tableau-2/src/measure.rs`), taken whenever `Z_q` is already a
stabilizer — is a structurally different formula, and `xor_flips_pivot` even
carries the hypothesis `s q = 1`, so the `s = 0` regime sits outside the
development so far. It is also the *majority* of the 85 measurements in the MSD
read-out sweep once the frame collapses, which is why the crate keeps it as an
explicit fast path rather than folding it into case a.

Specializing `proj` to `s = 0` settles the two behavioural contracts that arm
makes:

* `selfInverse_zero_phase_even` — the `SelfInverse` hypothesis degenerates to
  `φ k + φ k = 0`, i.e. `φ k ∈ {0, 2}`. This is the crate's `debug_assert!(
  phase_decomp == 0 || phase_decomp == 2, "Measurement result cannot be
  imaginary!")`, and the `if !phase.is_multiple_of(2) { continue }` guard in the
  overlap loop is then vacuous.
* `proj_zero_apply` — the projector is `c j ↦ c j` on the surviving set and `0`
  off it: **factor 1**, not the factor `2` of `projectRaw_eq_two_proj`. That is
  precisely what licenses case b applying *no* magnitude filter at all and
  calling `normalize()` only when the support shrank (`proj_zero_eq_self`: a
  case-b projection that drops nothing is the identity, hence exactly
  norm-preserving). A factor-2 slip here would be invisible in outcomes — the
  unconditional `normalize()` of case a would hide it — and would surface only
  much later as a wrong `⟨Z⟩`.
* `proj_zero_eq_caseB_retain` — the crate's `retain` predicate
  `(popcount(α ∧ destab_anticomm) odd) ⊕ outcome == (phase_decomp == 2)` is the
  indicator of exactly that surviving set. -/

theorem zmod4_add_self_eq_zero : ∀ a : ZMod 4, a + a = 0 → a = 0 ∨ a = 2 := by decide

/-- **With `s = 0` the phase is `ℤ/2`-valued**: `SelfInverse` degenerates to
`φ k + φ k = 0`, so `φ k ∈ {0, 2}` — the crate's "measurement result cannot be
imaginary" assertion. -/
theorem selfInverse_zero_phase_even {φ : Bitstring n → ZMod 4} (hφ : SelfInverse 0 φ)
    (k : Bitstring n) : φ k = 0 ∨ φ k = 2 :=
  zmod4_add_self_eq_zero (φ k) (by simpa using hφ k)

@[simp] theorem shiftOp_zero (φ : Bitstring n → ZMod 4) (c : Amp n) (j : Bitstring n) :
    shiftOp 0 φ c j = iPow (φ j) * c j := by
  simp only [shiftOp, add_zero]

/-- **Case b is the projector with factor 1.** For `s = 0`, `P_b` keeps `c j`
*untouched* when `(φ j = 2) = b` and kills it otherwise — no rescaling of the
survivors, unlike the factor `2` of `projectRaw_eq_two_proj` in case a. -/
theorem proj_zero_apply {φ : Bitstring n → ZMod 4} (hφ : SelfInverse 0 φ) (b : Bool)
    (c : Amp n) (j : Bitstring n) :
    proj b 0 φ c j = if decide (φ j = 2) = b then c j else 0 := by
  have h02 : ¬((0 : ZMod 4) = 2) := by decide
  rw [proj_apply, shiftOp_zero]
  rcases selfInverse_zero_phase_even hφ j with h | h <;> rw [h] <;> cases b <;>
    norm_num [sgn, h02]

/-- **A case-b projection that drops nothing is the identity** — hence exactly
norm-preserving, which is why the crate calls `normalize()` only when the support
actually shrank (behavioural contract #4) and applies no magnitude filter at all
(contract #3). -/
theorem proj_zero_eq_self {φ : Bitstring n → ZMod 4} (hφ : SelfInverse 0 φ) (b : Bool)
    (c : Amp n) (h : ∀ j, decide (φ j = 2) = b) : proj b 0 φ c = c := by
  funext j
  rw [proj_zero_apply hφ, if_pos (h j)]

/-- The crate's case-b keep-rule, as a `ℤ/2` identity: with
`φ α = phase_decomp + 2·⟨α, destab_anticomm⟩`, `zsign = (phase_decomp == 2)` and
`par = (popcount(α ∧ destab_anticomm) odd)`, the test `φ α = 2` holds iff
`zsign ⊕ par`, so `(φ α = 2) = outcome` iff `par ⊕ outcome = zsign`. -/
theorem caseB_retain_iff (zsign par outcome : Bool) :
    (decide (((if zsign then (2 : ZMod 4) else 0) + (if par then 2 else 0)) = 2) = outcome)
      ↔ (xor par outcome = zsign) := by
  cases zsign <;> cases par <;> cases outcome <;> decide

/-- **The crate's case-b `retain` predicate *is* the projector's support.** -/
theorem proj_zero_eq_caseB_retain {φ : Bitstring n → ZMod 4} (hφ : SelfInverse 0 φ)
    (zsign outcome : Bool) (par : Bitstring n → Bool) (c : Amp n)
    (hφα : ∀ α, φ α = (if zsign then 2 else 0) + (if par α then 2 else 0))
    (α : Bitstring n) :
    proj outcome 0 φ c α = if xor (par α) outcome = zsign then c α else 0 := by
  rw [proj_zero_apply hφ, hφα]
  by_cases h : xor (par α) outcome = zsign
  · rw [if_pos h, if_pos ((caseB_retain_iff zsign (par α) outcome).mpr h)]
  · rw [if_neg h, if_neg fun hc => h ((caseB_retain_iff zsign (par α) outcome).mp hc)]

end PPVM.Tableau.Projection
