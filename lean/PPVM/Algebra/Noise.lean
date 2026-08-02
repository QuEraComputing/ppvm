/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Data.Real.Basic
import PPVM.Algebra.GradedMap
import PPVM.Pauli.Symplectic

/-!
# Noise channels and observable extraction

The Tier-3 targets from `sum/noise.rs` and the observable read-out, modeled on
the Pauli-basis coefficient vector.

* **Unital Pauli channel eigenvalue.** A Pauli channel acts diagonally in the
  Pauli basis: `P ↦ λ_P · P` with `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}`. Because
  `Σ_Q p_Q = 1`, this collapses to `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`.
  The Pauli-specific corollary ties `anticommute` to the actual symplectic form
  `PPVM.Symplectic.omega`.
* **Pauli-basis orthonormality.** In the `C[K]` model the L3 `overlap` pairing
  plays the role of the normalized trace `Tr(PQ)/2ⁿ`; the Pauli basis vectors
  `single P 1` are orthonormal under it. (The `2ⁿ`-normalized matrix trace itself
  is not constructed *here*; it is now built in `PPVM.PauliMatrix`, where
  `trace_tensorPauli_mul` proves the genuine matrix identity
  `Tr(g(p) g(q)) = 2ⁿ δ` and `overlap_eq_trace_div` proves `overlap` *is*
  `Tr(Â B̂)/2ⁿ`. This file's statement is the abstract-key model form.)
* **Zero-state read-out.** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P ∈ {I,Z}ⁿ} c_P`: only the diagonal
  (`X`-free) Paulis contribute; the corollary uses the concrete `X`-free predicate.
-/

namespace PPVM.Noise

open scoped BigOperators

/-! ### Unital Pauli channel eigenvalue -/

/-- **Unital Pauli channel eigenvalue** (general form). With `Σ_Q p_Q = 1`, the
Pauli-transfer eigenvalue `λ_P = Σ_Q p_Q (−1)^{[anti]}` equals
`1 − 2·Σ_{anti} p_Q`. -/
theorem pauli_channel_eigenvalue {K : Type*} [Fintype K] (anti : K → K → Prop)
    [∀ P Q, Decidable (anti P Q)] (p : K → ℝ) (hp : ∑ Q, p Q = 1) (P : K) :
    ∑ Q, p Q * (if anti P Q then -1 else 1)
      = 1 - 2 * ∑ Q ∈ Finset.univ.filter (anti P ·), p Q := by
  have h : ∀ Q, p Q * (if anti P Q then (-1 : ℝ) else 1)
             = p Q - 2 * (if anti P Q then p Q else 0) := by
    intro Q; split <;> ring
  simp_rw [h]
  rw [Finset.sum_sub_distrib, hp, ← Finset.mul_sum, Finset.sum_filter]

/-- **Pauli channel eigenvalue formula, tied to the symplectic form.** Here `anti`
is genuine anticommutation `ω(P,Q) = 1`. This is the arithmetic identity behind
`λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`; the channel superoperator and its
diagonalization are not constructed here — only the eigenvalue's algebraic form. -/
theorem pauli_channel_eigenvalue_omega {m : ℕ} (p : Symplectic.Sp m → ℝ)
    (hp : ∑ Q, p Q = 1) (P : Symplectic.Sp m) :
    ∑ Q, p Q * (if Symplectic.omega P Q = 1 then -1 else 1)
      = 1 - 2 * ∑ Q ∈ Finset.univ.filter (fun Q => Symplectic.omega P Q = 1), p Q :=
  pauli_channel_eigenvalue (fun P Q => Symplectic.omega P Q = 1) p hp P

/-! ### The eigenvalue is **contractive**, and the channel never grows a key

`pauli_channel_eigenvalue` above gives the *formula* for `λ_P` but says nothing
about its size. The size is what the implementation actually cites: the
`Sum::scale_by_key` fast path (`crates/ppvm-pauli-sum-2/src/sum.rs`,
`PauliError`/`Depolarizing`/…) runs **no** truncation pass and **no** Pauli-weight
re-check after a channel, on the grounds that "the channel is contractive
(`|λ_P| ≤ 1`), so it can never grow a key's Pauli weight". Two separate facts are
being claimed there, and both are proved below:

* `|λ_P| ≤ 1` for a *sub-stochastic* probability vector (`p ≥ 0`, `Σ p ≤ 1` — the
  single-qubit `PauliError` takes `[p_X, p_Y, p_Z]` with an implicit
  `p_I = 1 − Σ`), hence the diagonal channel is an `ℓ¹` **contraction**
  (`l1_contractive`). That is also the missing hypothesis for composing
  `PPVM.Truncation.l1_bound` across a *noisy* circuit: the per-truncation `ℓ¹`
  error bound only telescopes over a long propagation if every intervening
  channel is non-expanding in `ℓ¹`.
* The diagonal channel's support only shrinks (`scaleByKey_support_subset`) and
  never moves a key, so every surviving key's Pauli weight is one it already had
  — nothing for a weight policy to re-check.

`eigenvalue_abs_le_one_needs_substochastic` shows the hypothesis is load-bearing:
a caller passing an over-normalized `[p_X, p_Y, p_Z]` breaks the bound (and with
it both the skipped truncation and the composed error bound), so `Σ p ≤ 1` is a
real precondition on the Rust channel constructors, not a formality. -/

/-- **The Pauli-transfer eigenvalue is contractive.** For a sub-stochastic
probability vector (`p Q ≥ 0` and `Σ_Q p Q ≤ 1`), the eigenvalue
`λ_P = 1 − 2·Σ_{Q anti P} p_Q` of `pauli_channel_eigenvalue` satisfies
`|λ_P| ≤ 1`. (The anticommuting mass is between `0` and `1`, so `λ_P ∈ [−1, 1]`.) -/
theorem pauli_channel_eigenvalue_abs_le_one {K : Type*} [Fintype K] (anti : K → K → Prop)
    [∀ P Q, Decidable (anti P Q)] (p : K → ℝ) (hp0 : ∀ Q, 0 ≤ p Q) (hp : ∑ Q, p Q ≤ 1)
    (P : K) :
    |1 - 2 * ∑ Q ∈ Finset.univ.filter (anti P ·), p Q| ≤ 1 := by
  have h0 : 0 ≤ ∑ Q ∈ Finset.univ.filter (anti P ·), p Q :=
    Finset.sum_nonneg fun Q _ => hp0 Q
  have h1 : ∑ Q ∈ Finset.univ.filter (anti P ·), p Q ≤ 1 :=
    le_trans (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      (fun Q _ _ => hp0 Q)) hp
  rw [abs_le]
  constructor <;> linarith

/-- **The same bound, tied to the symplectic form** — the eigenvalue of
`pauli_channel_eigenvalue_omega` is contractive. -/
theorem pauli_channel_eigenvalue_omega_abs_le_one {m : ℕ} (p : Symplectic.Sp m → ℝ)
    (hp0 : ∀ Q, 0 ≤ p Q) (hp : ∑ Q, p Q ≤ 1) (P : Symplectic.Sp m) :
    |1 - 2 * ∑ Q ∈ Finset.univ.filter (fun Q => Symplectic.omega P Q = 1), p Q| ≤ 1 :=
  pauli_channel_eigenvalue_abs_le_one (fun P Q => Symplectic.omega P Q = 1) p hp0 hp P

/-- **A contractive diagonal channel is an `ℓ¹` contraction on `C[K]`.** If every
transfer eigenvalue satisfies `|λ_k| ≤ 1` then `‖Λ A‖₁ ≤ ‖A‖₁` on any finite set
of keys. This is what lets the `ℓ¹` truncation bound (`PPVM.Truncation.l1_bound`)
compose across a noisy propagation: the channels between two truncations cannot
inflate the mass a later truncation is measured against. -/
theorem l1_contractive {K : Type*} (lam c : K → ℝ) (hlam : ∀ k, |lam k| ≤ 1)
    (D : Finset K) :
    ∑ k ∈ D, |lam k * c k| ≤ ∑ k ∈ D, |c k| := by
  refine Finset.sum_le_sum fun k _ => ?_
  rw [abs_mul]
  exact mul_le_of_le_one_left (abs_nonneg _) (hlam k)

/-- The diagonal channel on `C[K]`: rescale each coefficient by its own key's
eigenvalue, leaving the key alone. This is `Sum::scale_by_key`'s action
(`f(&k, &mut c)` over `iter_mut`: no key ever moves). -/
noncomputable def scaleByKey {K C : Type*} [DecidableEq K] [Semiring C] (lam : K → C)
    (A : GradedMap.CMap K C) : GradedMap.CMap K C :=
  A.sum fun k a => Finsupp.single k (lam k * a)

/-- **A diagonal channel never introduces a key.** Its support is contained in the
input's, so every surviving term's Pauli weight is one the map already carried —
the formal content of "`pauli_error` runs no weight re-check". (It is only a
*subset*: a zero eigenvalue kills a coefficient. The Rust backend deliberately
keeps such a term in its map with coefficient `0`, which refines this ideal
`Finsupp` model, where a zero coordinate is structurally absent —
`GradedMap.reduce_structural`.) -/
theorem scaleByKey_support_subset {K C : Type*} [DecidableEq K] [Semiring C] (lam : K → C)
    (A : GradedMap.CMap K C) : (scaleByKey lam A).support ⊆ A.support := by
  classical
  refine (Finsupp.support_sum).trans ?_
  refine Finset.biUnion_subset.2 fun k hk => ?_
  exact (Finsupp.support_single_subset).trans (by simpa using hk)

/-- **The sub-stochastic hypothesis is load-bearing.** Drop `Σ_Q p_Q ≤ 1` and the
bound fails: the nonnegative vector `p ≡ 1` on a two-element alphabet gives
`|1 − 2·2| = 3 > 1`. So an over-normalized `[p_X, p_Y, p_Z]` silently breaks both
the skipped weight re-check and the composed `ℓ¹` truncation bound; the Rust
channel constructors owe the precondition. -/
theorem eigenvalue_abs_le_one_needs_substochastic :
    ¬ ∀ p : Bool → ℝ, (∀ Q, 0 ≤ p Q) → |1 - 2 * ∑ Q, p Q| ≤ 1 := by
  intro h
  have := h (fun _ => 1) (fun _ => by norm_num)
  rw [Fintype.sum_bool] at this
  norm_num at this

/-! ### Pauli-basis orthonormality (the `overlap` pairing) -/

/-- **Pauli-basis orthonormality** in the `C[K]` model: the basis vectors
`single P 1` are orthonormal under the L3 `overlap` pairing, `⟪P, Q⟫ = δ_{PQ}`.
This is the model-level form of `Tr(PQ)/2ⁿ = δ_{PQ}` over an *abstract* key set,
with `overlap` standing in for the normalized trace. On the concrete Pauli key
the matrix identity itself is now proved: `PPVM.PauliMatrix.trace_tensorPauli_mul`
(and `overlap_eq_trace_div`, `Tr(Â B̂) = 2ⁿ ⟪A,B⟫`). -/
theorem overlap_single_single {K : Type*} [DecidableEq K] (P Q : K) :
    GradedMap.overlap (Finsupp.single P (1 : ℝ)) (Finsupp.single Q 1)
      = if P = Q then 1 else 0 := by
  rw [GradedMap.overlap, Finsupp.sum_single_index (by simp)]
  simp [Finsupp.single_apply, eq_comm]

/-! ### Zero-state read-out -/

/-- **Zero-state read-out** (general form): with a per-Pauli `diag` indicator,
`Σ_P c_P·[diag P] = Σ_{diag} c_P`. -/
theorem overlap_with_zero {K : Type*} [Fintype K] (c : K → ℝ) (diag : K → Prop)
    [DecidablePred diag] :
    ∑ P, c P * (if diag P then 1 else 0) = ∑ P ∈ Finset.univ.filter diag, c P := by
  simp_rw [mul_ite, mul_one, mul_zero]
  rw [Finset.sum_filter]

/-- **Zero-state read-out, concrete.** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P X-free} c_P`. The physics
input — a Pauli's zero-state expectation `⟨0|P|0⟩` is `1` exactly for the diagonal
(`X`-free, `I`/`Z`-only) Paulis — is the *modeling assumption* (the `if … then 1
else 0` factor); the theorem is the resulting collapse to the coefficient sum over
the concrete `X`-free sector (`diag P := ∀ i, (P i).1 = 0`). -/
theorem overlap_with_zero_xfree {m : ℕ} (c : Symplectic.Sp m → ℝ) :
    ∑ P, c P * (if (∀ i, (P i).1 = 0) then 1 else 0)
      = ∑ P ∈ Finset.univ.filter (fun P => ∀ i, (P i).1 = 0), c P :=
  overlap_with_zero c (fun P => ∀ i, (P i).1 = 0)

end PPVM.Noise
