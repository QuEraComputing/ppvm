/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib
import PPVM.Algebra.GradedMap

/-!
# Noise channels and observable extraction

The Tier-3 targets from `sum/noise.rs` and the observable read-out, modeled on
the Pauli-basis coefficient vector.

* **Unital Pauli channel eigenvalue.** A Pauli channel acts diagonally in the
  Pauli basis: `P ↦ λ_P · P` with `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}`. Because
  `Σ_Q p_Q = 1`, this collapses to `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`.
* **Pauli-basis orthonormality.** `Tr(PQ)/2ⁿ = δ_{PQ}`: in the `C[K]` model the
  normalized trace pairing is `GradedMap.overlap`, and the Pauli basis vectors
  `single P 1` are orthonormal under it.
* **Zero-state read-out.** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P ∈ {I,Z}ⁿ} c_P`: only the diagonal
  (`X`-free) Paulis contribute, and each contributes its coefficient.
-/

namespace PPVM.Noise

open scoped BigOperators

variable {K : Type*} [Fintype K] [DecidableEq K]

/-- **Unital Pauli channel eigenvalue.** With `Σ_Q p_Q = 1`, the Pauli-transfer
eigenvalue `λ_P = Σ_Q p_Q (−1)^{ω(P,Q)}` equals `1 − 2·Σ_{Q anticommutes} p_Q`.
`anti P Q` is `ω(P,Q) = 1` (they anticommute), so `(−1)^{ω} = if anti then −1
else 1`. -/
theorem pauli_channel_eigenvalue (anti : K → K → Prop) [∀ P Q, Decidable (anti P Q)]
    (p : K → ℝ) (hp : ∑ Q, p Q = 1) (P : K) :
    ∑ Q, p Q * (if anti P Q then -1 else 1)
      = 1 - 2 * ∑ Q ∈ Finset.univ.filter (anti P ·), p Q := by
  have h : ∀ Q, p Q * (if anti P Q then (-1 : ℝ) else 1)
             = p Q - 2 * (if anti P Q then p Q else 0) := by
    intro Q; split <;> ring
  simp_rw [h]
  rw [Finset.sum_sub_distrib, hp, ← Finset.mul_sum, Finset.sum_filter]

/-- **Pauli-basis orthonormality** `Tr(PQ)/2ⁿ = δ_{PQ}`. In the `C[K]` model the
normalized trace is the `L3` pairing `GradedMap.overlap`, and the Pauli basis
`single P 1` is orthonormal: `⟪P, Q⟫ = δ_{PQ}`. -/
theorem pauli_orthonormal (P Q : K) :
    GradedMap.overlap (Finsupp.single P (1 : ℝ)) (Finsupp.single Q 1)
      = if P = Q then 1 else 0 := by
  rw [GradedMap.overlap, Finsupp.sum_single_index (by simp)]
  simp [Finsupp.single_apply, eq_comm]

/-- **Zero-state read-out** `⟨0ⁿ|ρ|0ⁿ⟩ = Σ_{P diagonal} c_P`. Each Pauli's
zero-state expectation is `1` for a diagonal (`X`-free, `I`/`Z`-only) Pauli and
`0` otherwise, so the expectation collapses to the coefficient sum over the
diagonal sector. -/
theorem overlap_with_zero (c : K → ℝ) (diag : K → Prop) [DecidablePred diag] :
    ∑ P, c P * (if diag P then 1 else 0) = ∑ P ∈ Finset.univ.filter diag, c P := by
  simp_rw [mul_ite, mul_one, mul_zero]
  rw [Finset.sum_filter]

end PPVM.Noise
