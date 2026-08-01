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
