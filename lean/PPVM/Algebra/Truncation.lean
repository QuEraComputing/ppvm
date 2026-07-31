/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# Truncation error bounds

`GradedMap.truncMag_not_additive` showed truncation is not an algebra operation.
Here are the quantitative error bounds it satisfies — the guarantees that justify
dropping terms (`docs/design/…`, Tier-3 truncation targets).

A sparse observable is `⟨O⟩ = Σ_P c_P ⟨P⟩` where each Pauli expectation obeys
`|⟨P⟩| ≤ 1`. Truncation drops a set `D` of terms; the incurred error is the
dropped contribution `Σ_{P∈D} c_P ⟨P⟩`. We bound it two ways:

* **L1 (`PauliSum`):** `|error| ≤ Σ_{P∈D} |c_P|` — triangle inequality.
* **L2 (tableau):** `error² ≤ (Σ_{P∈D} c_P²)(Σ_{P∈D} ⟨P⟩²)` — Cauchy–Schwarz.

We also pin down the design-noted **`<` vs `≥` cutoff mismatch** between the two
backends: at a coefficient exactly equal to the threshold, the two keep-rules
disagree.
-/

namespace PPVM.Truncation

variable {K : Type*}

/-- **L1 truncation bound.** With per-key expectations `e` bounded by `1`, the
error from dropping the terms in `D` is at most the `ℓ¹` mass of their
coefficients: `|Σ_{k∈D} c_k e_k| ≤ Σ_{k∈D} |c_k|`. -/
theorem l1_bound (e : K → ℝ) (he : ∀ k, |e k| ≤ 1) (c : K → ℝ) (D : Finset K) :
    |∑ k ∈ D, c k * e k| ≤ ∑ k ∈ D, |c k| := by
  calc |∑ k ∈ D, c k * e k|
      ≤ ∑ k ∈ D, |c k * e k| := Finset.abs_sum_le_sum_abs _ _
    _ = ∑ k ∈ D, |c k| * |e k| := by simp_rw [abs_mul]
    _ ≤ ∑ k ∈ D, |c k| * 1 :=
        Finset.sum_le_sum fun k _ =>
          mul_le_mul_of_nonneg_left (he k) (abs_nonneg _)
    _ = ∑ k ∈ D, |c k| := by simp

/-- **L2 truncation bound.** The squared error from dropping the terms in `D` is
bounded by the product of the dropped coefficients' and expectations' `ℓ²`
masses. This is Cauchy–Schwarz specialized to the dropped set `D` — that is
exactly the `ℓ²` truncation bound the stabilizer-tableau path uses. -/
theorem l2_bound (e : K → ℝ) (c : K → ℝ) (D : Finset K) :
    (∑ k ∈ D, c k * e k) ^ 2 ≤ (∑ k ∈ D, c k ^ 2) * (∑ k ∈ D, e k ^ 2) :=
  Finset.sum_mul_sq_le_sq_mul_sq D c e

/-- With expectations bounded (`e² ≤ 1`, from `|⟨P⟩| ≤ 1`), the squared error is
bounded by the dropped coefficients' `ℓ²` mass times the number of dropped terms
— the truncation-specific form of the `ℓ²` bound. -/
theorem l2_bound_normalized (e : K → ℝ) (he : ∀ k, e k ^ 2 ≤ 1) (c : K → ℝ) (D : Finset K) :
    (∑ k ∈ D, c k * e k) ^ 2 ≤ (∑ k ∈ D, c k ^ 2) * D.card := by
  refine (l2_bound e c D).trans (mul_le_mul_of_nonneg_left ?_ ?_)
  · calc ∑ k ∈ D, e k ^ 2 ≤ ∑ _k ∈ D, (1 : ℝ) := Finset.sum_le_sum fun k _ => he k
      _ = D.card := by simp
  · exact Finset.sum_nonneg fun k _ => sq_nonneg _

/-- **The backend cutoff mismatch, for every threshold.** `PauliSum` keeps a term
when `threshold ≤ |c|`; the tableau keeps it when `threshold < |c|`. At *any*
threshold `t ≥ 0`, a coefficient with `|c| = t` (take `c = t`) is kept by the
first rule and dropped by the second — the two backends disagree on the boundary. -/
theorem cutoff_mismatch (t : ℝ) (ht : 0 ≤ t) : (t ≤ |t|) ≠ (t < |t|) := by
  simp [abs_of_nonneg ht]

end PPVM.Truncation
