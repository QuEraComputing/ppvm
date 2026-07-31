/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib

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

/-- **L2 truncation bound.** The squared error is bounded by the product of the
dropped coefficients' and expectations' `ℓ²` masses (Cauchy–Schwarz) — the bound
the stabilizer-tableau truncation uses. -/
theorem l2_bound (e : K → ℝ) (c : K → ℝ) (D : Finset K) :
    (∑ k ∈ D, c k * e k) ^ 2 ≤ (∑ k ∈ D, c k ^ 2) * (∑ k ∈ D, e k ^ 2) :=
  Finset.sum_mul_sq_le_sq_mul_sq D c e

/-- **The backend cutoff mismatch.** `PauliSum` keeps a term when
`threshold ≤ |c|`; the tableau keeps it when `threshold < |c|`. At a coefficient
exactly at the threshold the two rules disagree — a term kept by one backend is
dropped by the other. -/
theorem cutoff_mismatch : ∃ t c : ℝ, (t ≤ |c|) ≠ (t < |c|) :=
  ⟨1, 1, by norm_num⟩

end PPVM.Truncation
