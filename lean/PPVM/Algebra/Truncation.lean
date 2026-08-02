/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Normed.Unbundled.RingSeminorm
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import PPVM.Algebra.GradedMap

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

The `ℓ¹` bound is stated twice: once concretely over `ℝ`/`|·|`, and once over an
**arbitrary coefficient ring with an absolute value** `N` — which is the law the
Rust `Coefficient::magnitude` must satisfy for `CoefficientThreshold` to have any
error guarantee at all. `Complex<f64>`'s shipped `magnitude() = norm()` is an
instance.

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

/-! ### The `ℓ¹` bound over an arbitrary coefficient ring

`l1_bound` above is stated at `C = ℝ` with `N = |·|`. The Rust `Coefficient`
trait is generic in the coefficient ring and exposes only `magnitude() -> f64`,
documented as "a nonnegative magnitude … a property for a `Policy` to threshold"
— *stating no law*. But `CoefficientThreshold`'s keep-rule is only meaningful,
and the `ℓ¹` truncation error is only bounded, if `magnitude` is an **absolute
value** on the coefficient ring:

* `N x ≥ 0`,
* `N x = 0 ↔ x = 0`,
* `N (x + y) ≤ N x + N y` (subadditive),
* `N (x * y) = N x * N y` (multiplicative).

That bundle is Mathlib's `AbsoluteValue C ℝ`. A conforming-but-degenerate
`magnitude` (constant, or merely nonnegative) satisfies the *trait* while voiding
the bound, so this is a genuine implementation obligation, not a restatement.
-/

/-- **L1 truncation bound over any coefficient ring with an absolute value.**
For a coefficient ring `C` and an absolute value `N : AbsoluteValue C ℝ` — the
law `Coefficient::magnitude` must satisfy — with per-key expectations bounded by
`1` in `N`, the error from dropping the terms in `D` is at most the `ℓ¹` mass of
their coefficients: `N (Σ_{k∈D} c_k e_k) ≤ Σ_{k∈D} N (c_k)`.

`l1_bound` is the `C = ℝ`, `N = |·|` case; `l1_bound_norm` below is the
`Complex<f64>` case. Only subadditivity, multiplicativity, and nonnegativity are
used — exactly the absolute-value laws. -/
theorem l1_bound_abv {C : Type*} [Semiring C] (N : AbsoluteValue C ℝ)
    (e : K → C) (he : ∀ k, N (e k) ≤ 1) (c : K → C) (D : Finset K) :
    N (∑ k ∈ D, c k * e k) ≤ ∑ k ∈ D, N (c k) := by
  calc N (∑ k ∈ D, c k * e k)
      ≤ ∑ k ∈ D, N (c k * e k) := N.sum_le _ _
    _ = ∑ k ∈ D, N (c k) * N (e k) := by simp_rw [N.map_mul]
    _ ≤ ∑ k ∈ D, N (c k) * 1 :=
        Finset.sum_le_sum fun k _ =>
          mul_le_mul_of_nonneg_left (he k) (N.nonneg _)
    _ = ∑ k ∈ D, N (c k) := by simp

/-- **L1 truncation bound for a normed coefficient field**, e.g. `ℂ`. This is
`l1_bound_abv` at the absolute value `‖·‖`, which is what the shipped
`impl Coefficient for Complex<f64> { fn magnitude(&self) = self.norm() }`
computes — so `PauliSum<Complex<f64>>` + `CoefficientThreshold` is covered by a
machine-checked bound, not only the real case. -/
theorem l1_bound_norm {C : Type*} [NormedField C]
    (e : K → C) (he : ∀ k, ‖e k‖ ≤ 1) (c : K → C) (D : Finset K) :
    ‖∑ k ∈ D, c k * e k‖ ≤ ∑ k ∈ D, ‖c k‖ :=
  l1_bound_abv (NormedField.toAbsoluteValue C) e he c D

/-- The complex instantiation spelled out: the shipped `Complex<f64>`
coefficient domain has the `ℓ¹` truncation bound. -/
theorem l1_bound_complex (e : K → ℂ) (he : ∀ k, ‖e k‖ ≤ 1) (c : K → ℂ)
    (D : Finset K) : ‖∑ k ∈ D, c k * e k‖ ≤ ∑ k ∈ D, ‖c k‖ :=
  l1_bound_norm e he c D

/-- **A merely nonnegative (even multiplicative) `magnitude` is not enough.**
`N x = x²` is nonnegative, vanishes only at `0`, and is multiplicative — it meets
every property the Rust trait doc currently states — yet it is *not subadditive*
and the `ℓ¹` bound fails for it: two unit expectations with unit coefficients give
a true "error" of `N 2 = 4` against a claimed bound of `2`. So subadditivity is
load-bearing, and `Coefficient::magnitude` must be documented as an absolute
value, not merely as "a nonnegative magnitude". -/
theorem l1_bound_needs_subadditive :
    ¬ ∀ (N : ℝ → ℝ), (∀ x, 0 ≤ N x) → (∀ x y, N (x * y) = N x * N y) →
      ∀ (e : Bool → ℝ), (∀ k, N (e k) ≤ 1) → ∀ (c : Bool → ℝ) (D : Finset Bool),
        N (∑ k ∈ D, c k * e k) ≤ ∑ k ∈ D, N (c k) := by
  intro h
  have := h (fun x => x ^ 2) (fun x => sq_nonneg x) (fun x y => by ring)
    (fun _ => 1) (by intro k; norm_num) (fun _ => 1) Finset.univ
  norm_num [Fintype.sum_bool] at this

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

/-! ### The preserved-key post-filter is a **widened** keep-rule

`Sum::truncate` (`ppvm-pauli-sum-2/src/sum.rs`, ported from
`ppvm-pauli-sum/src/sum/data.rs`) is not just the policy: with a non-empty
keep-set `P` (`preserve_strings`) it runs three steps —

1. **snapshot** the *pre-truncate* coefficients of the keys in `P` that are in
   the support;
2. run the configured policy **verbatim** (a `retain`);
3. **restore**, through `AddTerm::add_term`, every snapshotted key the policy
   dropped, guarded by a membership test.

Both guards are load-bearing and neither is visible in the types: `add_term`
*accumulates*, so without the membership test a survivor would have its
coefficient **doubled**, and if the snapshot were taken after the policy ran it
would restore a post-truncate residue instead of the original coefficient.

`truncate_preserve_eq_widened_retain` collapses the whole composite to a single
pass with the **widened** keep-rule `keep k c ∨ k ∈ P`. That is the specification
`preserve_strings` advertises, and it has two immediate corollaries:

* `truncatePreserve_apply_of_mem` — a preserved key keeps **exactly** its
  pre-truncate coefficient (this is old's conservation test: `Σᵢ Zᵢ` under
  repeated `rxx`/`ryy` with aggressive truncation holds every single-`Z`
  coefficient at exactly `1.0`);
* the dropped set is `D \ P`, so the `ℓ¹` bound above applies to `D \ P` — the
  error guarantee a caller who passes a keep-set actually gets.

As everywhere in this model the support is the canonical (zero-free) one
(`GradedMap.reduce_structural`), so the statement is about **coefficients**: a
backend that additionally stores explicit zeros — which `Sum` does, since it
never drops them on its own — agrees with it pointwise, while whether such a key
is *listed* is that backend's `reduce` question, not this one.
-/

section Preserve

open PPVM.GradedMap

variable {C : Type*} [AddCommMonoid C] [DecidableEq K]

/-- Step 1 — the snapshot: the preserved keys' **pre-truncate** coefficients
(only those actually in the support, exactly as the Rust `storage.get` scan). -/
noncomputable def snapshot (P : Finset K) (A : CMap K C) : CMap K C :=
  Finsupp.filter (· ∈ P) A

/-- The whole three-step `Sum::truncate`: policy `keep`, then restore each
snapshotted key the policy dropped, at its pre-truncate coefficient. The
difference `\ (retain keep A).support` is the `storage.get(&key).is_none()`
guard — without it the restore would *add* to a survivor. -/
noncomputable def truncatePreserve (keep : K → C → Bool) (P : Finset K) (A : CMap K C) :
    CMap K C :=
  retain keep A
    + ∑ k ∈ (snapshot P A).support \ (retain keep A).support, Finsupp.single k (A k)

/-- **Snapshot–truncate–restore is one widened `retain`.**
`truncatePreserve keep P A = retain (fun k c => keep k c ∨ k ∈ P) A`: the
restored coefficient is the pre-truncate one, and no surviving key is doubled. -/
theorem truncate_preserve_eq_widened_retain (keep : K → C → Bool) (P : Finset K)
    (A : CMap K C) :
    truncatePreserve keep P A = retain (fun k c => keep k c || decide (k ∈ P)) A := by
  classical
  ext j
  have hsum : (∑ k ∈ (snapshot P A).support \ (retain keep A).support,
      Finsupp.single k (A k)) j
      = if j ∈ (snapshot P A).support \ (retain keep A).support then A j else 0 := by
    rw [Finsupp.finsetSum_apply]
    simp only [Finsupp.single_apply]
    exact Finset.sum_ite_eq' _ _ _
  rw [truncatePreserve, Finsupp.add_apply, hsum]
  simp only [retain_apply, Finset.mem_sdiff, Finsupp.mem_support_iff, snapshot,
    Finsupp.filter_apply]
  by_cases hA : A j = 0
  · simp [hA]
  · by_cases hP : j ∈ P <;> cases hk : keep j (A j) <;> simp [hA, hP]

/-- **A preserved key keeps its pre-truncate coefficient exactly.** Whatever the
policy decides, `truncate()` leaves every key of the keep-set at the value it had
before the call — old's `tests/preserve.rs` conservation property. -/
theorem truncatePreserve_apply_of_mem (keep : K → C → Bool) (P : Finset K)
    (A : CMap K C) {j : K} (hj : j ∈ P) : truncatePreserve keep P A j = A j := by
  rw [truncate_preserve_eq_widened_retain, retain_apply, if_pos (by simp [hj])]

/-- **An empty keep-set is the policy, verbatim** — the hot-path short-circuit
(`if self.preserve.is_empty()`) is exact, not an approximation. -/
theorem truncatePreserve_empty (keep : K → C → Bool) (A : CMap K C) :
    truncatePreserve keep (∅ : Finset K) A = retain keep A := by
  rw [truncate_preserve_eq_widened_retain]
  congr 1
  funext k c
  simp

end Preserve

end PPVM.Truncation
