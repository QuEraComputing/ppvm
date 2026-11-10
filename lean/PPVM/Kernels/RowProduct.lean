/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Kernels.Model

/-!
The phase of the packed Pauli-word product in `ppvm-pauli-kernels`
(`crates/ppvm-pauli-kernels/src/product.rs`, `row_product`).

For words `P = ⊗ p_k` and `Q = ⊗ q_k` the product is `P · Q = iᵉ · ⊗ (p_k ⊕ q_k)`.
The kernel computes two masks per site, `sign` and `imag`, and returns
`e = (2 · popcount(sign) + popcount(imag)) mod 4`.
-/

namespace PPVM.Kernels

/-- The `sign` mask at one site, for `p = (a, b)` and `q = (c, d)`. -/
def siteSign (a b c d : Bool) : Bool :=
  (a && b && c && !d) || (a && !b && !c && d) || (!a && b && c && d)

/-- The `imag` mask at one site. -/
def siteImag (a b c d : Bool) : Bool :=
  (a && !b && d) || (a && !c && d) || (!a && b && c) || (b && c && !d)

/-- The phase exponent one site contributes. -/
def siteExp (a b c d : Bool) : Nat := 2 * (siteSign a b c d).toNat + (siteImag a b c d).toNat

/-- One site: `p · q = i^{siteExp} · (p ⊕ q)`, against the Pauli matrices. -/
theorem site_product : ∀ a b c d : Bool,
    Mat.Eq (Mat.mul (pauli a b) (pauli c d))
      (Mat.smul (GI.iPow (siteExp a b c d)) (pauli (a ^^ c) (b ^^ d))) := by
  decide

/-- Two sites: the exponents add across a tensor product (all 256 cases). -/
theorem two_site_product : ∀ a b a' b' c d c' d' : Bool,
    Mat.Eq (Mat.mul (pauli2 a b a' b') (pauli2 c d c' d'))
      (Mat.smul (GI.iPow (siteExp a b c d + siteExp a' b' c' d'))
        (pauli2 (a ^^ c) (b ^^ d) (a' ^^ c') (b' ^^ d'))) := by
  decide

/-- Any width: summing the per-site exponents is what the popcounts compute.
A site is `(a, b, c, d)`; `countP` is the popcount of the mask. -/
theorem popcount_exponent (sites : List (Bool × Bool × Bool × Bool)) :
    2 * sites.countP (fun s => siteSign s.1 s.2.1 s.2.2.1 s.2.2.2)
      + sites.countP (fun s => siteImag s.1 s.2.1 s.2.2.1 s.2.2.2)
    = (sites.map fun s => siteExp s.1 s.2.1 s.2.2.1 s.2.2.2).sum := by
  induction sites with
  | nil => simp
  | cons s rest ih =>
    rw [List.countP_cons, List.countP_cons, List.map_cons, List.sum_cons, ← ih]
    unfold siteExp
    cases siteSign s.1 s.2.1 s.2.2.1 s.2.2.2 <;>
      cases siteImag s.1 s.2.1 s.2.2.1 s.2.2.2 <;> simp <;> omega

end PPVM.Kernels
