/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.NumberTheory.Zsqrtd.GaussianInt
import Mathlib.LinearAlgebra.Matrix.Notation
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Algebra.BigOperators.Ring.Finset
import PPVM.Pauli.Phase
import PPVM.Pauli.Word
import PPVM.Algebra.Twisted

/-!
# The Pauli phase, validated against actual `ℤ[i]` matrices

`PPVM.PauliPhase.phaseExp_eq_ref` proves the packed `sign`/`imag` booleans equal
a *formula* `phaseRef`. That formula was derived by hand, so on its own it is only
an *algebraically* independent reference. This file removes that caveat by
building the genuine single-qubit Pauli matrices over the Gaussian integers
`ℤ[i]` (which have decidable equality) and proving `phaseExp` is exactly the
base-`i` exponent of the real matrix product — a *model*-level check, closed by
`decide`.
-/

namespace PPVM.PauliMatrix

open PPVM.PauliPhase

/-- The imaginary unit `i ∈ ℤ[i]`. -/
def iU : GaussianInt := ⟨0, 1⟩

@[simp] theorem iU_sq : iU * iU = -1 := by decide

/-- The physical single-qubit Pauli `g(x,z) = iˣᶻ · Xˣ · Zᶻ` as a genuine 2×2
matrix over `ℤ[i]`: `I`, `X`, `Z`, and `Y = iXZ`. -/
def pauliMat : Bool → Bool → Matrix (Fin 2) (Fin 2) GaussianInt
  | false, false => 1
  | true,  false => !![0, 1; 1, 0]
  | false, true  => !![1, 0; 0, -1]
  | true,  true  => !![0, -iU; iU, 0]

/-- Sanity check on the convention: `Y = i · X · Z`. -/
theorem pauliMat_Y : pauliMat true true = iU • (pauliMat true false * pauliMat false true) := by
  decide

/-- **`phaseExp` is the exponent of the real matrix product.** For all 16 bit
patterns, `g(a,b) · g(c,d) = i^{phaseExp(a,b,c,d)} · g(a⊕c, b⊕d)` as honest 2×2
`ℤ[i]` matrices. Combined with `phaseExp_eq_ref`, this makes `phaseRef` a
*model*-grounded reference, not just a second formula. -/
theorem pauliMat_mul : ∀ a b c d : Bool,
    pauliMat a b * pauliMat c d
      = iU ^ (phaseExp a b c d).val • pauliMat (xor a c) (xor b d) := by
  decide

/-! ### The loosened L4 coefficient bound is realizable over an exact ring

Formalizing L4's twisted product (`PPVM.Twisted.tmul_assoc`) showed it is
associative over *any* commutative ring with a fourth root of unity — the design
therefore loosens its `Multiply` bound from `ComplexCoefficient` to a primitive
fourth root of unity. Here is the point made concrete: `ℤ[i]` is an **exact**
ring (no floats), `iU` is its fourth root, and the twisted Pauli product is
associative over it. So exact / symbolic Pauli multiplication is admissible. -/

/-- `iU` is a genuine fourth root of unity in `ℤ[i]`. -/
theorem iU_pow_four : iU ^ 4 = 1 := by decide

/-- The twisted `key_mul` product is associative over the exact ring `ℤ[i]`. -/
example (a b c : Twisted.Mono GaussianInt) :
    Twisted.tmul iU (Twisted.tmul iU a b) c = Twisted.tmul iU a (Twisted.tmul iU b c) :=
  Twisted.tmul_assoc iU iU_pow_four a b c

/-! ### `ℤ[i]` also realizes the L3 `Conjugate` capability

The sesquilinear `Pair::hermitian_overlap` (`PPVM.GradedMap.hermitianOverlap`)
needs the design's `Conjugate` (`*`-ring) capability on the coefficient ring.
`ℤ[i]` supplies it: `StarRing GaussianInt` is genuine complex conjugation, an
involutive ring involution whose action on the imaginary unit is `conj(i) = −i` —
the design's `Conjugate` law, again on an exact ring (no floats). -/

/-- `ℤ[i]` is a `*`-ring: the `Conjugate` capability is `StarRing.star`. -/
example : StarRing GaussianInt := inferInstance

/-- **The `Conjugate` law `conj(i) = −i` on `ℤ[i]`.** So `hermitian_overlap`'s
conjugation is realized by an exact ring, matching the design's `Conjugate`
requirement `conj(i) == −i` under `ImaginaryUnit`. -/
theorem star_iU : star iU = -iU := by decide

/-! ### The `n`-qubit phase, grounded against genuine tensor-product matrices

`PPVM.PauliWord.phaseExpN` is the packed multi-qubit phase exponent — the
`(2*sign_count + imag_count) mod 4` that `KeyProduct::key_mul` accumulates — and
it is *defined* as the **sum** over qubits of the single-qubit `phaseExp`. On its
own that leaves a gap: the object `key_mul` actually computes is an `n`-fold
tensor-product Pauli operator, and the equality between "sum of per-qubit
exponents" and "base-`i` exponent of the real `2ⁿ×2ⁿ` matrix product" was
*assumed* (the phase-multiplicativity step `∏ᵢ iᵏⁱ = i^{Σ kᵢ}`). This section
closes it: it builds the genuine `n`-fold tensor product of the single-qubit
`ℤ[i]` matrices and proves `phaseExpN` is exactly the base-`i` exponent of the
real matrix product, end-to-end. -/

open PPVM.PauliWord

variable {n : ℕ}

/-- The genuine `n`-fold tensor (Kronecker) product of the single-qubit Pauli
matrices `g(p i)`, as a `2ⁿ×2ⁿ` matrix over `ℤ[i]` indexed by multi-indices
`Fin n → Fin 2`: the `(r, c)` entry is the product of the per-qubit entries. -/
def tensorPauli (p : Word n) : Matrix (Fin n → Fin 2) (Fin n → Fin 2) GaussianInt :=
  fun r c => ∏ i, pauliMat (p i).1 (p i).2 (r i) (c i)

/-- The phase-multiplicativity of the tensor product, `∏ᵢ iᵏⁱ = i^{Σ kᵢ}`, made
precise as the monoid hom `(ℤ/4ℤ, +) → (ℤ[i], ·)` sending `k ↦ iᵏ`. Well defined
because `iU` is a fourth root of unity (`iU_pow_four`); the two homomorphism laws
close by `decide` over the four residues. -/
def iuPow : Multiplicative (ZMod 4) →* GaussianInt where
  toFun k := iU ^ (Multiplicative.toAdd k).val
  map_one' := by decide
  map_mul' x y := by
    suffices h : ∀ a b : ZMod 4, iU ^ (a + b).val = iU ^ a.val * iU ^ b.val from
      h (Multiplicative.toAdd x) (Multiplicative.toAdd y)
    decide

/-- **Phase-multiplicativity.** The product of the per-qubit phases `iᵏⁱ` equals
`i` raised to the packed sum `phaseExpN` — the step the `Word.lean` scope note
had left assumed. This is `iuPow`'s `map_prod` combined with the definition of
`phaseExpN` as a sum. -/
theorem prod_iuPow (p q : Word n) :
    (∏ i, iU ^ (phaseExp (p i).1 (p i).2 (q i).1 (q i).2).val)
      = iU ^ (phaseExpN p q).val := by
  have h : ∀ i, iU ^ (phaseExp (p i).1 (p i).2 (q i).1 (q i).2).val
      = iuPow (Multiplicative.ofAdd (phaseExp (p i).1 (p i).2 (q i).1 (q i).2)) :=
    fun _ => rfl
  simp only [h]
  rw [← map_prod, ← ofAdd_sum]
  rfl

/-- **`phaseExpN` is the exponent of the real tensor-product matrix product.**
For all `n`-qubit words `p q`, the genuine `2ⁿ×2ⁿ` product
`g(p)·g(q) = i^{phaseExpN p q} · g(p·q)` as honest `ℤ[i]` matrices. This ties the
packed multi-qubit phase `KeyProduct::key_mul` computes back to real multi-qubit
Pauli operators, using the single-qubit grounding `pauliMat_mul` and the
tensor-product phase-multiplicativity `prod_iuPow`. -/
theorem tensorPauli_mul (p q : Word n) :
    tensorPauli p * tensorPauli q
      = iU ^ (phaseExpN p q).val • tensorPauli (mulWord p q) := by
  refine Matrix.ext fun r c => ?_
  -- Reduce the LHS entry to the product of per-qubit `2×2` matrix products
  -- (the tensor mixed-product property).
  have hL : (tensorPauli p * tensorPauli q) r c
      = ∏ i, (pauliMat (p i).1 (p i).2 * pauliMat (q i).1 (q i).2) (r i) (c i) := by
    simp only [tensorPauli, Matrix.mul_apply]
    rw [Fintype.prod_sum]
    simp only [Finset.prod_mul_distrib]
  rw [hL]
  -- Substitute the single-qubit grounding, pull out the scalar per qubit, then
  -- collapse the product of phases via `prod_iuPow`.
  simp only [pauliMat_mul, Matrix.smul_apply, smul_eq_mul]
  rw [Finset.prod_mul_distrib, prod_iuPow]
  simp only [tensorPauli, mulWord, mulBits]

/-! ### `Pair::overlap` really is the normalized trace `Tr(A B)/2ⁿ`

`PPVM.GradedMap.overlap` is the formal bilinear form `∑ₖ aₖ bₖ`, and the design
(and `graded.rs`) call it "the symmetric bilinear Hilbert–Schmidt trace pairing
`⟨A,B⟩ = Tr(A B)/2ⁿ`" — but that identification was, until here, *asserted*: the
`GradedMap` docstring conceded that the trace "is not itself constructed in
Lean", so `Noise.overlap_single_single` was only the *model* form of
`Tr(P Q)/2ⁿ = δ`. With the genuine `2ⁿ×2ⁿ` tensor-product matrices above, the
trace is now constructible, and this section closes the claim end-to-end:
`Tr(Â B̂) = 2ⁿ · ⟪A, B⟫` for the honest operators `Â, B̂` the coefficient maps
denote. That is what licenses `Sum::overlap` being read as a physical
expectation value. (`PPVM.Twisted.twistedConv_apply_id` is the same statement one
level up in the abstraction, purely inside `C[K]`: `⟪A,B⟫` is the identity
coefficient of the L4 twisted product.) -/

/-- **The trace of a tensor-product Pauli is `2ⁿ` on the identity and `0`
otherwise.** Each qubit contributes `tr(g(x,z)) = 2·δ`, and the tensor trace is
the product of the per-qubit traces — so a single non-identity qubit kills it.
This is the tracelessness of `X`, `Y`, `Z` lifted to `n` qubits. -/
theorem trace_tensorPauli (p : Word n) :
    Matrix.trace (tensorPauli p)
      = if p = (fun _ => (false, false)) then (2 : GaussianInt) ^ n else 0 := by
  classical
  have tr1 : ∀ x z : Bool, (∑ j, pauliMat x z j j)
      = if (x, z) = ((false, false) : Bool × Bool) then (2 : GaussianInt) else 0 := by
    intro x z
    rw [Fin.sum_univ_two]
    cases x <;> cases z <;> decide
  have hprod : Matrix.trace (tensorPauli p)
      = ∏ i, (if p i = ((false, false) : Bool × Bool) then (2 : GaussianInt) else 0) := by
    simp only [Matrix.trace, Matrix.diag, tensorPauli]
    rw [← Fintype.prod_sum fun (i : Fin n) (j : Fin 2) => pauliMat (p i).1 (p i).2 j j]
    exact Finset.prod_congr rfl fun i _ => by rw [tr1]
  rw [hprod]
  by_cases hp : p = fun _ => (false, false)
  · subst hp
    simp
  · obtain ⟨i, hi⟩ := Function.ne_iff.mp hp
    rw [if_neg hp, Finset.prod_eq_zero (Finset.mem_univ i) (if_neg hi)]

/-- **Pauli-basis Hilbert–Schmidt orthonormality, on genuine matrices.**
`Tr(g(p) · g(q)) = 2ⁿ · δ_{p,q}` — the real-matrix statement that
`PPVM.Noise.overlap_single_single` was only the model form of. Off the diagonal
`mulWord p q ≠ I` (`mulWord_eq_id_iff`) makes the product traceless; on it the
phase twist vanishes (`phaseExpN_self`), leaving `Tr(I) = 2ⁿ`. -/
theorem trace_tensorPauli_mul (p q : Word n) :
    Matrix.trace (tensorPauli p * tensorPauli q)
      = if q = p then (2 : GaussianInt) ^ n else 0 := by
  classical
  rw [tensorPauli_mul, Matrix.trace_smul, trace_tensorPauli]
  by_cases h : q = p
  · subst h
    rw [if_pos (mulWord_self q), if_pos rfl, phaseExpN_self]
    simp
  · rw [if_neg h, if_neg fun hc => h ((mulWord_eq_id_iff p q).mp hc), smul_zero]

/-- The genuine `2ⁿ×2ⁿ` operator a coefficient map denotes: `Â = ∑ₚ aₚ g(p)`.
This is the semantic function `C[PauliWord] → operators` that the whole `Sum`
abstraction refines. -/
def toOperator (A : Word n → GaussianInt) :
    Matrix (Fin n → Fin 2) (Fin n → Fin 2) GaussianInt :=
  ∑ p, A p • tensorPauli p

/-- **`overlap` is the normalized Hilbert–Schmidt trace pairing.**
`Tr(Â B̂) = 2ⁿ · ∑ₚ aₚ bₚ`, i.e. `∑ₚ aₚ bₚ = Tr(Â B̂)/2ⁿ` on honest `2ⁿ×2ⁿ`
matrices over the exact ring `ℤ[i]`. This removes the design's caveat that the
trace "is not constructed in Lean". -/
theorem trace_toOperator_mul (A B : Word n → GaussianInt) :
    Matrix.trace (toOperator A * toOperator B) = 2 ^ n * ∑ p, A p * B p := by
  classical
  have hterm : ∀ p q : Word n,
      Matrix.trace ((A p • tensorPauli p) * (B q • tensorPauli q))
        = if q = p then A p * B q * 2 ^ n else 0 := by
    intro p q
    rw [Matrix.smul_mul, Matrix.mul_smul, Matrix.trace_smul, Matrix.trace_smul,
      trace_tensorPauli_mul]
    by_cases h : q = p <;> simp [h, mul_assoc]
  simp only [toOperator, Matrix.sum_mul, Matrix.mul_sum, Matrix.trace_sum, hterm]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [Finset.sum_ite_eq Finset.univ p fun q => A q * B p * 2 ^ n, if_pos (Finset.mem_univ p)]
  ring

/-- The same statement against the L3 `Pair::overlap` of `PPVM.GradedMap`:
`⟪A, B⟫ = Tr(Â B̂)/2ⁿ` exactly as `graded.rs` documents it. -/
theorem overlap_eq_trace_div (A B : GradedMap.CMap (Word n) GaussianInt) :
    Matrix.trace (toOperator A * toOperator B) = 2 ^ n * GradedMap.overlap A B := by
  rw [trace_toOperator_mul, GradedMap.overlap_eq_fintype_sum]

end PPVM.PauliMatrix
