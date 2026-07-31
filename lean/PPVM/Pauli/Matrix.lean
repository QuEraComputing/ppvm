/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.NumberTheory.Zsqrtd.GaussianInt
import Mathlib.LinearAlgebra.Matrix.Notation
import PPVM.Pauli.Phase
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

end PPVM.PauliMatrix
