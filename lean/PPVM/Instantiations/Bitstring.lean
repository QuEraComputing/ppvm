/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Data.Finsupp.Basic
import Mathlib.Data.ZMod.Basic

/-!
# Third instantiation: the generalized-tableau bitstring algebra

The design's third `Sum` key type
(`traits-2-configuration-and-hashing.md#a-third-instantiation-the-generalized-tableau`):
a `GeneralizedTableau` is `U|c⟩` — a Clifford frame `U` times a sparse
superposition `|c⟩ = ∑_b c_b |b⟩`, whose amplitude vector is `C[Bitstring]`, the
*same* graded free module `C[K]` as `PauliSum`, keyed on `Bitstring`.

Two design claims are validated here:

* **Clifford gate → the frame only.** `G · U|c⟩ = (GU)|c⟩`, so a Clifford leaves
  the amplitude vector `|c⟩` untouched — the identity on `C[Bitstring]`.
* **The XOR relabel used in non-Clifford branching is a bijection.** A branch
  `c_b|b⟩ ↦ cos·c_b|b⟩ + sin·c_b|b⊕s⟩` relabels keys by `b ↦ b ⊕ s`; that this
  is a bijection is Yoder-2012 Lemma 5 / SOFT Eq. 4's requirement
  (`data.rs:803`), and it is what makes the relabel lift to an isomorphism of
  the amplitude module.
-/

namespace PPVM.GenTableau

variable {n : ℕ} {C : Type*}

/-- A computational-basis bitstring on `n` qubits: `Fin n → 𝔽₂`. XOR of
bitstrings is pointwise addition. -/
abbrev Bitstring (n : ℕ) := Fin n → ZMod 2

/-- The amplitude vector of a `GeneralizedTableau`: `C[Bitstring]`, the same free
`C`-module `C[K]` as `PauliSum`, with `K = Bitstring n`. -/
abbrev Amplitudes (n : ℕ) (C : Type*) [Zero C] := Bitstring n →₀ C

/-! ### The XOR relabel is a bijection -/

theorem zmod2_add_self : ∀ x : ZMod 2, x + x = 0 := by decide

/-- `s + s = 0` for bitstrings (each qubit is `𝔽₂`). -/
theorem add_self (s : Bitstring n) : s + s = 0 := by
  funext i; simp only [Pi.add_apply, Pi.zero_apply]; exact zmod2_add_self (s i)

/-- **The XOR relabel `b ↦ b ⊕ s` is a bijection** of bitstrings, presented as an
involutive `Equiv`. This is the branching relabel the generalized tableau uses;
being an `Equiv` is exactly "it is a bijection." -/
def xorRelabel (s : Bitstring n) : Bitstring n ≃ Bitstring n where
  toFun b := b + s
  invFun b := b + s
  left_inv b := by simp [add_assoc, add_self]
  right_inv b := by simp [add_assoc, add_self]

@[simp] theorem xorRelabel_apply (s b : Bitstring n) : xorRelabel s b = b + s := rfl

/-- The relabel is its own inverse (`⊕ s` twice is identity). -/
theorem xorRelabel_involutive (s : Bitstring n) :
    Function.Involutive (xorRelabel s) := by
  intro b; simp [add_assoc, add_self]

/-- The relabel is genuinely a bijection (the design's stated requirement). -/
theorem xorRelabel_bijective (s : Bitstring n) :
    Function.Bijective (xorRelabel s) := (xorRelabel s).bijective

/-! ### Lifting the relabel to the amplitude module

A bijection on keys lifts to an isomorphism of `C[Bitstring]` — `Finsupp`'s
`equivMapDomain`. So the branch's key relabel is a linear reindexing of the
amplitude vector, never a lossy operation; in particular it preserves the number
of stored terms. -/

/-- Reducing `(xorRelabel s).symm` — it is again `· + s` (the relabel is
involutive), which the amplitude proofs below use definitionally. -/
theorem xorRelabel_symm_apply (s b : Bitstring n) :
    (xorRelabel s).symm b = b + s := rfl

/-- The branch relabel lifted to amplitudes: reindex `|b⟩ ↦ |b ⊕ s⟩` via
`Finsupp.equivMapDomain`, packaged as an involutive `Equiv` of `C[Bitstring]`. -/
def relabelAmp [Zero C] (s : Bitstring n) : Amplitudes n C ≃ Amplitudes n C where
  toFun := Finsupp.equivMapDomain (xorRelabel s)
  invFun := Finsupp.equivMapDomain (xorRelabel s)
  left_inv c := by
    ext b
    simp only [Finsupp.equivMapDomain_apply, xorRelabel_symm_apply]
    rw [add_assoc, add_self, add_zero]
  right_inv c := by
    ext b
    simp only [Finsupp.equivMapDomain_apply, xorRelabel_symm_apply]
    rw [add_assoc, add_self, add_zero]

/-- Relabeling reads off the amplitude at the XOR-shifted key. -/
@[simp] theorem relabelAmp_apply [Zero C] (s : Bitstring n) (c : Amplitudes n C)
    (b : Bitstring n) : relabelAmp s c b = c (b + s) := by
  simp only [relabelAmp, Equiv.coe_fn_mk, Finsupp.equivMapDomain_apply,
    xorRelabel_symm_apply]

/-- The lifted relabel is a bijection of the amplitude module (it is an
`Equiv`), so branching only *moves* amplitude weight, never loses it. -/
theorem relabelAmp_bijective [Zero C] (s : Bitstring n) :
    Function.Bijective (relabelAmp (C := C) s) := (relabelAmp (C := C) s).bijective

/-! ### A Clifford gate touches the frame only

A `GeneralizedTableau` state is `U|c⟩` — a Clifford frame `U` paired with the
amplitude vector `|c⟩`. We model that as a pair `GenState = F × Amplitudes` (with
`F` the frame, kept abstract), and a Clifford gate as the map that updates the
frame factor and leaves the amplitude factor alone — the design's
`G · U|c⟩ = (GU)|c⟩`. The two projections below make that factorization explicit:
the gate is `frameUpdate` on the first component and `id` on the second. This is
a *modeling* statement (it holds by how the state factorizes), not an emergent
theorem — its value is pinning down that amplitudes and frame are independent
factors and the Clifford acts only on the frame. -/

/-- A generalized-tableau state: a Clifford frame (any type `F`) and the amplitude
vector. -/
abbrev GenState (F : Type*) (n : ℕ) (C : Type*) [Zero C] := F × Amplitudes n C

/-- A Clifford gate acts on the frame factor via `frameUpdate`, leaving amplitudes
untouched. -/
def cliffordStep {F : Type*} [Zero C] (frameUpdate : F → F) (s : GenState F n C) :
    GenState F n C := (frameUpdate s.1, s.2)

/-- **A Clifford gate leaves the amplitude vector fixed** (`G·U|c⟩ = (GU)|c⟩`):
its action projects to the identity on the amplitude factor. -/
theorem cliffordStep_amplitudes {F : Type*} [Zero C] (frameUpdate : F → F)
    (s : GenState F n C) : (cliffordStep frameUpdate s).2 = s.2 := rfl

/-- …and it updates the frame factor by exactly the given frame map. -/
theorem cliffordStep_frame {F : Type*} [Zero C] (frameUpdate : F → F)
    (s : GenState F n C) : (cliffordStep frameUpdate s).1 = frameUpdate s.1 := rfl

end PPVM.GenTableau
