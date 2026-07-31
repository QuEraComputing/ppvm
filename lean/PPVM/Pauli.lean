/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import Mathlib.Data.ZMod.Basic

/-!
# The Pauli group in the symplectic representation

This file starts the formalization of ppvm's Pauli-word algebra. A single-qubit
Pauli operator, taken up to phase, is a symplectic pair of bits
`(x, z) ∈ 𝔽₂ × 𝔽₂`:

| Pauli | `(x, z)` | Rust discriminant (`ppvm-traits/src/char.rs`) |
| :---- | :------: | :-------------------------------------------- |
| `I`   | `(0, 0)` | `0b00` |
| `X`   | `(1, 0)` | `0b01` |
| `Z`   | `(0, 1)` | `0b10` |
| `Y`   | `(1, 1)` | `0b11` |

so the physical operator is `iˣᶻ · Xˣ · Zᶻ` (in particular `Y = i·X·Z`).

What is proved here (all by `decide` over the finite type):

* `mul` is the phase-free group law — symplectic bits add in `𝔽₂` — mirroring
  `crates/ppvm-pauli-word/src/word/mul.rs`. It is commutative, associative,
  unital, and every element is its own inverse.
* `omega` is the symplectic form; `omega p q = 1` iff `p` and `q` anticommute,
  mirroring `PauliWordTrait::anticommutes_at`
  (`crates/ppvm-traits/src/traits/word_trait.rs:113`).

The phase cocycle (the `sign`/`imag` boolean formulas in
`crates/ppvm-pauli-word/src/phase/mul.rs`) is the next target and is described
in `lean/README.md`; it refines this phase-free law to the full Pauli group over
`ℤ/4ℤ`.
-/

namespace PPVM

/-- A single-qubit Pauli operator up to phase, as a symplectic pair of bits
`(x, z) ∈ 𝔽₂ × 𝔽₂`. -/
abbrev Pauli := ZMod 2 × ZMod 2

namespace Pauli

/-- The identity `I = (0, 0)`. -/
def I : Pauli := (0, 0)

/-- `X = (1, 0)`. -/
def X : Pauli := (1, 0)

/-- `Z = (0, 1)`. -/
def Z : Pauli := (0, 1)

/-- `Y = (1, 1)` (recall `Y = i·X·Z`, so its symplectic slot carries both bits). -/
def Y : Pauli := (1, 1)

/-- Phase-free product of two Paulis: symplectic bits add in `𝔽₂`. This is the
group law of `Pₙ / phases ≅ (𝔽₂², ⊕)` and mirrors the `MulAssign` in
`crates/ppvm-pauli-word/src/word/mul.rs`. -/
def mul (p q : Pauli) : Pauli := (p.1 + q.1, p.2 + q.2)

/-- The symplectic form on Paulis. `omega p q = 1` exactly when `p` and `q`
anticommute; the phase law `P·Q = (-1)^{omega p q} · Q·P` is proved (with the
phase tracked) as `PPVM.PauliPhase.phaseExp_sub_comm`. Mirrors `anticommutes_at`
(`crates/ppvm-traits/src/traits/word_trait.rs:113`). -/
def omega (p q : Pauli) : ZMod 2 := p.1 * q.2 + p.2 * q.1

/-- The support/weight of a single-qubit Pauli: `0` for `I`, `1` otherwise. The
n-qubit `weight()` (`crates/ppvm-pauli-word/src/word/data.rs:150`) is the sum of
these over all slots. -/
def weight (p : Pauli) : ℕ := if p = I then 0 else 1

/-! ### The group law (phase-free) -/

theorem mul_comm : ∀ p q : Pauli, mul p q = mul q p := by decide

theorem mul_assoc : ∀ p q r : Pauli, mul (mul p q) r = mul p (mul q r) := by decide

theorem mul_I : ∀ p : Pauli, mul p I = p := by decide

theorem I_mul : ∀ p : Pauli, mul I p = p := by decide

/-- Every Pauli is an involution up to phase: `P·P = I` at the bit level. -/
theorem mul_self : ∀ p : Pauli, mul p p = I := by decide

/-! ### The single-qubit product table (bit level, phase not yet tracked) -/

theorem X_mul_Z : mul X Z = Y := by decide
theorem Z_mul_X : mul Z X = Y := by decide
theorem Y_mul_Z : mul Y Z = X := by decide
theorem X_mul_Y : mul X Y = Z := by decide

/-! ### The symplectic form and (anti)commutation -/

theorem omega_self : ∀ p : Pauli, omega p p = 0 := by decide

theorem omega_comm : ∀ p q : Pauli, omega p q = omega q p := by decide

/-- `X` and `Z` anticommute. -/
theorem omega_X_Z : omega X Z = 1 := by decide

/-- **Distinct nonidentity Paulis on the same qubit anticommute.** (The
complementary cases all commute: `p = q` is `omega_self`; a factor of `I`
commutes by direct computation.) This is the substantive half: `omega p q = 1`
whenever `p ≠ q` and neither is `I`. -/
theorem omega_of_ne : ∀ p q : Pauli, p ≠ q → p ≠ I → q ≠ I → omega p q = 1 := by
  decide

end Pauli

end PPVM
