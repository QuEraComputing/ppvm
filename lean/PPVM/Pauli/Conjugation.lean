/-
Copyright (c) 2026 The PPVM Authors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: The PPVM Authors
-/
import PPVM.Pauli.Phase

/-!
# Clifford conjugation as signed symplectic automorphisms

The design's role-independent `Sp` part
(`traits-2-configuration-and-hashing.md`, "Role-independent"): conjugating a
Pauli by a Clifford gate `G` does the same bit-plane algebra as `G·P·G†`, which
for each Clifford `G` is a **signed symplectic automorphism**

  `G · P · G† = (−1)^{β(G,P)} · φ_G(P)`,

with `φ_G` a symplectic linear map on the `(x,z)` bits and `β` a sign predicate.
Because `PhasedPauli` (`𝒫₁`) is a genuine `Group`, this is literally a group
homomorphism of `𝒫₁`, and we exhibit it for the single-qubit generators:

* `conjH` — `H`-conjugation: the symplectic swap `(x,z) ↦ (z,x)` with sign
  `β = x∧z` (so `HXH=Z`, `HZH=X`, `HYH=−Y`); an *involutive* automorphism
  (`H² = I`);
* `conjS` — `S`-conjugation: the transvection `(x,z) ↦ (x, x⊕z)` with the same
  sign (so `SXS†=Y`, `SZS†=Z`, `SYS†=−X`); an automorphism of order 4.

Each is a `MonoidHom` (`map_mul` by `decide`) and injective, hence an
automorphism — the design's "signed symplectic automorphism," machine-checked.
-/

namespace PPVM.PauliPhase.PhasedPauli

/-- `H`-conjugation on `𝒫₁`: swap the `(x,z)` bits, and add the sign `2·(x∧z)`
(base-`i` exponent), i.e. `HXH=Z`, `HZH=X`, `HYH=−Y`. -/
def conjH (p : PhasedPauli) : PhasedPauli :=
  ⟨p.phase + (if p.x && p.z then 2 else 0), p.z, p.x⟩

/-- `S`-conjugation on `𝒫₁`: `(x,z) ↦ (x, x⊕z)`, sign `2·(x∧z)`, i.e. `SXS†=Y`,
`SZS†=Z`, `SYS†=−X`. -/
def conjS (p : PhasedPauli) : PhasedPauli :=
  ⟨p.phase + (if p.x && p.z then 2 else 0), p.x, xor p.x p.z⟩

-- Raw (over the underlying `mul`) so `decide` never touches the bundled `*`.
private theorem conjH_mul_raw : ∀ p q, conjH (mul p q) = mul (conjH p) (conjH q) := by decide
private theorem conjS_mul_raw : ∀ p q, conjS (mul p q) = mul (conjS p) (conjS q) := by decide

/-- **`H`-conjugation is a group homomorphism of `𝒫₁`.** -/
def conjHHom : PhasedPauli →* PhasedPauli where
  toFun := conjH
  map_one' := by decide
  map_mul' := conjH_mul_raw

/-- **`S`-conjugation is a group homomorphism of `𝒫₁`.** -/
def conjSHom : PhasedPauli →* PhasedPauli where
  toFun := conjS
  map_one' := by decide
  map_mul' := conjS_mul_raw

/-- `conjH` is injective, hence (on the finite group) an automorphism. -/
theorem conjH_injective : Function.Injective conjH := by decide

/-- `conjS` is injective, hence an automorphism. -/
theorem conjS_injective : Function.Injective conjS := by decide

/-! ### The signed-symplectic decomposition (the `Sp` part + the sign `β`) -/

/-- `conjH` acts on the symplectic bits by the swap `φ_H`. -/
theorem conjH_bits (p : PhasedPauli) : (conjH p).x = p.z ∧ (conjH p).z = p.x :=
  ⟨rfl, rfl⟩

/-- The `H`-conjugation sign `β(H,P) = x∧z`: only `Y ↦ −Y` picks up a sign. -/
theorem conjH_sign (p : PhasedPauli) :
    (conjH p).phase = p.phase + (if p.x && p.z then 2 else 0) := rfl

/-- `conjS` acts on the symplectic bits by the transvection `φ_S : (x,z) ↦ (x, x⊕z)`. -/
theorem conjS_bits (p : PhasedPauli) : (conjS p).x = p.x ∧ (conjS p).z = xor p.x p.z :=
  ⟨rfl, rfl⟩

/-- The `S`-conjugation sign `β(S,P) = x∧z`. -/
theorem conjS_sign (p : PhasedPauli) :
    (conjS p).phase = p.phase + (if p.x && p.z then 2 else 0) := rfl

/-! ### The gate laws `H² = I`, and the basis conjugations -/

/-- `H² = I`: `H`-conjugation is an involution. -/
theorem conjH_involutive : Function.Involutive conjH := by
  intro p; revert p; decide

/-- `HXH = Z`, `HZH = X`, `HYH = −Y`, `HIH = I` (phase-explicit). -/
theorem conjH_X : conjH ⟨0, true, false⟩ = ⟨0, false, true⟩ := by decide
theorem conjH_Z : conjH ⟨0, false, true⟩ = ⟨0, true, false⟩ := by decide
theorem conjH_Y : conjH ⟨0, true, true⟩ = ⟨2, true, true⟩ := by decide

/-- `SXS† = Y`, `SZS† = Z`, `SYS† = −X`. -/
theorem conjS_X : conjS ⟨0, true, false⟩ = ⟨0, true, true⟩ := by decide
theorem conjS_Z : conjS ⟨0, false, true⟩ = ⟨0, false, true⟩ := by decide
theorem conjS_Y : conjS ⟨0, true, true⟩ = ⟨2, true, false⟩ := by decide

end PPVM.PauliPhase.PhasedPauli
