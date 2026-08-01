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

The phase-tracking simulator runs the **backward** direction `S†·P·S` (its `s`
gate), whose sign is the sole convention-sensitive Clifford sign; `conjSdag`
(sign `x∧¬z`, so `S†XS=−Y`) pins the exact `ℤ/4ℤ` delta `PhaseTrack::s_phase`
emits, and `conjS_conjSdag` shows it is the inverse of `conjS`.

Each is a `MonoidHom` (`map_mul` by `decide`) and injective, hence an
automorphism — the design's "signed symplectic automorphism," machine-checked.

The two-qubit generators `CNOT`/`CZ` are handled in the second half of the file
(`TwoPauli`, the group `𝒫₂`): they are independent generators whose conjugation
*phase* is not derivable from the single-qubit `conjH_sign`/`conjS_sign` and is
not covered by the bit-only `cnotAct_isometry`/`czAct_isometry` of
`Symplectic.lean`. `conjCNOT_sign`/`conjCZ_sign` pin the exact `ℤ/4ℤ` phase delta
the Rust kernel commits to (`crates/ppvm-pauli-word/src/phase/clifford.rs`).
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

/-- **`S`-conjugation has order 4** (`S⁴` is a global phase, so conjugation by it
is trivial): `conjS⁴ = id` while `conjS² ≠ id` (conjugation by `S² = Z`). -/
theorem conjS_iterate_four : conjS^[4] = id := by
  funext p; revert p; decide

theorem conjS_iterate_two_ne_id : conjS^[2] ≠ id := by decide

/-! ### `S†`-conjugation: the backward (Heisenberg) direction the simulator runs

`conjS` above is the **forward** `S·P·S†`. The phase-tracking simulator conjugates
in the **backward** direction `S†·P·S` (its `s` gate advances an operator by `S†`
in the Heisenberg picture). On the bits this is the *same* transvection
`(x,z) ↦ (x, x⊕z)` (an involution), but the sign fires on `x∧¬z` instead of `x∧z`.
`S` is the sole Clifford generator whose conjugation sign is convention-sensitive
(the self-adjoint `H`/`X`/`Y`/`Z`/`CNOT`/`CZ` signs coincide forward and backward),
so the backward sign the phased word emits needs its own witness. This is exactly
what `PhaseTrack::s_phase` commits to
(`crates/ppvm-phased-pauli-word-2/src/clifford.rs:104`). -/

/-- `S†`-conjugation on `𝒫₁`: `(x,z) ↦ (x, x⊕z)` with sign `2·(x∧¬z)`, i.e.
`S†XS = −Y`, `S†YS = X`, `S†ZS = Z`. -/
def conjSdag (p : PhasedPauli) : PhasedPauli :=
  ⟨p.phase + (if p.x && !p.z then 2 else 0), p.x, xor p.x p.z⟩

private theorem conjSdag_mul_raw :
    ∀ p q, conjSdag (mul p q) = mul (conjSdag p) (conjSdag q) := by decide

/-- **`S†`-conjugation is a group homomorphism of `𝒫₁`.** -/
def conjSdagHom : PhasedPauli →* PhasedPauli where
  toFun := conjSdag
  map_one' := by decide
  map_mul' := conjSdag_mul_raw

/-- `conjSdag` acts on the bits by the same transvection `φ_S` as `conjS`. -/
theorem conjSdag_bits (p : PhasedPauli) :
    (conjSdag p).x = p.x ∧ (conjSdag p).z = xor p.x p.z := ⟨rfl, rfl⟩

/-- **The `S†`-conjugation sign** `β(S†,P) = x∧¬z` — the exact `ℤ/4ℤ` delta
`PhaseTrack::s_phase` applies (`crates/ppvm-phased-pauli-word-2/src/clifford.rs:104`).
It differs from the forward `conjS_sign` (`x∧z`) precisely on `X` and `Y`, the
convention-sensitive cases. -/
theorem conjSdag_sign (p : PhasedPauli) :
    (conjSdag p).phase = p.phase + (if p.x && !p.z then 2 else 0) := rfl

/-- `S†XS = −Y`, `S†YS = X`, `S†ZS = Z` (phase-explicit) — the simulator's `s`
table (`clifford.rs` `single_qubit_gates_track_sign`: `+X ↦ −Y`, `+Y ↦ +X`,
`+Z ↦ +Z`). -/
theorem conjSdag_X : conjSdag ⟨0, true, false⟩ = ⟨2, true, true⟩ := by decide
theorem conjSdag_Y : conjSdag ⟨0, true, true⟩ = ⟨0, true, false⟩ := by decide
theorem conjSdag_Z : conjSdag ⟨0, false, true⟩ = ⟨0, false, true⟩ := by decide

/-- `conjSdag` is the two-sided **inverse** of `conjS` (conjugating by `S` then by
`S†`, i.e. by `S†·S = I`, is the identity). This pins `conjSdag` to genuine
`S†`-conjugation rather than an arbitrary sign choice, and shows the two signs
`x∧z` and `x∧¬z` are the forward/backward pair. -/
theorem conjS_conjSdag : Function.LeftInverse conjS conjSdag := by
  intro p; revert p; decide
theorem conjSdag_conjS : Function.LeftInverse conjSdag conjS := by
  intro p; revert p; decide

/-! ### Pure-sign conjugation by the Pauli generators `X`, `Y`, `Z`

Unlike `H`/`S` (which move the word), conjugating by a Pauli generator `G ∈
{X,Y,Z}` **fixes the word** and applies only a `±1 = (−1)^{ω(G,P)}` sign. Because
each generator is Hermitian and squares to `+I`, `G⁻¹ = G` and `G·P·G = G·P·G⁻¹`
is literal group conjugation in `𝒫₁`; expanding it (`mul (mul G P) G`, so `decide`
never touches the bundled `*`) gives `word(P)` back with the phase advanced by `2`
exactly on the anticommuting Paulis. These are the exact signs the `PauliSum`
`x`/`y`/`z` **fast path** folds straight into the coefficient
(`ppvm-pauli-sum-2/src/clifford.rs`, `Sum::x`/`y`/`z` via `flip_sign_by_key`),
*without* routing through the audited phased-word `Clifford` impl — so the
duplicated fast-path signs are pinned here directly against the group product. -/

/-- **`X·P·X = (−1)^{z} P`.** The word is fixed; the phase advances by `2` iff
`p.z`. Matches `Sum::x` (`flip iff z_bit`, `clifford.rs:143`). -/
theorem conjX (p : PhasedPauli) :
    mul (mul ⟨0, true, false⟩ p) ⟨0, true, false⟩
      = ⟨p.phase + (if p.z then 2 else 0), p.x, p.z⟩ := by
  revert p; decide

/-- **`Y·P·Y = (−1)^{x⊕z} P`.** Word fixed; phase advances by `2` iff `p.x ⊕ p.z`.
Matches `Sum::y` (`flip iff x_bit ⊕ z_bit`, `clifford.rs:158`). -/
theorem conjY (p : PhasedPauli) :
    mul (mul ⟨0, true, true⟩ p) ⟨0, true, true⟩
      = ⟨p.phase + (if xor p.x p.z then 2 else 0), p.x, p.z⟩ := by
  revert p; decide

/-- **`Z·P·Z = (−1)^{x} P`.** Word fixed; phase advances by `2` iff `p.x`. Matches
`Sum::z` (`flip iff x_bit`, `clifford.rs:173`). -/
theorem conjZ (p : PhasedPauli) :
    mul (mul ⟨0, false, true⟩ p) ⟨0, false, true⟩
      = ⟨p.phase + (if p.x then 2 else 0), p.x, p.z⟩ := by
  revert p; decide

/-! ### A Clifford conjugation never emits `±i`: the phase stays real

`ppvm-pauli-sum-2/src/clifford.rs` conjugates each Pauli on a `Phased<PauliWord>`
started at phase `+1`, then drains the resulting phase to a coefficient `±1`
(`clifford_sign`), treating a `±i` result as a bug. That drain is total precisely
because each single-qubit generator's conjugation *delta* is the real
`if … then 2 else 0` (`conjH_sign`/`conjS_sign`/`conjSdag_sign`), so a real input
phase stays real (`IsRealPhase` is closed under adding such a delta). Starting
from `+1` (`isRealPhase_zero`), the phase is therefore `Pos1`/`Neg1` at every
step — the `PosI`/`NegI` branch of `clifford_sign` is unreachable. -/

/-- **`H`-conjugation keeps the phase real** (`±1`, never `±i`). -/
theorem conjH_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjH p).phase := by
  rw [conjH_sign]; exact hp.add (isRealPhase_ite_two _)

/-- **`S`-conjugation keeps the phase real.** -/
theorem conjS_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjS p).phase := by
  rw [conjS_sign]; exact hp.add (isRealPhase_ite_two _)

/-- **`S†`-conjugation (the backward direction the simulator runs) keeps the phase
real.** -/
theorem conjSdag_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjSdag p).phase := by
  rw [conjSdag_sign]; exact hp.add (isRealPhase_ite_two _)

/-! ### The derived Clifford extensions as products of the audited generators

`crates/ppvm-traits-2/src/pauli.rs` implements `CliffordExtensions` **once**, for
every `BlanketClifford` opt-in (including the future phase-carrying `Tableau`),
by writing each derived gate as a *product of the audited generators* instead of
a fresh hand-written bit-and-sign rule. Its justification — "the bit rule and the
sign then follow from generators whose signs are already machine-checked" — is
only as good as the composition, so the composition is discharged here.

The crate conjugates in the **backward** Heisenberg direction `P ↦ U†PU`, so its
`s` gate is `conjSdag` (not `conjS`), while its self-adjoint `h`/`z` gates are
`conjH`/`conjPauliZ` either way. Applying `A` then `B` therefore conjugates by
the operator product, giving the five single-qubit rows of the crate's table:

| crate gate | operator identity | call sequence | composite below |
|:---|:---|:---|:---|
| `s_dag` | `S† = S·Z` | `s`, `z` | `extSdag` |
| `sqrt_x` | `√X ≃ H·S·H` | `h`, `s`, `h` | `extSqrtX` |
| `sqrt_x_dag` | `√X† ≃ H·S†·H` | `h`, `s_dag`, `h` | `extSqrtXdag` |
| `sqrt_y` | `√Y ≃ H·Z` | `h`, `z` | `extSqrtY` |
| `sqrt_y_dag` | `√Y† ≃ Z·H` | `z`, `h` | `extSqrtYdag` |

Each composite is a `MonoidHom` *for free* (a composite of the audited
generator homs), and the theorems below check that it reproduces the standard
conjugation table the old crate hand-wrote
(`crates/ppvm-traits/src/traits/clifford.rs`, replayed on a ℤ₄-phased stub by
`crates/ppvm-traits-2/tests/phase1_gate_surface.rs`
`blanket_clifford_extensions_match_old_conjugation_table`). Corollary at the end:
every composite delta is still real, so `ppvm-pauli-sum-2`'s `±1` drain stays
total on the extension gates too. -/

/-- `Z`-conjugation packaged as a *map* on `𝒫₁` (the crate's `z` gate): the word
is fixed and the phase advances by `2·x`. This is the closed form `conjZ` proves,
named so it can appear inside a composite. -/
def conjPauliZ (p : PhasedPauli) : PhasedPauli :=
  ⟨p.phase + (if p.x then 2 else 0), p.x, p.z⟩

/-- `conjPauliZ` really is conjugation by `Z` in the group `𝒫₁` (`conjZ`). -/
theorem conjPauliZ_eq_conj (p : PhasedPauli) :
    conjPauliZ p = mul (mul ⟨0, false, true⟩ p) ⟨0, false, true⟩ := (conjZ p).symm

private theorem conjPauliZ_mul_raw :
    ∀ p q, conjPauliZ (mul p q) = mul (conjPauliZ p) (conjPauliZ q) := by decide

/-- **`Z`-conjugation is a group homomorphism of `𝒫₁`** (conjugation by a fixed
group element always is). -/
def conjPauliZHom : PhasedPauli →* PhasedPauli where
  toFun := conjPauliZ
  map_one' := by decide
  map_mul' := conjPauliZ_mul_raw

/-- **`Z`-conjugation keeps the phase real.** -/
theorem conjPauliZ_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjPauliZ p).phase :=
  hp.add (isRealPhase_ite_two _)

/-- The crate's `s_dag` gate: `s` then `z`, i.e. conjugation by `S·Z = S†` in the
backward convention. -/
def extSdag (p : PhasedPauli) : PhasedPauli := conjPauliZ (conjSdag p)

/-- The crate's `sqrt_x` gate: `h`, `s`, `h`. -/
def extSqrtX (p : PhasedPauli) : PhasedPauli := conjH (conjSdag (conjH p))

/-- The crate's `sqrt_x_dag` gate: `h`, `s_dag`, `h`. -/
def extSqrtXdag (p : PhasedPauli) : PhasedPauli := conjH (extSdag (conjH p))

/-- The crate's `sqrt_y` gate: `h`, `z`. -/
def extSqrtY (p : PhasedPauli) : PhasedPauli := conjPauliZ (conjH p)

/-- The crate's `sqrt_y_dag` gate: `z`, `h`. -/
def extSqrtYdag (p : PhasedPauli) : PhasedPauli := conjH (conjPauliZ p)

/-! Each extension gate is a `MonoidHom` **by construction** — the composite of
the audited generator homs — so no new sign algebra is introduced. -/

/-- `s_dag` as a generator product. -/
def extSdagHom : PhasedPauli →* PhasedPauli := conjPauliZHom.comp conjSdagHom
/-- `sqrt_x` as a generator product. -/
def extSqrtXHom : PhasedPauli →* PhasedPauli :=
  conjHHom.comp (conjSdagHom.comp conjHHom)
/-- `sqrt_x_dag` as a generator product. -/
def extSqrtXdagHom : PhasedPauli →* PhasedPauli :=
  conjHHom.comp (extSdagHom.comp conjHHom)
/-- `sqrt_y` as a generator product. -/
def extSqrtYHom : PhasedPauli →* PhasedPauli := conjPauliZHom.comp conjHHom
/-- `sqrt_y_dag` as a generator product. -/
def extSqrtYdagHom : PhasedPauli →* PhasedPauli := conjHHom.comp conjPauliZHom

@[simp] theorem extSdagHom_apply (p : PhasedPauli) : extSdagHom p = extSdag p := rfl
@[simp] theorem extSqrtXHom_apply (p : PhasedPauli) : extSqrtXHom p = extSqrtX p := rfl
@[simp] theorem extSqrtXdagHom_apply (p : PhasedPauli) :
    extSqrtXdagHom p = extSqrtXdag p := rfl
@[simp] theorem extSqrtYHom_apply (p : PhasedPauli) : extSqrtYHom p = extSqrtY p := rfl
@[simp] theorem extSqrtYdagHom_apply (p : PhasedPauli) :
    extSqrtYdagHom p = extSqrtYdag p := rfl

/-- **The crate's `s_dag` gate is exactly the forward `S`-conjugation.** The
product `s ; z` (backward `S†` then `Z`) collapses to `conjS`, so its table is
`conjS_X`/`conjS_Y`/`conjS_Z`: `X ↦ Y`, `Y ↦ −X`, `Z ↦ Z`. -/
theorem extSdag_eq_conjS : extSdag = conjS := by
  funext p; revert p; decide

/-- **`sqrt_x` table**: `X ↦ X`, `Y ↦ −Z`, `Z ↦ Y`. -/
theorem extSqrtX_X : extSqrtX ⟨0, true, false⟩ = ⟨0, true, false⟩ := by decide
theorem extSqrtX_Y : extSqrtX ⟨0, true, true⟩ = ⟨2, false, true⟩ := by decide
theorem extSqrtX_Z : extSqrtX ⟨0, false, true⟩ = ⟨0, true, true⟩ := by decide

/-- **`sqrt_x_dag` table**: `X ↦ X`, `Y ↦ Z`, `Z ↦ −Y`. -/
theorem extSqrtXdag_X : extSqrtXdag ⟨0, true, false⟩ = ⟨0, true, false⟩ := by decide
theorem extSqrtXdag_Y : extSqrtXdag ⟨0, true, true⟩ = ⟨0, false, true⟩ := by decide
theorem extSqrtXdag_Z : extSqrtXdag ⟨0, false, true⟩ = ⟨2, true, true⟩ := by decide

/-- **`sqrt_y` table**: `X ↦ Z`, `Y ↦ Y`, `Z ↦ −X`. -/
theorem extSqrtY_X : extSqrtY ⟨0, true, false⟩ = ⟨0, false, true⟩ := by decide
theorem extSqrtY_Y : extSqrtY ⟨0, true, true⟩ = ⟨0, true, true⟩ := by decide
theorem extSqrtY_Z : extSqrtY ⟨0, false, true⟩ = ⟨2, true, false⟩ := by decide

/-- **`sqrt_y_dag` table**: `X ↦ −Z`, `Y ↦ Y`, `Z ↦ X`. -/
theorem extSqrtYdag_X : extSqrtYdag ⟨0, true, false⟩ = ⟨2, false, true⟩ := by decide
theorem extSqrtYdag_Y : extSqrtYdag ⟨0, true, true⟩ = ⟨0, true, true⟩ := by decide
theorem extSqrtYdag_Z : extSqrtYdag ⟨0, false, true⟩ = ⟨0, true, false⟩ := by decide

/-- **Each dagger really is the inverse conjugation** — the property the `S³ = S†`
/ `H·S·H` products must not silently break
(`phase1_gate_surface.rs::blanket_extension_daggers_are_inverse_conjugations`). -/
theorem extSqrtXdag_extSqrtX : Function.LeftInverse extSqrtXdag extSqrtX := by
  intro p; revert p; decide
theorem extSqrtX_extSqrtXdag : Function.LeftInverse extSqrtX extSqrtXdag := by
  intro p; revert p; decide
theorem extSqrtYdag_extSqrtY : Function.LeftInverse extSqrtYdag extSqrtY := by
  intro p; revert p; decide
theorem extSqrtY_extSqrtYdag : Function.LeftInverse extSqrtY extSqrtYdag := by
  intro p; revert p; decide
theorem extSdag_conjSdag : Function.LeftInverse extSdag conjSdag := by
  intro p; revert p; decide

/-- **`√X² = X`**: applying `sqrt_x` twice is `X`-conjugation, the pure sign
`(−1)^z` of `conjX` (same test). -/
theorem extSqrtX_sq (p : PhasedPauli) :
    extSqrtX (extSqrtX p) = ⟨p.phase + (if p.z then 2 else 0), p.x, p.z⟩ := by
  revert p; decide

/-- **`√Y² = Y`**: applying `sqrt_y` twice is `Y`-conjugation, the sign
`(−1)^{x⊕z}` of `conjY`. -/
theorem extSqrtY_sq (p : PhasedPauli) :
    extSqrtY (extSqrtY p) = ⟨p.phase + (if xor p.x p.z then 2 else 0), p.x, p.z⟩ := by
  revert p; decide

/-! **Corollary — the composite phases stay real.** Each extension gate is a
composite of generators whose deltas are the real `if … then 2 else 0`, so it
maps a real phase to a real phase; the `±1` drain of
`ppvm-pauli-sum-2/src/clifford.rs` is therefore total on the extension gates too,
not just on `H`/`S`/`CNOT`/`CZ`. -/

theorem extSdag_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (extSdag p).phase :=
  conjPauliZ_isRealPhase (conjSdag_isRealPhase hp)

theorem extSqrtX_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (extSqrtX p).phase :=
  conjH_isRealPhase (conjSdag_isRealPhase (conjH_isRealPhase hp))

theorem extSqrtXdag_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (extSqrtXdag p).phase :=
  conjH_isRealPhase (extSdag_isRealPhase (conjH_isRealPhase hp))

theorem extSqrtY_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (extSqrtY p).phase :=
  conjPauliZ_isRealPhase (conjH_isRealPhase hp)

theorem extSqrtYdag_isRealPhase {p : PhasedPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (extSqrtYdag p).phase :=
  conjH_isRealPhase (conjPauliZ_isRealPhase hp)

end PPVM.PauliPhase.PhasedPauli

/-! ## Two-qubit conjugation: the `CNOT`/`CZ` phase rule

`H` and `S` are single-qubit generators, so their conjugation is a signed
automorphism of `𝒫₁` (above). `CNOT`/`CZ` are *independent* two-qubit
generators: `Symplectic.lean` proves only the **bit** half (`cnotAct_isometry`/
`czAct_isometry` — that they land in `Sp(2n,2)`), and the single-qubit
`conjH_sign`/`conjS_sign` say nothing about them. The **phase** half — the exact
`ℤ/4ℤ` delta the blanket `impl Clifford` in `ppvm-traits-2/src/pauli.rs`
commits `cnot_phase`/`cz_phase` to (concretely
`crates/ppvm-pauli-word/src/phase/clifford.rs:80,95`) — is pinned here.

We build the two-qubit phased Pauli group `𝒫₂` as the external central product
of two copies of `𝒫₁` (independent qubits, one shared global phase in `ℤ/4ℤ`, so
the two per-qubit cocycles simply add), and exhibit `CNOT`/`CZ` conjugation as
**signed symplectic automorphisms** of it. Because conjugation by any fixed
unitary is a group homomorphism, proving each map is a `MonoidHom` *and* checking
it on the four generators `X_c, Z_c, X_t, Z_t` (the standard tableau table)
pins it to genuine `G·P·G†`; the phase deltas on the remaining Paulis are then
forced, and `conjCNOT_sign`/`conjCZ_sign` read off that the forced delta is
exactly the Rust boolean formula. -/

namespace PPVM.PauliPhase

/-- A two-qubit phased Pauli: a global `ℤ/4ℤ` phase and the `(x,z)` bits of the
control (`c`) and target (`t`) qubits — a point of `ℤ₄ × 𝔽₂⁴`. -/
@[ext]
structure TwoPauli where
  /-- Base-`i` global phase exponent in `ℤ/4ℤ`. -/
  phase : ZMod 4
  /-- Control X bit. -/
  xc : Bool
  /-- Control Z bit. -/
  zc : Bool
  /-- Target X bit. -/
  xt : Bool
  /-- Target Z bit. -/
  zt : Bool
deriving DecidableEq, Fintype

namespace TwoPauli

/-- The two-qubit twisted product: per-qubit bits `⊕`, phases add, plus **both**
per-qubit cocycles (`𝒫₂` is the external central product of two `𝒫₁`'s). -/
def mul (p q : TwoPauli) : TwoPauli where
  phase := p.phase + q.phase
    + phaseExp p.xc p.zc q.xc q.zc + phaseExp p.xt p.zt q.xt q.zt
  xc := xor p.xc q.xc
  zc := xor p.zc q.zc
  xt := xor p.xt q.xt
  zt := xor p.zt q.zt

/-- The identity `+1 · I⊗I`. -/
def one : TwoPauli := ⟨0, false, false, false, false⟩

/-- Inverse: same bits, negated phase (each Pauli squares to `+I`). -/
def inv (p : TwoPauli) : TwoPauli := ⟨-p.phase, p.xc, p.zc, p.xt, p.zt⟩

theorem mul_assoc' (p q r : TwoPauli) : mul (mul p q) r = mul p (mul q r) := by
  have hc := phaseExp_cocycle p.xc p.zc q.xc q.zc r.xc r.zc
  have ht := phaseExp_cocycle p.xt p.zt q.xt q.zt r.xt r.zt
  ext
  · simp only [mul]; linear_combination hc + ht
  · simp only [mul, Bool.xor_assoc]
  · simp only [mul, Bool.xor_assoc]
  · simp only [mul, Bool.xor_assoc]
  · simp only [mul, Bool.xor_assoc]

theorem one_mul' (p : TwoPauli) : mul one p = p := by
  ext <;>
    simp only [mul, one, phaseExp_id_left, Bool.false_xor, zero_add, add_zero]

theorem mul_one' (p : TwoPauli) : mul p one = p := by
  ext <;>
    simp only [mul, one, phaseExp_id_right, Bool.xor_false, add_zero]

theorem inv_mul_cancel' (p : TwoPauli) : mul (inv p) p = one := by
  ext <;>
    simp only [mul, inv, one, phaseExp_self, Bool.xor_self, neg_add_cancel, add_zero]

/-- **The two-qubit phased Pauli group `𝒫₂`.** -/
instance : Group TwoPauli where
  mul := mul
  one := one
  inv := inv
  mul_assoc := mul_assoc'
  one_mul := one_mul'
  mul_one := mul_one'
  inv_mul_cancel := inv_mul_cancel'

@[simp] theorem mul_def (p q : TwoPauli) :
    p * q = ⟨p.phase + q.phase
      + phaseExp p.xc p.zc q.xc q.zc + phaseExp p.xt p.zt q.xt q.zt,
      xor p.xc q.xc, xor p.zc q.zc, xor p.xt q.xt, xor p.zt q.zt⟩ := rfl

/-! ### `CNOT` conjugation

Bit half (`Symplectic.cnotAct`, control `c` target `t`): `x_t ⊕= x_c`,
`z_c ⊕= z_t`. Phase half (`phase/clifford.rs:80`): a `−1 = i²` exactly when
`x_c ∧ z_t ∧ (x_t = z_c)`. -/

/-- The `CNOT`-conjugation phase delta `β(CNOT,P)` as a `ℤ/4ℤ` value
(`phase/clifford.rs:80`): a `−1 = i²` exactly when `x_c ∧ z_t ∧ (x_t = z_c)`.
Kept as a named function so the homomorphism proof treats it as one atom. -/
def cnotDelta (xc zt xt zc : Bool) : ZMod 4 :=
  if xc && zt && (xt == zc) then 2 else 0

/-- `CNOT`-conjugation on `𝒫₂`. -/
def conjCNOT (p : TwoPauli) : TwoPauli where
  phase := p.phase + cnotDelta p.xc p.zt p.xt p.zc
  xc := p.xc
  zc := xor p.zc p.zt
  xt := xor p.xt p.xc
  zt := p.zt

/-! ### `CZ` conjugation

Bit half (`Symplectic.czAct`): `z_c ⊕= x_t`, `z_t ⊕= x_c`. Phase half
(`phase/clifford.rs:95`): a `−1` exactly when `x_c ∧ x_t ∧ (z_c ≠ z_t)`. -/

/-- The `CZ`-conjugation phase delta `β(CZ,P)` as a `ℤ/4ℤ` value
(`phase/clifford.rs:95`): a `−1` exactly when `x_c ∧ x_t ∧ (z_c ≠ z_t)`. -/
def czDelta (xc xt zc zt : Bool) : ZMod 4 :=
  if xc && xt && xor zc zt then 2 else 0

/-- `CZ`-conjugation on `𝒫₂`. -/
def conjCZ (p : TwoPauli) : TwoPauli where
  phase := p.phase + czDelta p.xc p.xt p.zc p.zt
  xc := p.xc
  zc := xor p.zc p.xt
  xt := p.xt
  zt := xor p.zt p.xc

-- The hom property is proved structurally, not by `decide` over the 64-element
-- type: the global phases `φ, ψ` cancel, leaving a pure 8-boolean identity for
-- the phase deltas (`*_phase_delta`, 256 cases) and `xor`-associativity for the
-- bits. A direct `decide` would enumerate all `64² = 4096` operand pairs
-- (including redundant `ℤ/4ℤ` phase values) and blow the heartbeat budget.

/-- The `CNOT` phase-delta 2-cocycle-compatibility, over the 8 operand bits. -/
private theorem conjCNOT_phase_delta : ∀ a b c d e f g h : Bool,
    phaseExp a b e f + phaseExp c d g h
      + cnotDelta (xor a e) (xor d h) (xor c g) (xor b f)
    = cnotDelta a d c b + cnotDelta e h g f
      + phaseExp a (xor b d) e (xor f h) + phaseExp (xor c a) d (xor g e) h := by decide

/-- The `CZ` phase-delta 2-cocycle-compatibility, over the 8 operand bits. -/
private theorem conjCZ_phase_delta : ∀ a b c d e f g h : Bool,
    phaseExp a b e f + phaseExp c d g h
      + czDelta (xor a e) (xor c g) (xor b f) (xor d h)
    = czDelta a c b d + czDelta e g f h
      + phaseExp a (xor b c) e (xor f g) + phaseExp c (xor d a) g (xor h e) := by decide

private theorem conjCNOT_mul_raw :
    ∀ p q, conjCNOT (mul p q) = mul (conjCNOT p) (conjCNOT q) := by
  rintro ⟨φ, a, b, c, d⟩ ⟨ψ, e, f, g, h⟩
  ext
  · simp only [conjCNOT, mul]
    linear_combination conjCNOT_phase_delta a b c d e f g h
  all_goals (simp only [conjCNOT, mul]; try (revert a b c d e f g h; decide))

private theorem conjCZ_mul_raw :
    ∀ p q, conjCZ (mul p q) = mul (conjCZ p) (conjCZ q) := by
  rintro ⟨φ, a, b, c, d⟩ ⟨ψ, e, f, g, h⟩
  ext
  · simp only [conjCZ, mul]
    linear_combination conjCZ_phase_delta a b c d e f g h
  all_goals (simp only [conjCZ, mul]; try (revert a b c d e f g h; decide))

/-- **`CNOT`-conjugation is a group homomorphism of `𝒫₂`** — hence, being
injective, a signed symplectic automorphism. -/
def conjCNOTHom : TwoPauli →* TwoPauli where
  toFun := conjCNOT
  map_one' := by decide
  map_mul' := conjCNOT_mul_raw

/-- **`CZ`-conjugation is a group homomorphism of `𝒫₂`.** -/
def conjCZHom : TwoPauli →* TwoPauli where
  toFun := conjCZ
  map_one' := by decide
  map_mul' := conjCZ_mul_raw

/-- `CNOT` is its own inverse (`CNOT² = I`), so its conjugation is an involution,
hence injective — an automorphism. -/
theorem conjCNOT_involutive : Function.Involutive conjCNOT := by
  intro p; revert p; decide
theorem conjCNOT_injective : Function.Injective conjCNOT := conjCNOT_involutive.injective

/-- `CZ` is its own inverse (`CZ² = I`); its conjugation is an involution. -/
theorem conjCZ_involutive : Function.Involutive conjCZ := by
  intro p; revert p; decide
theorem conjCZ_injective : Function.Injective conjCZ := conjCZ_involutive.injective

/-! ### The signed-symplectic decomposition (the `Sp` part + the sign `β`) -/

/-- `conjCNOT` acts on the bits by `Symplectic.cnotAct`: `x_t ⊕= x_c`,
`z_c ⊕= z_t` (control X and target Z unchanged). -/
theorem conjCNOT_bits (p : TwoPauli) :
    (conjCNOT p).xc = p.xc ∧ (conjCNOT p).zc = xor p.zc p.zt
      ∧ (conjCNOT p).xt = xor p.xt p.xc ∧ (conjCNOT p).zt = p.zt :=
  ⟨rfl, rfl, rfl, rfl⟩

/-- **The `CNOT`-conjugation phase delta** `β(CNOT,P) = x_c ∧ z_t ∧ (x_t = z_c)`,
matching `cnot_phase` (`phase/clifford.rs:80`). -/
theorem conjCNOT_sign (p : TwoPauli) :
    (conjCNOT p).phase
      = p.phase + (if p.xc && p.zt && (p.xt == p.zc) then 2 else 0) := rfl

/-- `conjCZ` acts on the bits by `Symplectic.czAct`: `z_c ⊕= x_t`, `z_t ⊕= x_c`. -/
theorem conjCZ_bits (p : TwoPauli) :
    (conjCZ p).xc = p.xc ∧ (conjCZ p).zc = xor p.zc p.xt
      ∧ (conjCZ p).xt = p.xt ∧ (conjCZ p).zt = xor p.zt p.xc :=
  ⟨rfl, rfl, rfl, rfl⟩

/-- **The `CZ`-conjugation phase delta** `β(CZ,P) = x_c ∧ x_t ∧ (z_c ≠ z_t)`,
matching `cz_phase` (`phase/clifford.rs:95`). -/
theorem conjCZ_sign (p : TwoPauli) :
    (conjCZ p).phase
      = p.phase + (if p.xc && p.xt && xor p.zc p.zt then 2 else 0) := rfl

/-! ### Fixing the maps to genuine conjugation: the generator tables

Two homomorphisms that agree on a generating set are equal, so these basis
identities (the standard `CNOT`/`CZ` tableau tables) fix `conjCNOT`/`conjCZ` as
`G·P·G†`; the phase deltas above are then forced, not chosen. -/

/-- `CNOT`: `X_c → X_cX_t`, `X_t → X_t`, `Z_c → Z_c`, `Z_t → Z_cZ_t` (all sign-free
on the generators). -/
theorem conjCNOT_Xc :
    conjCNOT ⟨0, true, false, false, false⟩ = ⟨0, true, false, true, false⟩ := by decide
theorem conjCNOT_Xt :
    conjCNOT ⟨0, false, false, true, false⟩ = ⟨0, false, false, true, false⟩ := by decide
theorem conjCNOT_Zc :
    conjCNOT ⟨0, false, true, false, false⟩ = ⟨0, false, true, false, false⟩ := by decide
theorem conjCNOT_Zt :
    conjCNOT ⟨0, false, false, false, true⟩ = ⟨0, false, true, false, true⟩ := by decide

/-- A forced sign: `CNOT · (Y_c ⊗ Y_t) · CNOT = −(X_c Z_c)⊗(X_t Z_t)` bit-pattern
carries the `i²` the delta predicts (`Y_c Y_t → −X_c Y_t Z_c`, i.e. the phase
advances by `2`). -/
theorem conjCNOT_YcYt :
    conjCNOT ⟨0, true, true, true, true⟩ = ⟨2, true, false, false, true⟩ := by decide

/-- `CZ`: `X_c → X_cZ_t`, `X_t → Z_cX_t`, `Z_c → Z_c`, `Z_t → Z_t`. -/
theorem conjCZ_Xc :
    conjCZ ⟨0, true, false, false, false⟩ = ⟨0, true, false, false, true⟩ := by decide
theorem conjCZ_Xt :
    conjCZ ⟨0, false, false, true, false⟩ = ⟨0, false, true, true, false⟩ := by decide
theorem conjCZ_Zc :
    conjCZ ⟨0, false, true, false, false⟩ = ⟨0, false, true, false, false⟩ := by decide
theorem conjCZ_Zt :
    conjCZ ⟨0, false, false, false, true⟩ = ⟨0, false, false, false, true⟩ := by decide

/-- A forced `CZ` sign on `Y_c ⊗ X_t`: `CZ·(Y_c X_t)·CZ = −X_c Y_t` — the delta
fires (`z_c ≠ z_t`), advancing the phase by `2`. -/
theorem conjCZ_YcXt :
    conjCZ ⟨0, true, true, true, false⟩ = ⟨2, true, false, true, true⟩ := by decide

/-! ### The two-qubit conjugations never emit `±i` either

Same reachability invariant as the single-qubit generators: the `CNOT`/`CZ`
conjugation deltas are the real `if … then 2 else 0` (`conjCNOT_sign`/
`conjCZ_sign`), so a real input phase stays real. This is what makes
`ppvm-pauli-sum-2`'s `clifford_sign` drain total on the two-qubit gates too. -/

/-- **`CNOT`-conjugation keeps the phase real** (`±1`, never `±i`). -/
theorem conjCNOT_isRealPhase {p : TwoPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjCNOT p).phase := by
  rw [conjCNOT_sign]; exact hp.add (isRealPhase_ite_two _)

/-- **`CZ`-conjugation keeps the phase real.** -/
theorem conjCZ_isRealPhase {p : TwoPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjCZ p).phase := by
  rw [conjCZ_sign]; exact hp.add (isRealPhase_ite_two _)

/-! ### `CY` as a product of the audited generators

`crates/ppvm-traits-2/src/pauli.rs` emits `cy(c,t)` as `s(t) ; cnot(c,t) ;
s_dag(t)` — and `s_dag(t)` is itself `s(t) ; z(t)` — i.e. the operator identity
`CY = (I⊗S)·CNOT·(I⊗S†)`. Nothing about that composite's *phase* follows from
`conjCNOT_sign` alone, so it is discharged here: we lift the single-qubit `S`,
`S†`, `Z` conjugations to the **target factor** of `𝒫₂` (the control factor is
untouched), form the crate's exact call sequence, and check the resulting map
against the old crate's 16-entry `CY` table
(`crates/ppvm-traits/src/traits/clifford.rs`, replayed in
`crates/ppvm-traits-2/tests/phase1_gate_surface.rs`
`blanket_cy_matches_old_two_qubit_table`).

As with the single-qubit extensions, the crate's `s` gate is the **backward**
`conjSdag`; its `z` gate is self-adjoint. -/

/-- The `S†`-conjugation phase delta `x ∧ ¬z` (`conjSdag_sign`), named so the
homomorphism proofs treat it as one atom (as `cnotDelta` is). -/
def sdagDelta (x z : Bool) : ZMod 4 := if x && !z then 2 else 0

/-- The `S`-conjugation phase delta `x ∧ z` (`conjS_sign`), likewise named. -/
def sDelta (x z : Bool) : ZMod 4 := if x && z then 2 else 0

/-- The crate's `s(t)`: single-qubit `S†`-conjugation on the target factor of
`𝒫₂`, the control factor untouched. -/
def conjSdagT (p : TwoPauli) : TwoPauli :=
  ⟨p.phase + sdagDelta p.xt p.zt, p.xc, p.zc, p.xt, xor p.xt p.zt⟩

/-- Single-qubit `S`-conjugation on the target factor. -/
def conjST (p : TwoPauli) : TwoPauli :=
  ⟨p.phase + sDelta p.xt p.zt, p.xc, p.zc, p.xt, xor p.xt p.zt⟩

/-- The crate's `z(t)`: single-qubit `Z`-conjugation on the target factor (word
fixed, phase advanced by `2·x_t`). -/
def conjZT (p : TwoPauli) : TwoPauli :=
  ⟨p.phase + (if p.xt then 2 else 0), p.xc, p.zc, p.xt, p.zt⟩

-- Same structural route as `conjCNOT`/`conjCZ`: the global phases cancel and the
-- control cocycle is untouched, leaving a 4-boolean identity on the target.

/-- The `S†`-on-target delta is compatible with the target cocycle (4 booleans). -/
private theorem conjSdagT_phase_delta : ∀ c d g h : Bool,
    phaseExp c d g h + sdagDelta (xor c g) (xor d h)
      = sdagDelta c d + sdagDelta g h + phaseExp c (xor c d) g (xor g h) := by decide

/-- The `S`-on-target delta is compatible with the target cocycle. -/
private theorem conjST_phase_delta : ∀ c d g h : Bool,
    phaseExp c d g h + sDelta (xor c g) (xor d h)
      = sDelta c d + sDelta g h + phaseExp c (xor c d) g (xor g h) := by decide

private theorem conjSdagT_mul_raw :
    ∀ p q, conjSdagT (mul p q) = mul (conjSdagT p) (conjSdagT q) := by
  rintro ⟨φ, a, b, c, d⟩ ⟨ψ, e, f, g, h⟩
  ext
  · simp only [conjSdagT, mul]
    linear_combination conjSdagT_phase_delta c d g h
  all_goals (simp only [conjSdagT, mul]; try (revert a b c d e f g h; decide))

private theorem conjST_mul_raw :
    ∀ p q, conjST (mul p q) = mul (conjST p) (conjST q) := by
  rintro ⟨φ, a, b, c, d⟩ ⟨ψ, e, f, g, h⟩
  ext
  · simp only [conjST, mul]
    linear_combination conjST_phase_delta c d g h
  all_goals (simp only [conjST, mul]; try (revert a b c d e f g h; decide))

/-- **`S†`-on-target is a group homomorphism of `𝒫₂`.** -/
def conjSdagTHom : TwoPauli →* TwoPauli where
  toFun := conjSdagT
  map_one' := by decide
  map_mul' := conjSdagT_mul_raw

/-- **`S`-on-target is a group homomorphism of `𝒫₂`.** -/
def conjSTHom : TwoPauli →* TwoPauli where
  toFun := conjST
  map_one' := by decide
  map_mul' := conjST_mul_raw

/-- The crate's `s_dag(t)` (`s(t)` then `z(t)`) collapses to `S`-conjugation on
the target — the two-qubit form of `extSdag_eq_conjS`. -/
theorem conjZT_conjSdagT_eq_conjST : ∀ p, conjZT (conjSdagT p) = conjST p := by decide

/-- **`CY`-conjugation**, defined as the crate's generator product
`s(t) ; cnot(c,t) ; s_dag(t)`. -/
def conjCY (p : TwoPauli) : TwoPauli := conjST (conjCNOT (conjSdagT p))

/-- The *literal* four-primitive call sequence `ppvm-traits-2` emits
(`s(t)`, `cnot`, `s(t)`, `z(t)` — the last two being the inlined `s_dag(t)`)
equals `conjCY`. -/
theorem conjCY_calls : ∀ p, conjZT (conjSdagT (conjCNOT (conjSdagT p))) = conjCY p := by
  decide

/-- **`CY`-conjugation is a group homomorphism of `𝒫₂`** — by construction, the
composite of three audited generator homs. No new sign algebra is introduced. -/
def conjCYHom : TwoPauli →* TwoPauli := conjSTHom.comp (conjCNOTHom.comp conjSdagTHom)

@[simp] theorem conjCYHom_apply (p : TwoPauli) : conjCYHom p = conjCY p := rfl

/-- **The `CY` bit rule falls out of the product**: `z_c ⊕= x_t ⊕ z_t`,
`x_t ⊕= x_c`, `z_t ⊕= x_c`, control X untouched — exactly the old crate's atomic
`fn cy` (`crates/ppvm-traits/src/traits/clifford.rs`) and the `cyAct` of
`lean/PPVM/Pauli/Symplectic.lean`. -/
theorem conjCY_bits : ∀ p : TwoPauli,
    (conjCY p).xc = p.xc ∧ (conjCY p).zc = xor p.zc (xor p.xt p.zt)
      ∧ (conjCY p).xt = xor p.xt p.xc ∧ (conjCY p).zt = xor p.zt p.xc := by decide

/-- **The `CY` sign falls out of the product**: starting from `+1`, the composite
phase is `−1` on exactly two of the sixteen input Paulis, `X⊗X` and `Y⊗Z`. With
`conjCY_bits` this *is* the old crate's 16-entry table (`XX ↦ −YZ`,
`YZ ↦ −XX`, the other fourteen sign-free). -/
theorem conjCY_sign : ∀ xc zc xt zt : Bool,
    (conjCY ⟨0, xc, zc, xt, zt⟩).phase
      = (if (xc && !zc && xt && !zt) || (xc && zc && !xt && zt) then 2 else 0) := by
  decide

/-- `CY` on the four generators: `X_c ↦ X_cY_t`, `Z_c ↦ Z_c`, `X_t ↦ Z_cX_t`,
`Z_t ↦ Z_cZ_t` — all sign-free, the standard tableau table. -/
theorem conjCY_Xc :
    conjCY ⟨0, true, false, false, false⟩ = ⟨0, true, false, true, true⟩ := by decide
theorem conjCY_Zc :
    conjCY ⟨0, false, true, false, false⟩ = ⟨0, false, true, false, false⟩ := by decide
theorem conjCY_Xt :
    conjCY ⟨0, false, false, true, false⟩ = ⟨0, false, true, true, false⟩ := by decide
theorem conjCY_Zt :
    conjCY ⟨0, false, false, false, true⟩ = ⟨0, false, true, false, true⟩ := by decide

/-- The two forced signs: `X_cX_t ↦ −Y_cZ_t` and `Y_cZ_t ↦ −X_cX_t`. -/
theorem conjCY_XcXt :
    conjCY ⟨0, true, false, true, false⟩ = ⟨2, true, true, false, true⟩ := by decide
theorem conjCY_YcZt :
    conjCY ⟨0, true, true, false, true⟩ = ⟨2, true, false, true, false⟩ := by decide

/-- **`CY`-conjugation keeps the phase real** — the composite of three real-delta
maps, so `ppvm-pauli-sum-2`'s `±1` drain stays total on `cy` too. -/
theorem conjCY_isRealPhase {p : TwoPauli} (hp : IsRealPhase p.phase) :
    IsRealPhase (conjCY p).phase :=
  (conjCNOT_isRealPhase (hp.add (isRealPhase_ite_two _))).add (isRealPhase_ite_two _)

end TwoPauli

end PPVM.PauliPhase
