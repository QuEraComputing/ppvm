// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Clifford`] propagation for a Pauli-keyed [`Sum`] (the `PauliSum` alias),
//! draining each conjugation's `±1` sign to the coefficient. Two fast paths:
//!
//! * The **pure-sign** single-qubit gates `X`/`Y`/`Z` leave the Pauli word fixed
//!   (`XPX = (−1)^z P`, `YPY = (−1)^{x⊕z} P`, `ZPZ = (−1)^x P` — proven as the
//!   group conjugation `G·P·G⁻¹` in `lean/PPVM/Pauli/Conjugation.lean`,
//!   `conjX`/`conjY`/`conjZ`), so they take the in-place
//!   [`Sum::flip_sign_by_key`](crate::Sum) path — walk the existing entries and
//!   scale each coefficient by the `±1` its own bits demand, **no** map rebuild,
//!   no key movement, no reallocation (the old crate's in-place `scale`, restored).
//! * The **word-changing** gates `H`/`S`/`CNOT`/`CZ` re-key every term via the
//!   move-based [`Sum::rekey_bijective`](crate::Sum) fast path.
//!
//! # The Clifford-sign subtlety
//!
//! A Clifford conjugates each Pauli to `±` another Pauli. The **bare**
//! [`PauliWord`] Clifford (the blanket over `SymplecticColumns` + `PhaseTrack`)
//! is *bit-only*: its `PhaseTrack` is a no-op, so it computes the resulting
//! Pauli's bits but **drops the `±` sign**. The design's generic
//! `impl<S> Clifford for Sum<S, P> where S::Key: Clifford` would therefore lose
//! every conjugation sign for a `PauliSum` — a silent correctness bug.
//!
//! So this impl does **not** dispatch to the key's bare `Clifford`. For each key
//! `w` it wraps `w` in a [`Phased`]`<PauliWord>` at phase `+1`, applies the gate
//! through the phased word's **audited fused** `Clifford` (which *does* track the
//! `ℤ₄` sign — `ppvm-phased-pauli-word-2`), extracts the resulting phase (a
//! Clifford never emits `i`, only `ℤ₄ ∈ {+1, −1}`), and multiplies the
//! coefficient by that `±1` via [`Coefficient::mul_sign`]. The re-keyed term is
//! keyed on the bare, phase-stripped word. No `ImaginaryUnit` capability is
//! needed — the phase is a pure real sign. That "never emits `i`" fact is
//! machine-checked: every generator's conjugation delta is even in `ℤ₄` (`∈ {0,
//! 2}`), so a phase starting at `+1` stays real, making `clifford_sign`'s `±1`
//! drain total and its `PosI`/`NegI` branch unreachable
//! (`lean/PPVM/Pauli/Conjugation.lean` `conjH_isRealPhase`/`conjS_isRealPhase`/
//! `conjSdag_isRealPhase`/`conjCNOT_isRealPhase`/`conjCZ_isRealPhase` over
//! `IsRealPhase`).
//!
//! A Clifford re-key is a bijection, so colliding re-keyed terms never occur —
//! which is what lets this path take the plain-`insert`
//! [`Sum::rekey_bijective`] fast path rather than [`Sum::apply`]'s batch
//! round-trip. Nothing here runs `reduce` (a `±1` sign cannot zero a
//! coefficient, and a zero term must survive regardless — old has no `reduce`)
//! and nothing here truncates (caller-driven, [`Sum::truncate`]). The bijection is
//! machine-checked at the phase-stripped word level: each generator's `Sp(2n, 2)`
//! bit map is an involution, hence bijective (`lean/PPVM/Pauli/Symplectic.lean`
//! `hAct_involutive`/`sAct_involutive`/`cnotAct_involutive`/`czAct_involutive`,
//! `*_bijective`).
//!
//! Composing that bijectivity with the `±1`-sign reality above gives the
//! semantic guarantee that ties this Clifford path to [`Sum::overlap`]: the
//! Heisenberg re-key **preserves the Hilbert–Schmidt trace pairing**
//! (`overlap(conj_G A, conj_G B) = overlap(A, B)`), because the drained sign
//! squares out (`s_P² = 1`) while the bijection only permutes the summands. This
//! is machine-checked as `clifford_conjugation_preserves_overlap` in
//! `lean/PPVM/Algebra/GradedMap.lean` (over `overlap_eq_fintype_sum`).
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits"
//! (`PauliSum` is deliberately *not* a `BlanketClifford` implementer: "the sum
//! applies the one-row action pointwise and drains each term's phase delta to its
//! coefficient") and §"apply". The conjugation signs are machine-checked in
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_Y`: `HYH = −Y`, `conjSdag_sign`,
//! `conjCNOT_sign`, `conjCZ_sign`); the underlying bit maps are the `Sp(2n, 2)`
//! isometries of `lean/PPVM/Pauli/Symplectic.lean`.

use std::hash::BuildHasher;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage, PauliWord};
use ppvm_phased_pauli_word_2::Phased;
use ppvm_traits_2::{Accumulate, Clifford, Coefficient, PauliBits, Phase, Retain};

use crate::policy::Policy;
use crate::store::{RekeyBijective, SignFlipByKey, StoreAlloc};
use crate::sum::Sum;

/// The `±1` a Clifford conjugation puts on the coefficient, read off the phase
/// the fused phased-word `Clifford` accumulated from a `+1`-phase input.
///
/// A Clifford maps a Pauli to `±` a Pauli, so the phase is always `Pos1` or
/// `Neg1`; an imaginary result would indicate a bug in the fused kernel.
#[inline]
fn clifford_sign(phase: Phase) -> i8 {
    match phase {
        Phase::Pos1 => 1,
        Phase::Neg1 => -1,
        Phase::PosI | Phase::NegI => {
            debug_assert!(false, "Clifford conjugation emitted an imaginary phase");
            1
        }
    }
}

impl<S, P, A, H, C> Sum<S, P>
where
    S: Accumulate<Key = PauliWord<A, H>, Coeff = C>
        + StoreAlloc
        + Retain<PauliWord<A, H>, C>
        + RekeyBijective<PauliWord<A, H>, C>,
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
    C: Coefficient,
    P: Policy<PauliWord<A, H>, C>,
{
    /// Re-key every term by conjugating its Pauli with `gate` (run on a phased
    /// wrapper), draining the resulting `±1` sign to the coefficient.
    ///
    /// Each key is **moved** into the [`Phased`] wrapper (`Phased::new` takes the
    /// word by value), conjugated in place, and moved back out — no key clone.
    /// The whole re-key runs through [`Sum::rekey_bijective`], the move-based
    /// fast path that reuses the support's allocation and skips the batch
    /// round-trip.
    #[inline]
    fn rekey_clifford<G>(&mut self, gate: G)
    where
        G: Fn(&mut Phased<PauliWord<A, H>>),
    {
        self.rekey_bijective(|k: PauliWord<A, H>, c: C| {
            let mut p = Phased::new(k);
            gate(&mut p);
            let (word, phase) = p.into_parts();
            (word, c.mul_sign(clifford_sign(phase)))
        });
    }
}

/// Clifford propagation on a Pauli-keyed `Sum`. Each gate re-keys the whole
/// support pointwise, folding the conjugation sign into the coefficient.
impl<S, P, A, H, C> Clifford for Sum<S, P>
where
    S: Accumulate<Key = PauliWord<A, H>, Coeff = C>
        + StoreAlloc
        + Retain<PauliWord<A, H>, C>
        + RekeyBijective<PauliWord<A, H>, C>
        + SignFlipByKey<PauliWord<A, H>, C>,
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
    C: Coefficient,
    P: Policy<PauliWord<A, H>, C>,
{
    /// `X` conjugation is a **pure sign**: `XPX = (−1)^z P`. The word is fixed, so
    /// this takes the in-place [`Sum::flip_sign_by_key`] fast path — flipping each
    /// term's coefficient iff its `z` bit at `qubit` is set — instead of rebuilding
    /// the map. Sign matches [`PhaseTrack::x_phase`](ppvm_traits_2::PhaseTrack) and
    /// the phased word's fused `Phased::x` (`ppvm-phased-pauli-word-2`).
    #[inline]
    fn x(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && k.z_bit(qubit) {
                -1
            } else {
                1
            }
        });
    }

    /// `Y` conjugation is a **pure sign**: `YPY = (−1)^{x⊕z} P`. Word fixed → the
    /// in-place fast path, flipping iff `x ⊕ z` at `qubit`. Sign matches
    /// [`PhaseTrack::y_phase`](ppvm_traits_2::PhaseTrack) and the phased word's
    /// fused `Phased::y`.
    #[inline]
    fn y(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && (k.x_bit(qubit) ^ k.z_bit(qubit)) {
                -1
            } else {
                1
            }
        });
    }

    /// `Z` conjugation is a **pure sign**: `ZPZ = (−1)^x P`. Word fixed → the
    /// in-place fast path, flipping iff the `x` bit at `qubit` is set. Sign matches
    /// [`PhaseTrack::z_phase`](ppvm_traits_2::PhaseTrack) and the phased word's
    /// fused `Phased::z`.
    #[inline]
    fn z(&mut self, qubit: usize) {
        self.flip_sign_by_key(move |k| {
            if !k.is_lost(qubit) && k.x_bit(qubit) {
                -1
            } else {
                1
            }
        });
    }

    #[inline]
    fn h(&mut self, qubit: usize) {
        self.rekey_clifford(move |p| p.h(qubit));
    }

    #[inline]
    fn s(&mut self, qubit: usize) {
        self.rekey_clifford(move |p| p.s(qubit));
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        self.rekey_clifford(move |p| p.cnot(control, target));
    }

    #[inline]
    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        self.rekey_clifford(move |p| p.cz(qubit0, qubit1));
    }
}
