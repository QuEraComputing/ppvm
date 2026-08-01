// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford conjugation for [`LossyPauliWord`], split into the role-independent
//! symplectic bit algebra ([`SymplecticColumns`]) and the phase extension
//! ([`PhaseTrack`]). Exactly like the bare `PauliWord`, a `LossyPauliWord`
//! carries **no phase**, so its [`PhaseTrack`] is the no-op that drops every
//! phase delta; the blanket `impl<T: SymplecticColumns + PhaseTrack> Clifford for
//! T` in `ppvm-traits-2` therefore realizes the pure `Sp(2n, 2)` bit map.
//!
//! # Loss guard (ported from `ppvm-pauli-word/src/loss/clifford.rs`)
//!
//! The one behavioral difference from the bare word: **a gate touching a lost
//! qubit is a no-op**, because a lost qubit no longer participates. In the old
//! crate this guard lived in the gate-level blanket `Clifford` impl (`if
//! get_lbit(q) { return }`); here the blanket `Clifford` is algebra-agnostic and
//! calls the column primitives directly, so the guard moves *into* each
//! [`SymplecticColumns`] primitive. `H`/`S` swap or XOR a lost site's already-zero
//! X/Z bits (harmless), but `CNOT`/`CZ` would otherwise write a present control's
//! bit onto a lost target — corrupting the canonical loss invariant `lost ⇒ X/Z
//! identity`. Guarding every primitive on `is_lost` reproduces the old
//! whole-gate skip and keeps the invariant.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits:
//! symplectic structure and phase"; `word-data-structures.md` §"Loss-specific
//! behavior" (loss is preserved/skipped by propagation) and §"Canonical loss
//! invariant". Lean spec: each generator is an `Sp(2n, 2)` isometry on the bit
//! planes (`lean/PPVM/Pauli/Symplectic.lean`
//! `hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/`czAct_isometry`); the
//! conjugation signs a *phased* word would pick up
//! (`lean/PPVM/Pauli/Conjugation.lean`) are deliberately discarded here. The
//! lost-qubit no-op is a hardware-loss modeling choice, but the algebraic
//! invariant it upholds is machine-checked: the loss-guarded generators
//! (`…ActL` in `lean/PPVM/Pauli/Symplectic.lean`) preserve the canonical loss
//! invariant `lost ⇒ X/Z identity` (`hActL_preserves_loss`/…/
//! `czActL_preserves_loss`, critical `CNOT` case
//! `cnotActL_lost_target_stays_identity`) and coincide with the `Sp(2n, 2)`
//! isometry on the present-qubit sub-block (`hActL_present_isometry`/…/
//! `czActL_present_isometry`). The `CNOT` decomposition into two independently
//! guarded columns is checked at primitive granularity: each column preserves
//! the invariant alone (`xorXColL_preserves_loss`/`xorZColL_preserves_loss`) and
//! their composition equals the whole-gate skip
//! (`xorZColL_xorXColL_eq_cnotActL`), so this "reproduces the old whole-gate
//! skip" is a theorem, not just prose. `CY` — which the blanket
//! `CliffordExtensions` decomposes into `s(t); cnot(c,t); s_dag(t)`, three
//! primitives whose guards differ, so a lost control with a present target must
//! *cancel* rather than skip — gets the same treatment:
//! `sActL_cnotActL_sActL_eq_cyActL`, plus `cyActL_preserves_loss` and
//! `cyActL_present_isometry`.

use ppvm_pauli_word_2::PauliStorage;
use ppvm_traits_2::{BlanketClifford, PhaseTrack, SymplecticColumns};

use crate::data::LossyPauliWord;

/// Opt into the single audited blanket `Clifford` (`ppvm-traits-2`): the
/// loss-guarded [`SymplecticColumns`] bit map composed with the phase-discarding
/// [`PhaseTrack`] below.
impl<A: PauliStorage, H> BlanketClifford for LossyPauliWord<A, H> {}

impl<A: PauliStorage, H> SymplecticColumns for LossyPauliWord<A, H> {
    #[inline]
    fn n_qubits(&self) -> usize {
        self.nqubits
    }

    /// `H` on `q`: swap the X and Z bits (no-op on a lost qubit).
    #[inline]
    fn swap_xz(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.is_lost(q) {
            return;
        }
        let xb = self.xbits[q];
        let zb = self.zbits[q];
        self.xbits.set(q, zb);
        self.zbits.set(q, xb);
        self.invalidate_xz();
    }

    /// `S` on `q`: `z_q ⊕= x_q` (no-op on a lost qubit).
    #[inline]
    fn xor_z_from_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.is_lost(q) {
            return;
        }
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.invalidate_xz();
    }

    /// `CNOT` bit rule, part one: `x_tgt ⊕= x_ctrl`. No-op if either qubit is
    /// lost, so the whole `CNOT` skips a lost pair (preserving `lost ⇒ X/Z
    /// identity`).
    #[inline]
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        debug_assert!(
            ctrl < self.nqubits && tgt < self.nqubits,
            "qubit out of bounds"
        );
        if self.is_lost(ctrl) || self.is_lost(tgt) {
            return;
        }
        let x = self.xbits[tgt] ^ self.xbits[ctrl];
        self.xbits.set(tgt, x);
        self.invalidate_xz();
    }

    /// `CNOT` bit rule, part two: `z_ctrl ⊕= z_tgt` (same loss guard as part one).
    #[inline]
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        debug_assert!(
            ctrl < self.nqubits && tgt < self.nqubits,
            "qubit out of bounds"
        );
        if self.is_lost(ctrl) || self.is_lost(tgt) {
            return;
        }
        let z = self.zbits[ctrl] ^ self.zbits[tgt];
        self.zbits.set(ctrl, z);
        self.invalidate_xz();
    }

    /// `CZ` on `(a, b)`: `z_a ⊕= x_b` and `z_b ⊕= x_a` (no-op if either is lost).
    #[inline]
    fn cz_bits(&mut self, a: usize, b: usize) {
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        if self.is_lost(a) || self.is_lost(b) {
            return;
        }
        let za = self.zbits[a] ^ self.xbits[b];
        let zb = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, za);
        self.zbits.set(b, zb);
        self.invalidate_xz();
    }
}

/// A `LossyPauliWord` stores no phase, so every phase delta is a no-op — the
/// same phase-discarding `PhaseTrack` as the bare `PauliWord`.
impl<A: PauliStorage, H> PhaseTrack for LossyPauliWord<A, H> {
    #[inline]
    fn flip_phase_where_xz(&mut self, _q: usize) {}
    #[inline]
    fn s_phase(&mut self, _q: usize) {}
    #[inline]
    fn cnot_phase(&mut self, _ctrl: usize, _tgt: usize) {}
    #[inline]
    fn cz_phase(&mut self, _a: usize, _b: usize) {}
    #[inline]
    fn x_phase(&mut self, _q: usize) {}
    #[inline]
    fn y_phase(&mut self, _q: usize) {}
    #[inline]
    fn z_phase(&mut self, _q: usize) {}
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::{Clifford, PauliBits};

    fn conj(input: &str, gate: impl Fn(&mut LossyPauliWord)) -> String {
        let mut w: LossyPauliWord = input.into();
        gate(&mut w);
        w.to_string()
    }

    #[test]
    fn hadamard_bit_map_and_loss_noop() {
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X"), ("L", "L")] {
            assert_eq!(conj(input, |w| w.h(0)), target, "H {input}");
        }
    }

    #[test]
    fn phase_gate_bit_map_and_loss_noop() {
        for (input, target) in [("I", "I"), ("X", "Y"), ("Y", "X"), ("Z", "Z"), ("L", "L")] {
            assert_eq!(conj(input, |w| w.s(0)), target, "S {input}");
        }
    }

    #[test]
    fn cnot_bit_map_and_loss_noop() {
        for (input, target) in [
            ("II", "II"),
            ("IZ", "ZZ"),
            ("XI", "XX"),
            ("YZ", "XY"),
            // Any lost qubit skips the whole gate.
            ("IL", "IL"),
            ("XL", "XL"),
            ("LI", "LI"),
            ("LX", "LX"),
            ("LL", "LL"),
        ] {
            assert_eq!(conj(input, |w| w.cnot(0, 1)), target, "CNOT {input}");
        }
    }

    #[test]
    fn cz_bit_map_and_loss_noop() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "ZX"),
            ("XX", "YY"),
            ("YZ", "YI"),
            ("XL", "XL"),
            ("LY", "LY"),
            ("LL", "LL"),
        ] {
            assert_eq!(conj(input, |w| w.cz(0, 1)), target, "CZ {input}");
        }
    }

    #[test]
    fn pauli_gates_are_bit_noops() {
        for input in ["I", "X", "Y", "Z", "L"] {
            assert_eq!(conj(input, |w| w.x(0)), input);
            assert_eq!(conj(input, |w| w.y(0)), input);
            assert_eq!(conj(input, |w| w.z(0)), input);
        }
    }

    #[test]
    fn cnot_present_control_lost_target_preserves_invariant() {
        // The critical case: a present X control and a lost target must NOT write
        // the target's X bit (which would break `lost ⇒ X/Z identity`).
        let mut w: LossyPauliWord = "XL".into();
        w.cnot(0, 1);
        assert_eq!(w.to_string(), "XL");
        assert!(
            !w.x_bit(1) && !w.z_bit(1),
            "lost target stayed X/Z identity"
        );
    }
}
