// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford conjugation for [`PauliWord`]. Public gates use fused packed-bit
//! kernels so each operation refreshes the eager structural digest exactly once.
//! The role-independent [`SymplecticColumns`] and phase-discarding [`PhaseTrack`]
//! primitives remain available for generic algebra, but this type intentionally
//! does not opt into the blanket Clifford implementation: decomposing `CNOT`,
//! `√X`, and `CY` into independently eager primitives would recompute the same
//! width-sensitive digest multiple times.
//!
//! The single-qubit bit rules are ported from that old kernel.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits:
//! symplectic structure and phase". Lean spec: each generator is an `Sp(2n, 2)`
//! isometry on the bit planes (`lean/PPVM/Pauli/Symplectic.lean`
//! `hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/`czAct_isometry`); the
//! conjugation signs a *phased* word would pick up are in
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_Y`: `HYH = −Y`, …), which a bare
//! word deliberately discards.
//!
//! # The phaseless word's no-op `PhaseTrack`
//!
//! The bare word's packed layout has no phase field, so `PhaseTrack` here is a
//! total no-op. The fused public kernels apply the same pure `Sp(2n, 2)` maps and
//! likewise drop every Clifford sign; recovering it is the job of the phased
//! wrapper / the sum's phase-draining path (Phase 3).

use std::hash::BuildHasher;

use ppvm_traits_2::{Clifford, CliffordExtensions, PhaseTrack, SymplecticColumns};

use crate::data::PauliWord;
use crate::storage::{HashFinalize, PauliStorage};

impl<A, H> Clifford for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn x(&mut self, _q: usize) {}

    #[inline]
    fn y(&mut self, _q: usize) {}

    #[inline]
    fn z(&mut self, _q: usize) {}

    #[inline]
    fn h(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let x = self.xbits[q];
        self.xbits.set(q, self.zbits[q]);
        self.zbits.set(q, x);
        self.invalidate_hash();
    }

    #[inline]
    fn s(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.invalidate_hash();
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        debug_assert!(
            control < self.nqubits && target < self.nqubits,
            "qubit out of bounds"
        );
        let xt = self.xbits[target] ^ self.xbits[control];
        let zc = self.zbits[control] ^ self.zbits[target];
        self.xbits.set(target, xt);
        self.zbits.set(control, zc);
        self.invalidate_hash();
    }

    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        let za = self.zbits[a] ^ self.xbits[b];
        let zb = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, za);
        self.zbits.set(b, zb);
        self.invalidate_hash();
    }
}

impl<A, H> CliffordExtensions for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn s_dag(&mut self, q: usize) {
        self.s(q);
    }

    #[inline]
    fn sqrt_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let x = self.xbits[q] ^ self.zbits[q];
        self.xbits.set(q, x);
        self.invalidate_hash();
    }

    #[inline]
    fn sqrt_x_dag(&mut self, q: usize) {
        self.sqrt_x(q);
    }

    #[inline]
    fn sqrt_y(&mut self, q: usize) {
        self.h(q);
    }

    #[inline]
    fn sqrt_y_dag(&mut self, q: usize) {
        self.h(q);
    }

    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        debug_assert!(
            control < self.nqubits && target < self.nqubits,
            "qubit out of bounds"
        );
        let xc = self.xbits[control];
        let zc = self.zbits[control];
        let xt = self.xbits[target];
        let zt = self.zbits[target];
        self.zbits.set(control, zc ^ xt ^ zt);
        self.xbits.set(target, xt ^ xc);
        self.zbits.set(target, zt ^ xc);
        self.invalidate_hash();
    }
}

impl<A, H> SymplecticColumns for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn n_qubits(&self) -> usize {
        self.nqubits
    }

    /// `H` on `q`: swap the X and Z bits. Touches only in-range slot `q`, so the
    /// canonical-unused-bits invariant is preserved.
    #[inline]
    fn swap_xz(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let xb = self.xbits[q];
        let zb = self.zbits[q];
        self.xbits.set(q, zb);
        self.zbits.set(q, xb);
        self.invalidate_hash();
    }

    /// `S` on `q`: `z_q ⊕= x_q` (maps `X → Y`).
    #[inline]
    fn xor_z_from_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.invalidate_hash();
    }

    /// `CNOT` bit rule, part one: `x_tgt ⊕= x_ctrl`.
    #[inline]
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize) {
        debug_assert!(
            ctrl < self.nqubits && tgt < self.nqubits,
            "qubit out of bounds"
        );
        let x = self.xbits[tgt] ^ self.xbits[ctrl];
        self.xbits.set(tgt, x);
        self.invalidate_hash();
    }

    /// `CNOT` bit rule, part two: `z_ctrl ⊕= z_tgt`.
    #[inline]
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize) {
        debug_assert!(
            ctrl < self.nqubits && tgt < self.nqubits,
            "qubit out of bounds"
        );
        let z = self.zbits[ctrl] ^ self.zbits[tgt];
        self.zbits.set(ctrl, z);
        self.invalidate_hash();
    }

    /// `CZ` on `(a, b)`: `z_a ⊕= x_b` and `z_b ⊕= x_a`.
    #[inline]
    fn cz_bits(&mut self, a: usize, b: usize) {
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        let za = self.zbits[a] ^ self.xbits[b];
        let zb = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, za);
        self.zbits.set(b, zb);
        self.invalidate_hash();
    }
}

/// A bare `PauliWord` stores no phase, so every phase delta is a no-op. See the
/// module-level friction note.
impl<A: PauliStorage, H> PhaseTrack for PauliWord<A, H> {
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
    use ppvm_traits_2::Clifford;

    fn conj(input: &str, gate: impl Fn(&mut PauliWord)) -> String {
        let mut w: PauliWord = input.into();
        gate(&mut w);
        w.to_string()
    }

    #[test]
    fn hadamard_bit_map() {
        // Bit-only: X↔Z, Y→Y (the −Y sign a phased word carries is dropped).
        for (input, target) in [("I", "I"), ("X", "Z"), ("Y", "Y"), ("Z", "X")] {
            assert_eq!(conj(input, |w| w.h(0)), target, "H {input}");
        }
    }

    #[test]
    fn phase_gate_bit_map() {
        for (input, target) in [("I", "I"), ("X", "Y"), ("Y", "X"), ("Z", "Z")] {
            assert_eq!(conj(input, |w| w.s(0)), target, "S {input}");
        }
    }

    #[test]
    fn cnot_bit_map() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "IX"),
            ("IZ", "ZZ"),
            ("IY", "ZY"),
            ("XI", "XX"),
            ("XX", "XI"),
            ("XY", "YZ"),
            ("XZ", "YY"),
            ("ZI", "ZI"),
            ("ZX", "ZX"),
            ("ZY", "IY"),
            ("ZZ", "IZ"),
            ("YI", "YX"),
            ("YX", "YI"),
            ("YY", "XZ"),
            ("YZ", "XY"),
        ] {
            assert_eq!(conj(input, |w| w.cnot(0, 1)), target, "CNOT {input}");
        }
    }

    #[test]
    fn cz_bit_map() {
        for (input, target) in [
            ("II", "II"),
            ("IX", "ZX"),
            ("IY", "ZY"),
            ("IZ", "IZ"),
            ("XI", "XZ"),
            ("XX", "YY"),
            ("XY", "YX"),
            ("XZ", "XI"),
            ("ZI", "ZI"),
            ("ZX", "IX"),
            ("ZY", "IY"),
            ("ZZ", "ZZ"),
            ("YI", "YZ"),
            ("YX", "XY"),
            ("YY", "XX"),
            ("YZ", "YI"),
        ] {
            assert_eq!(conj(input, |w| w.cz(0, 1)), target, "CZ {input}");
        }
    }

    #[test]
    fn pauli_gates_are_bit_noops() {
        // X/Y/Z conjugation is pure sign; on a phaseless word it changes nothing.
        for input in ["I", "X", "Y", "Z"] {
            assert_eq!(conj(input, |w| w.x(0)), input);
            assert_eq!(conj(input, |w| w.y(0)), input);
            assert_eq!(conj(input, |w| w.z(0)), input);
        }
    }
}
