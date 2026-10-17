// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! X/Z bit operations skip lost qubits and invalidate only the X/Z hash.
//! Lossy words store no phase, so phase updates do nothing.

use ppvm_pauli_word_2::PauliStorage;
use ppvm_traits_2::{PhaseTrack, SymplecticColumns};

use crate::LossyPauliWord;

impl<A: PauliStorage, H> SymplecticColumns for LossyPauliWord<A, H> {
    fn n_qubits(&self) -> usize {
        self.nqubits
    }

    fn swap_xz(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.is_lost(q) {
            return;
        }
        let (x, z) = (self.xbits[q], self.zbits[q]);
        self.xbits.set(q, z);
        self.zbits.set(q, x);
        self.invalidate_xz_hash();
    }

    fn xor_z_from_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.is_lost(q) {
            return;
        }
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.invalidate_xz_hash();
    }

    fn xor_x_col(&mut self, control: usize, target: usize) {
        debug_assert!(
            control < self.nqubits && target < self.nqubits,
            "qubit out of bounds"
        );
        if self.is_lost(control) || self.is_lost(target) {
            return;
        }
        let target_x = self.xbits[target] ^ self.xbits[control];
        self.xbits.set(target, target_x);
        self.invalidate_xz_hash();
    }

    fn xor_z_col(&mut self, target: usize, control: usize) {
        debug_assert!(
            control < self.nqubits && target < self.nqubits,
            "qubit out of bounds"
        );
        if self.is_lost(control) || self.is_lost(target) {
            return;
        }
        let control_z = self.zbits[control] ^ self.zbits[target];
        self.zbits.set(control, control_z);
        self.invalidate_xz_hash();
    }

    fn cz_bits(&mut self, a: usize, b: usize) {
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        if self.is_lost(a) || self.is_lost(b) {
            return;
        }
        let z_a = self.zbits[a] ^ self.xbits[b];
        let z_b = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, z_a);
        self.zbits.set(b, z_b);
        self.invalidate_xz_hash();
    }
}

/// Phase updates are no-ops because the word stores only X, Z, and loss bits.
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
