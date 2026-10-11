// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford conjugation of packed Pauli words, discarding all phase signs.
//! Fused gates refresh the eager hash once instead of hashing each primitive step.
//! Symplectic bit operations remain available; `PhaseTrack` is a no-op.
//! Phase signs must be handled separately by phased words or the sum engine.

use std::hash::BuildHasher;

use ppvm_traits_2::{Clifford, PhaseTrack, SymplecticColumns};

use crate::data::PauliWord;
use crate::{HashFinalize, PauliStorage};

impl<A, H> Clifford for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
    }

    #[inline]
    fn y(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
    }

    #[inline]
    fn z(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
    }

    #[inline]
    fn h(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let x = self.xbits[q];
        self.xbits.set(q, self.zbits[q]);
        self.zbits.set(q, x);
        self.refresh_hash();
    }

    #[inline]
    fn s(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.refresh_hash();
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        debug_assert_ne!(control, target, "CNOT requires distinct qubits");
        debug_assert!(
            control < self.nqubits && target < self.nqubits,
            "qubit out of bounds"
        );
        let xt = self.xbits[target] ^ self.xbits[control];
        let zc = self.zbits[control] ^ self.zbits[target];
        self.xbits.set(target, xt);
        self.zbits.set(control, zc);
        self.refresh_hash();
    }

    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        debug_assert_ne!(a, b, "CZ requires distinct qubits");
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        let za = self.zbits[a] ^ self.xbits[b];
        let zb = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, za);
        self.zbits.set(b, zb);
        self.refresh_hash();
    }
    #[inline]
    fn s_dag(&mut self, q: usize) {
        self.s(q);
    }

    #[inline]
    fn sqrt_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let x = self.xbits[q] ^ self.zbits[q];
        self.xbits.set(q, x);
        self.refresh_hash();
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
        debug_assert_ne!(control, target, "CY requires distinct qubits");
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
        self.refresh_hash();
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
        self.refresh_hash();
    }

    /// `S` on `q`: `z_q ⊕= x_q` (maps `X → Y`).
    #[inline]
    fn xor_z_from_x(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let z = self.zbits[q] ^ self.xbits[q];
        self.zbits.set(q, z);
        self.refresh_hash();
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
        self.refresh_hash();
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
        self.refresh_hash();
    }

    /// `CZ` on `(a, b)`: `z_a ⊕= x_b` and `z_b ⊕= x_a`.
    #[inline]
    fn cz_bits(&mut self, a: usize, b: usize) {
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        let za = self.zbits[a] ^ self.xbits[b];
        let zb = self.zbits[b] ^ self.xbits[a];
        self.zbits.set(a, za);
        self.zbits.set(b, zb);
        self.refresh_hash();
    }
}

/// Bare words store no phase, so phase updates deliberately do nothing.
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
    use ppvm_traits_2::{Clifford, Indexable};

    fn assert_rejected_unchanged(input: &str, operation: impl FnOnce(&mut PauliWord)) {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let mut word = PauliWord::from(input);
        let original = word;
        let hash = word.key_hash();
        assert!(catch_unwind(AssertUnwindSafe(|| operation(&mut word))).is_err());
        assert_eq!(word, original);
        assert_eq!(word.key_hash(), hash);
    }

    #[test]
    fn single_qubit_operations_enforce_logical_bounds() {
        let operations: [fn(&mut PauliWord, usize); 12] = [
            PauliWord::x,
            PauliWord::y,
            PauliWord::z,
            PauliWord::h,
            PauliWord::s,
            PauliWord::s_dag,
            PauliWord::sqrt_x,
            PauliWord::sqrt_x_dag,
            PauliWord::sqrt_y,
            PauliWord::sqrt_y_dag,
            PauliWord::swap_xz,
            PauliWord::xor_z_from_x,
        ];
        for operation in operations {
            for invalid in [1, 7, 63] {
                assert_rejected_unchanged("X", |word| operation(word, invalid));
            }
            assert_rejected_unchanged("", |word| operation(word, 0));
        }
    }

    #[test]
    fn two_qubit_operations_enforce_logical_bounds() {
        let operations: [fn(&mut PauliWord, usize, usize); 6] = [
            PauliWord::cnot,
            PauliWord::cy,
            PauliWord::cz,
            PauliWord::xor_x_col,
            PauliWord::xor_z_col,
            PauliWord::cz_bits,
        ];
        for operation in operations {
            for input in ["X", "Y", "Z"] {
                for invalid in [1, 7, 63] {
                    assert_rejected_unchanged(input, |word| operation(word, 0, invalid));
                    assert_rejected_unchanged(input, |word| operation(word, invalid, 0));
                }
            }
        }
    }

    #[test]
    fn two_qubit_gates_reject_identical_indices_before_mutation() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let gates: [fn(&mut PauliWord, usize, usize); 3] = [
            <PauliWord as Clifford>::cnot,
            <PauliWord as Clifford>::cy,
            <PauliWord as Clifford>::cz,
        ];
        for gate in gates {
            for input in ["I", "X", "Y", "Z"] {
                let mut word = PauliWord::from(input);
                let original = word;
                let hash = word.key_hash();
                assert!(catch_unwind(AssertUnwindSafe(|| gate(&mut word, 0, 0))).is_err());
                assert_eq!(word, original);
                assert_eq!(word.key_hash(), hash);
            }
        }
    }

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
