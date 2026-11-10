// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! X/Z bit operations skip lost qubits and invalidate only the X/Z hash.
//! Lossy words store no phase, so phase updates do nothing.
//!
//! `LossyPauliWord` implements `ApplyGateRule`, so it gets `Clifford` from the
//! blanket impl: a gate touching a lost site does nothing; otherwise the bits come
//! from `ppvm_pauli_kernels::heisenberg` and the sign flip is dropped.

use ppvm_pauli_word_2::{PauliStorage, xor_bits2};
use ppvm_traits_2::{ApplyGateRule, PhaseTrack, SymplecticColumns};

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

impl<A: PauliStorage, H> ApplyGateRule for LossyPauliWord<A, H> {
    type B = bool;

    #[inline]
    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        if self.is_lost(q) {
            return;
        }
        let (x, z) = (self.xbits[q], self.zbits[q]);
        let ([x2, z2], _) = rule(x, z);
        let (tx, tz) = (x ^ x2, z ^ z2);
        xor_bits2(&mut self.xbits, q, tx, q, false);
        xor_bits2(&mut self.zbits, q, tz, q, false);
        // Branch-free: the lazy hash stays valid only if nothing changed.
        *self.xz_hash_valid.get_mut() &= !(tx || tz);
    }

    #[inline]
    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        debug_assert_ne!(a, b, "two-qubit gates require distinct qubits");
        if self.is_lost(a) || self.is_lost(b) {
            return;
        }
        let old = [self.xbits[a], self.zbits[a], self.xbits[b], self.zbits[b]];
        let (new, _) = rule(old);
        let [txa, tza, txb, tzb] = [0, 1, 2, 3].map(|i| old[i] ^ new[i]);
        xor_bits2(&mut self.xbits, a, txa, b, txb);
        xor_bits2(&mut self.zbits, a, tza, b, tzb);
        *self.xz_hash_valid.get_mut() &= !(txa || tza || txb || tzb);
    }

    /// X, Y and Z only flip signs, which a lossy word does not store.
    #[inline]
    fn apply_sign1(&mut self, q: usize, _rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.check_qubit(q);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ppvm_pauli_word_2::PauliWord;
    use ppvm_traits_2::{Clifford, Indexable};

    fn random_gate<C: Clifford>(c: &mut C, g: u64, a: usize, b: usize) {
        match g % 13 {
            0 => c.x(a),
            1 => c.y(a),
            2 => c.z(a),
            3 => c.h(a),
            4 => c.s(a),
            5 => c.s_dag(a),
            6 => c.sqrt_x(a),
            7 => c.sqrt_x_dag(a),
            8 => c.sqrt_y(a),
            9 => c.sqrt_y_dag(a),
            10 => c.cnot(a, b),
            11 => c.cz(a, b),
            _ => c.cy(a, b),
        }
    }

    /// Xorshift64, so the tests need no dependency.
    fn rng(mut seed: u64) -> impl FnMut() -> u64 {
        move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        }
    }

    fn random_paulis(next: &mut impl FnMut() -> u64, n: usize) -> String {
        (0..n)
            .map(|_| ['I', 'X', 'Y', 'Z'][(next() % 4) as usize])
            .collect()
    }

    /// Without loss, a lossy word conjugates exactly like a plain `PauliWord`.
    #[test]
    fn loss_free_gates_match_pauli_word() {
        let mut next = rng(0xD1B5_4A32_D192_ED03);
        for _ in 0..64 {
            let start = random_paulis(&mut next, 64);
            let mut lossy = LossyPauliWord::<u64>::from(start.as_str());
            let mut plain = PauliWord::<u64>::from(start.as_str());
            for _ in 0..256 {
                let (g, a) = (next(), (next() % 64) as usize);
                let b = (a + 1 + (next() % 63) as usize) % 64;
                random_gate(&mut lossy, g, a, b);
                random_gate(&mut plain, g, a, b);
            }
            assert_eq!(lossy.to_string(), plain.to_string());
        }
    }

    /// Every gate on a lost site does nothing; gates elsewhere still act.
    #[test]
    fn gates_touching_lost_sites_do_nothing() {
        for g in 0..13 {
            for input in ["XL", "LY", "ZL"] {
                let mut word = LossyPauliWord::<u8>::from(input);
                let lost = if word.is_lost(0) { 0 } else { 1 };
                random_gate(&mut word, g, lost, 1 - lost);
                assert_eq!(word, LossyPauliWord::from(input), "gate {g} on {input}");
            }
        }
        // Gates away from the lost site act as on a plain word.
        let (mut word, mut plain) = (
            LossyPauliWord::<u8>::from("XLZ"),
            PauliWord::<u8>::from("XIZ"),
        );
        word.h(0);
        word.cnot(2, 0);
        plain.h(0);
        plain.cnot(2, 0);
        let expected: String = plain
            .to_string()
            .chars()
            .enumerate()
            .map(|(i, c)| if i == 1 { 'L' } else { c })
            .collect();
        assert_eq!(word.to_string(), expected);
    }

    /// After any circuit, with some sites lost, the cached hash matches a fresh word.
    #[test]
    fn gates_keep_the_hash_structural() {
        let mut next = rng(0x9E37_79B9_7F4A_7C15);
        for _ in 0..64 {
            let mut start: Vec<char> = random_paulis(&mut next, 40).chars().collect();
            for _ in 0..4 {
                start[(next() % 40) as usize] = 'L';
            }
            let start: String = start.into_iter().collect();
            let mut word = LossyPauliWord::<u64>::from(start.as_str());
            for _ in 0..64 {
                let (g, a) = (next(), (next() % 40) as usize);
                let b = (a + 1 + (next() % 39) as usize) % 40;
                random_gate(&mut word, g, a, b);
                let fresh = LossyPauliWord::<u64>::from(word.to_string().as_str());
                assert_eq!(word.key_hash(), fresh.key_hash());
            }
        }
    }

    #[test]
    fn single_qubit_operations_transform_present_sites() {
        let mut word = LossyPauliWord::<u8>::from("XYZL");

        word.swap_xz(0);
        word.swap_xz(2);
        assert_eq!(word.to_string(), "ZYXL");

        word.xor_z_from_x(1);
        word.xor_z_from_x(2);
        assert_eq!(word.to_string(), "ZXYL");
    }

    #[test]
    fn two_qubit_operations_transform_present_sites() {
        let mut word = LossyPauliWord::<u8>::from("XZ");

        word.xor_x_col(0, 1);
        word.xor_z_col(1, 0);
        assert_eq!(word.to_string(), "YY");

        word.cz_bits(0, 1);
        assert_eq!(word.to_string(), "XX");
    }

    #[test]
    fn operations_touching_lost_sites_leave_the_word_unchanged() {
        for input in ["XL", "LY", "LL"] {
            let mut word = LossyPauliWord::<u8>::from(input);
            let lost = if word.is_lost(0) { 0 } else { 1 };

            word.swap_xz(lost);
            word.xor_z_from_x(lost);
            word.xor_x_col(0, 1);
            word.xor_z_col(1, 0);
            word.cz_bits(0, 1);

            assert_eq!(word, LossyPauliWord::from(input));
        }
    }
}
