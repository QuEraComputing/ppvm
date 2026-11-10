// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford conjugation of packed Pauli words, discarding all phase signs.
//! `PauliWord` implements `ApplyGateRule`, so `Clifford` comes from the blanket impl
//! and each gate's bits from `ppvm_pauli_kernels::heisenberg`. The hash is
//! refreshed once per gate, and only when a bit changes.
//! Symplectic bit operations remain available; `PhaseTrack` is a no-op.
//! Phase signs must be handled separately by phased words or the sum engine.

use std::hash::BuildHasher;

use ppvm_traits_2::{ApplyGateRule, PhaseTrack, SymplecticColumns};

use crate::data::{PauliWord, xor_bits2};
use crate::{HashFinalize, PauliStorage};

/// Whether a gate that changes no bit skips the hash refresh. Measured on arm64:
/// for one storage word the hash is cheaper than the unpredictable branch, for
/// larger storage skipping it is faster.
const fn skips_unchanged<A>() -> bool {
    std::mem::size_of::<A>() > 8
}

impl<A, H> ApplyGateRule for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    type B = bool;

    /// The sign flip is dropped: an unsigned word has no sign. Only changed bits
    /// are written, and the hash is refreshed once if any changed.
    #[inline]
    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let (x, z) = (self.xbits[q], self.zbits[q]);
        let ([x2, z2], _) = rule(x, z);
        let (tx, tz) = (x ^ x2, z ^ z2);
        if skips_unchanged::<A>() && !(tx || tz) {
            return;
        }
        xor_bits2(&mut self.xbits, q, tx, q, false);
        xor_bits2(&mut self.zbits, q, tz, q, false);
        self.refresh_hash();
    }

    #[inline]
    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        debug_assert_ne!(a, b, "two-qubit gates require distinct qubits");
        debug_assert!(a < self.nqubits && b < self.nqubits, "qubit out of bounds");
        let old = [self.xbits[a], self.zbits[a], self.xbits[b], self.zbits[b]];
        let (new, _) = rule(old);
        let [txa, tza, txb, tzb] = [0, 1, 2, 3].map(|i| old[i] ^ new[i]);
        if skips_unchanged::<A>() && !(txa || tza || txb || tzb) {
            return;
        }
        xor_bits2(&mut self.xbits, a, txa, b, txb);
        xor_bits2(&mut self.zbits, a, tza, b, tzb);
        self.refresh_hash();
    }

    /// X, Y and Z only flip signs, so they leave an unsigned word unchanged.
    #[inline]
    fn apply_sign1(&mut self, q: usize, _rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
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

    /// After any gate sequence the cached hash equals that of a freshly built word.
    fn hash_stays_structural<A: PauliStorage>(n: usize)
    where
        PauliWord<A>: Clifford,
    {
        let mut next = rng(0x9E37_79B9_7F4A_7C15);
        for _ in 0..64 {
            let start: String = (0..n)
                .map(|_| ['I', 'X', 'Y', 'Z'][(next() % 4) as usize])
                .collect();
            let mut w = PauliWord::<A>::from(start.as_str());
            for _ in 0..64 {
                let a = (next() % n as u64) as usize;
                let b = (a + 1 + (next() % (n as u64 - 1)) as usize) % n;
                random_gate(&mut w, next(), a, b);
                let fresh = PauliWord::<A>::from(w.to_string().as_str());
                assert_eq!(w, fresh);
                assert_eq!(w.key_hash(), fresh.key_hash());
            }
        }
    }

    /// One-word storage always writes; larger storage skips gates that change no
    /// bit. Both must agree on every circuit.
    #[test]
    fn storage_sizes_agree() {
        let mut next = rng(0xD1B5_4A32_D192_ED03);
        for _ in 0..64 {
            let start: String = (0..64)
                .map(|_| ['I', 'X', 'Y', 'Z'][(next() % 4) as usize])
                .collect();
            let mut small = PauliWord::<u64>::from(start.as_str());
            let mut large = PauliWord::<[u64; 4]>::from(start.as_str());
            for _ in 0..256 {
                let (g, a) = (next(), (next() % 64) as usize);
                let b = (a + 1 + (next() % 63) as usize) % 64;
                random_gate(&mut small, g, a, b);
                random_gate(&mut large, g, a, b);
            }
            assert_eq!(small.to_string(), large.to_string());
        }
    }

    #[test]
    fn gates_keep_the_hash_structural() {
        // One storage word always refreshes; larger storage skips unchanged gates.
        hash_stays_structural::<u64>(5);
        hash_stays_structural::<[u64; 4]>(200);
    }
}
