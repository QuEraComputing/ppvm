// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_phased_pauli_word_2::Phased;
use ppvm_traits_2::{Clifford, CliffordBatch, PauliBits, Phase, Support, loss::LossInfo};

use crate::{RekeyBijective, SignFlipByKey, Sum};

impl<S: Support, P> Sum<S, P> {
    fn check_qubit(&self, qubit: usize) {
        assert!(qubit < self.n_sites(), "qubit out of bounds");
    }

    fn check_pair(&self, qubit0: usize, qubit1: usize) {
        self.check_qubit(qubit0);
        self.check_qubit(qubit1);
        assert_ne!(qubit0, qubit1, "two-qubit gates require distinct qubits");
    }
}

impl<S: RekeyBijective, P> Sum<S, P> {
    fn conjugate_terms(&mut self, gate: impl Fn(&mut Phased<S::Key>)) {
        self.storage.rekey_bijective(|key, coeff| {
            let mut phased_word = Phased::new(key);
            gate(&mut phased_word);
            let (key, phase) = phased_word.into_parts();
            // Keys carry no phase; absorb the conjugation sign into the coefficient.
            match phase {
                Phase::Pos1 => (key, coeff),
                Phase::Neg1 => (key, -coeff),
                _ => unreachable!("Clifford conjugation produces only real signs"),
            }
        });
    }
}

/// Conjugate each term, validating qubit indices before changing storage.
/// Two-qubit gates also require distinct indices, including on empty sums.
impl<S, P> Clifford for Sum<S, P>
where
    S: RekeyBijective + SignFlipByKey,
    S::Key: PauliBits + LossInfo,
{
    fn x(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.storage
            .flip_sign_by_key(|key| !key.is_lost(qubit) && key.z_bit(qubit));
    }

    fn y(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.storage
            .flip_sign_by_key(|key| !key.is_lost(qubit) && (key.x_bit(qubit) ^ key.z_bit(qubit)));
    }

    fn z(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.storage
            .flip_sign_by_key(|key| !key.is_lost(qubit) && key.x_bit(qubit));
    }

    fn h(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.h(qubit));
    }

    fn s(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.s(qubit));
    }

    fn s_dag(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.s_dag(qubit));
    }

    fn sqrt_x(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.sqrt_x(qubit));
    }

    fn sqrt_x_dag(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.sqrt_x_dag(qubit));
    }

    fn sqrt_y(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.sqrt_y(qubit));
    }

    fn sqrt_y_dag(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.conjugate_terms(|word| word.sqrt_y_dag(qubit));
    }

    fn cnot(&mut self, control: usize, target: usize) {
        self.check_pair(control, target);
        self.conjugate_terms(|word| word.cnot(control, target));
    }

    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        self.check_pair(qubit0, qubit1);
        self.conjugate_terms(|word| word.cz(qubit0, qubit1));
    }

    fn cy(&mut self, control: usize, target: usize) {
        self.check_pair(control, target);
        self.conjugate_terms(|word| word.cy(control, target));
    }
}

impl<S: Support, P> CliffordBatch for Sum<S, P> where Self: Clifford {}
