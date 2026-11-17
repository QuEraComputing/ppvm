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

    fn conjugate_single_batch(
        &mut self,
        indices: &[usize],
        gate: impl Fn(&mut Phased<S::Key>, usize),
    ) {
        for &qubit in indices {
            self.check_qubit(qubit);
        }
        match indices {
            [] => {}
            &[qubit] => self.conjugate_terms(|word| gate(word, qubit)),
            _ => self.conjugate_terms(|word| {
                for &qubit in indices {
                    gate(word, qubit);
                }
            }),
        }
    }

    fn conjugate_pair_batch(
        &mut self,
        pairs: &[(usize, usize)],
        gate: impl Fn(&mut Phased<S::Key>, usize, usize),
    ) {
        for &(a, b) in pairs {
            self.check_pair(a, b);
        }
        match pairs {
            [] => {}
            &[(a, b)] => self.conjugate_terms(|word| gate(word, a, b)),
            _ => self.conjugate_terms(|word| {
                for &(a, b) in pairs {
                    gate(word, a, b);
                }
            }),
        }
    }
}

impl<S: SignFlipByKey, P> Sum<S, P> {
    fn flip_signs_many(&mut self, indices: &[usize], negate: impl Fn(&S::Key, usize) -> bool) {
        if indices.is_empty() {
            return;
        }
        for &qubit in indices {
            self.check_qubit(qubit);
        }
        self.storage
            .flip_sign_by_key(|key| indices.iter().fold(false, |flip, &q| flip ^ negate(key, q)));
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

/// Apply gates in input order, rebuilding storage once per batch of changed keys.
/// Pauli batches only update signs. Validate the whole batch before changing terms.
impl<S, P> CliffordBatch for Sum<S, P>
where
    S: RekeyBijective + SignFlipByKey,
    S::Key: PauliBits + LossInfo,
{
    fn x_many(&mut self, indices: &[usize]) {
        if let &[qubit] = indices {
            self.x(qubit);
            return;
        }
        self.flip_signs_many(indices, |key, q| !key.is_lost(q) && key.z_bit(q));
    }

    fn y_many(&mut self, indices: &[usize]) {
        if let &[qubit] = indices {
            self.y(qubit);
            return;
        }
        self.flip_signs_many(indices, |key, q| {
            !key.is_lost(q) && (key.x_bit(q) ^ key.z_bit(q))
        });
    }

    fn z_many(&mut self, indices: &[usize]) {
        if let &[qubit] = indices {
            self.z(qubit);
            return;
        }
        self.flip_signs_many(indices, |key, q| !key.is_lost(q) && key.x_bit(q));
    }

    fn h_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::h);
    }

    fn s_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::s);
    }

    fn s_dag_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::s_dag);
    }

    fn sqrt_x_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::sqrt_x);
    }

    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::sqrt_x_dag);
    }

    fn sqrt_y_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::sqrt_y);
    }

    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        self.conjugate_single_batch(indices, Phased::sqrt_y_dag);
    }

    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        self.conjugate_pair_batch(pairs, Phased::cnot);
    }

    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        self.conjugate_pair_batch(pairs, Phased::cz);
    }

    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        self.conjugate_pair_batch(pairs, Phased::cy);
    }
}
