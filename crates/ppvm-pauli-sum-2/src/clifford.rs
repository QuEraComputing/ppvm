// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_phased_pauli_word_2::Phased;
use ppvm_traits_2::{Clifford, CliffordBatch, Coefficient, PauliBits, Phase, loss::LossInfo};

use crate::{RekeyBijective, Sum, UpdateCoefficients};

/// Whether a gate on `qubit` can change `key`: every Clifford fixes `I` and lost sites.
fn active<K: PauliBits + LossInfo>(key: &K, qubit: usize) -> bool {
    !key.is_lost(qubit) && (key.x_bit(qubit) || key.z_bit(qubit))
}

impl<S: RekeyBijective, P> Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
{
    /// Conjugate the terms for which `touched` is true; all others are left in place.
    fn conjugate_terms(
        &mut self,
        touched: impl Fn(&S::Key) -> bool,
        gate: impl Fn(&mut Phased<S::Key>),
    ) {
        self.storage.rekey_bijective(|key, coeff| {
            if !touched(key) {
                return None;
            }
            let mut phased_word = Phased::new(key.clone());
            gate(&mut phased_word);
            let (new_key, phase) = phased_word.into_parts();
            // Keys carry no phase; absorb the conjugation sign into the coefficient.
            match phase {
                Phase::Pos1 => {}
                Phase::Neg1 => coeff.mul_sign_assign(-1),
                _ => unreachable!("Clifford conjugation produces only real signs"),
            }
            (new_key != *key).then_some(new_key)
        });
    }

    fn conjugate_single(&mut self, qubit: usize, gate: impl Fn(&mut Phased<S::Key>)) {
        self.check_qubit(qubit);
        self.conjugate_terms(|key| active(key, qubit), gate);
    }

    fn conjugate_pair(&mut self, a: usize, b: usize, gate: impl Fn(&mut Phased<S::Key>)) {
        self.check_pair(a, b);
        self.conjugate_terms(|key| active(key, a) || active(key, b), gate);
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
            &[qubit] => self.conjugate_single(qubit, |word| gate(word, qubit)),
            _ => self.conjugate_terms(
                |key| indices.iter().any(|&q| active(key, q)),
                |word| {
                    for &qubit in indices {
                        gate(word, qubit);
                    }
                },
            ),
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
            &[(a, b)] => self.conjugate_pair(a, b, |word| gate(word, a, b)),
            _ => self.conjugate_terms(
                |key| pairs.iter().any(|&(a, b)| active(key, a) || active(key, b)),
                |word| {
                    for &(a, b) in pairs {
                        gate(word, a, b);
                    }
                },
            ),
        }
    }
}

impl<S: UpdateCoefficients, P> Sum<S, P> {
    fn flip_signs(&mut self, negate: impl Fn(&S::Key) -> bool) {
        self.storage.update_coefficients(|key, coeff| {
            if negate(key) {
                coeff.mul_sign_assign(-1);
            }
        });
    }

    fn flip_signs_many(&mut self, indices: &[usize], negate: impl Fn(&S::Key, usize) -> bool) {
        if indices.is_empty() {
            return;
        }
        for &qubit in indices {
            self.check_qubit(qubit);
        }
        self.flip_signs(|key| indices.iter().fold(false, |flip, &q| flip ^ negate(key, q)));
    }
}

/// Conjugate each term, validating qubit indices before changing storage.
/// Two-qubit gates also require distinct indices, including on empty sums.
impl<S, P> Clifford for Sum<S, P>
where
    S: RekeyBijective + UpdateCoefficients,
    S::Key: PauliBits + LossInfo,
    Phased<S::Key>: Clifford,
{
    fn x(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.flip_signs(|key| !key.is_lost(qubit) && key.z_bit(qubit));
    }

    fn y(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.flip_signs(|key| !key.is_lost(qubit) && (key.x_bit(qubit) ^ key.z_bit(qubit)));
    }

    fn z(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.flip_signs(|key| !key.is_lost(qubit) && key.x_bit(qubit));
    }

    fn h(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.h(qubit));
    }

    fn s(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.s(qubit));
    }

    fn s_dag(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.s_dag(qubit));
    }

    fn sqrt_x(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.sqrt_x(qubit));
    }

    fn sqrt_x_dag(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.sqrt_x_dag(qubit));
    }

    fn sqrt_y(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.sqrt_y(qubit));
    }

    fn sqrt_y_dag(&mut self, qubit: usize) {
        self.conjugate_single(qubit, |word| word.sqrt_y_dag(qubit));
    }

    fn cnot(&mut self, control: usize, target: usize) {
        self.conjugate_pair(control, target, |word| word.cnot(control, target));
    }

    fn cz(&mut self, qubit0: usize, qubit1: usize) {
        self.conjugate_pair(qubit0, qubit1, |word| word.cz(qubit0, qubit1));
    }

    fn cy(&mut self, control: usize, target: usize) {
        self.conjugate_pair(control, target, |word| word.cy(control, target));
    }
}

/// Apply gates in input order, rebuilding storage once per batch of changed keys.
/// Pauli batches only update signs. Validate the whole batch before changing terms.
impl<S, P> CliffordBatch for Sum<S, P>
where
    S: RekeyBijective + UpdateCoefficients,
    S::Key: PauliBits + LossInfo,
    Phased<S::Key>: Clifford,
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
