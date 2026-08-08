// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use ppvm_traits_2::{Depolarizing, Depolarizing2, Pauli, PauliError, TwoQubitPauliError};
use rand::RngExt;

use crate::mixture::equality::Mutation;
use crate::mixture::fingerprint::sign_mask;
use crate::mixture::{GeneralizedTableauMixture, LazyBranch};
use crate::{Bitstring, GeneralizedTableau, RowStorage};

fn pauli_delta<A: RowStorage, I, H>(
    tab: &GeneralizedTableau<A, I, H>,
    qubit: usize,
    pauli: Pauli,
) -> u64 {
    tab.tableau.rows().enumerate().fold(0, |delta, (row, _)| {
        let site = tab.tableau.row_site(row, qubit);
        let flip = matches!(
            (pauli, site),
            (Pauli::X, Pauli::Y | Pauli::Z)
                | (Pauli::Y, Pauli::X | Pauli::Z)
                | (Pauli::Z, Pauli::X | Pauli::Y)
        );
        if flip { delta ^ sign_mask(row) } else { delta }
    })
}

impl<A, I, H> PauliError<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn pauli_error(&mut self, qubit: usize, probabilities: [f64; 3]) {
        self.rebuild_buckets();
        let original_len = self.entries.len();
        let total: f64 = probabilities.iter().sum();
        let mut branches: Vec<LazyBranch> = Vec::with_capacity(3 * original_len);
        for parent in 0..original_len {
            if self.entries[parent].0.is_lost[qubit] {
                continue;
            }
            let base = self.fingerprints[parent];
            for (pauli, probability) in [Pauli::X, Pauli::Y, Pauli::Z]
                .into_iter()
                .zip(probabilities)
            {
                let delta = pauli_delta(&self.entries[parent].0, qubit, pauli);
                branches.push((
                    parent,
                    Mutation::Pauli { pauli, qubit },
                    self.entries[parent].1 * probability,
                    base ^ delta,
                ));
            }
            self.entries[parent].1 *= 1.0 - total;
        }
        if self.insert_lazy_branches(branches) {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<A, I, H> Depolarizing<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn depolarize1(&mut self, qubit: usize, probability: f64) {
        self.pauli_error(qubit, [probability / 3.0; 3]);
    }
}

impl<A, I, H> TwoQubitPauliError<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn two_qubit_pauli_error(&mut self, qubit0: usize, qubit1: usize, probabilities: [f64; 15]) {
        const PAIRS: [(Pauli, Pauli); 15] = [
            (Pauli::I, Pauli::X),
            (Pauli::I, Pauli::Y),
            (Pauli::I, Pauli::Z),
            (Pauli::X, Pauli::I),
            (Pauli::X, Pauli::X),
            (Pauli::X, Pauli::Y),
            (Pauli::X, Pauli::Z),
            (Pauli::Y, Pauli::I),
            (Pauli::Y, Pauli::X),
            (Pauli::Y, Pauli::Y),
            (Pauli::Y, Pauli::Z),
            (Pauli::Z, Pauli::I),
            (Pauli::Z, Pauli::X),
            (Pauli::Z, Pauli::Y),
            (Pauli::Z, Pauli::Z),
        ];
        self.rebuild_buckets();
        let original_len = self.entries.len();
        let total: f64 = probabilities.iter().sum();
        let mut branches = Vec::with_capacity(15 * original_len);
        for parent in 0..original_len {
            let tab = &self.entries[parent].0;
            if tab.is_lost[qubit0] || tab.is_lost[qubit1] {
                continue;
            }
            for ((first, second), probability) in PAIRS.into_iter().zip(probabilities) {
                let _: u64 = self.rng.random();
                let delta = pauli_delta(tab, qubit0, first) ^ pauli_delta(tab, qubit1, second);
                branches.push((
                    parent,
                    Mutation::Pauli2 {
                        first,
                        second,
                        qubit0,
                        qubit1,
                    },
                    self.entries[parent].1 * probability,
                    self.fingerprints[parent] ^ delta,
                ));
            }
            self.entries[parent].1 *= 1.0 - total;
        }
        if self.insert_lazy_branches(branches) {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<A, I, H> Depolarizing2<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn depolarize2(&mut self, qubit0: usize, qubit1: usize, probability: f64) {
        self.two_qubit_pauli_error(qubit0, qubit1, [probability / 15.0; 15]);
    }
}
