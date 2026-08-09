// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use ppvm_traits_2::{Depolarizing, Depolarizing2, Pauli, PauliError, TwoQubitPauliError};

use crate::mixture::equality::Mutation;
use crate::mixture::fingerprint::sign_mask;
use crate::mixture::{GeneralizedTableauMixture, LazyBranch};
use crate::{Bitstring, GeneralizedTableau, RowStorage};

fn pauli_deltas<A: RowStorage, I, H>(tab: &GeneralizedTableau<A, I, H>, qubit: usize) -> [u64; 4] {
    let (mut dx, mut dz) = (0, 0);
    for (row, _) in tab.tableau.rows().enumerate() {
        let site = tab.tableau.row_site(row, qubit);
        let mask = sign_mask(row);
        if matches!(site, Pauli::Y | Pauli::Z) {
            dx ^= mask;
        }
        if matches!(site, Pauli::X | Pauli::Y) {
            dz ^= mask;
        }
    }
    [0, dx, dx ^ dz, dz]
}

fn pauli_index(pauli: Pauli) -> usize {
    match pauli {
        Pauli::I => 0,
        Pauli::X => 1,
        Pauli::Y => 2,
        Pauli::Z => 3,
    }
}

fn two_qubit_pauli_deltas<A: RowStorage, I, H>(
    tab: &GeneralizedTableau<A, I, H>,
    qubit0: usize,
    qubit1: usize,
) -> [[u64; 4]; 2] {
    let ([mut dx0, mut dz0], [mut dx1, mut dz1]) = ([0, 0], [0, 0]);
    for (row, _) in tab.tableau.rows().enumerate() {
        let first = tab.tableau.row_site(row, qubit0);
        let second = tab.tableau.row_site(row, qubit1);
        let mask = sign_mask(row);
        if matches!(first, Pauli::Y | Pauli::Z) {
            dx0 ^= mask;
        }
        if matches!(first, Pauli::X | Pauli::Y) {
            dz0 ^= mask;
        }
        if matches!(second, Pauli::Y | Pauli::Z) {
            dx1 ^= mask;
        }
        if matches!(second, Pauli::X | Pauli::Y) {
            dz1 ^= mask;
        }
    }
    [[0, dx0, dx0 ^ dz0, dz0], [0, dx1, dx1 ^ dz1, dz1]]
}

impl<A, I, H> PauliError<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit: usize,
        probabilities: [f64; 3],
        _rng: &mut R,
    ) {
        self.rebuild_buckets();
        let original_len = self.entries.len();
        let total: f64 = probabilities.iter().sum();
        let mut branches: Vec<LazyBranch> = Vec::with_capacity(3 * original_len);
        for parent in 0..original_len {
            if self.entries[parent].0.is_lost[qubit] {
                continue;
            }
            let base = self.fingerprints[parent];
            let deltas = pauli_deltas(&self.entries[parent].0, qubit);
            for (pauli, probability) in [Pauli::X, Pauli::Y, Pauli::Z]
                .into_iter()
                .zip(probabilities)
            {
                branches.push((
                    parent,
                    Mutation::Pauli { pauli, qubit },
                    self.entries[parent].1 * probability,
                    base ^ deltas[pauli_index(pauli)],
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
    fn depolarize1<R: rand::Rng + ?Sized>(&mut self, qubit: usize, probability: f64, rng: &mut R) {
        self.pauli_error(qubit, [probability / 3.0; 3], rng);
    }
}

impl<A, I, H> TwoQubitPauliError<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn two_qubit_pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        qubit0: usize,
        qubit1: usize,
        probabilities: [f64; 15],
        _rng: &mut R,
    ) {
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
            if self.entries[parent].0.is_lost[qubit0] || self.entries[parent].0.is_lost[qubit1] {
                continue;
            }
            self.burn_legacy_tableau_seeds(15);
            let tab = &self.entries[parent].0;
            let deltas = two_qubit_pauli_deltas(tab, qubit0, qubit1);
            for ((first, second), probability) in PAIRS.into_iter().zip(probabilities) {
                let delta = deltas[0][pauli_index(first)] ^ deltas[1][pauli_index(second)];
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
    fn depolarize2<R: rand::Rng + ?Sized>(
        &mut self,
        qubit0: usize,
        qubit1: usize,
        probability: f64,
        rng: &mut R,
    ) {
        self.two_qubit_pauli_error(qubit0, qubit1, [probability / 15.0; 15], rng);
    }
}
