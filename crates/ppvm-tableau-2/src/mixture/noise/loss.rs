// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use ppvm_traits_2::{CorrelatedLossChannel, LossChannel, ResetLossChannel};
use rand::RngExt;

use crate::mixture::equality::Mutation;
use crate::mixture::fingerprint::loss_mask;
use crate::mixture::{Branch, GeneralizedTableauMixture, LazyBranch};
use crate::{Bitstring, RowStorage};

impl<A, I, H> LossChannel<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn loss_channel(&mut self, qubit: usize, probability: f64) {
        self.rebuild_buckets();
        let original_len = self.entries.len();
        let mut branches: Vec<LazyBranch> = Vec::with_capacity(original_len);
        for parent in 0..original_len {
            if self.entries[parent].0.is_lost[qubit] {
                continue;
            }
            branches.push((
                parent,
                Mutation::Loss { qubit },
                self.entries[parent].1 * probability,
                self.fingerprints[parent] ^ loss_mask(qubit),
            ));
            self.entries[parent].1 *= 1.0 - probability;
        }
        if self.insert_lazy_branches(branches) {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<A, I, H> CorrelatedLossChannel<f64> for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn correlated_loss_channel(&mut self, qubit0: usize, qubit1: usize, probabilities: [f64; 3]) {
        self.rebuild_buckets();
        let original_len = self.entries.len();
        let mut branches = Vec::with_capacity(3 * original_len);
        for parent in 0..original_len {
            let lost0 = self.entries[parent].0.is_lost[qubit0];
            let lost1 = self.entries[parent].0.is_lost[qubit1];
            if lost0 || lost1 {
                let qubit = if lost0 { qubit1 } else { qubit0 };
                if !self.entries[parent].0.is_lost[qubit] {
                    let _: u64 = self.rng.random();
                    push_loss(self, &mut branches, parent, qubit, probabilities[2]);
                }
                continue;
            }
            let _: [u64; 3] = std::array::from_fn(|_| self.rng.random());
            let base = self.fingerprints[parent];
            let p = self.entries[parent].1;
            branches.push((
                parent,
                Mutation::Loss2 { qubit0, qubit1 },
                p * probabilities[0],
                base ^ loss_mask(qubit0) ^ loss_mask(qubit1),
            ));
            branches.push((
                parent,
                Mutation::Loss { qubit: qubit0 },
                p * probabilities[1] / 2.0,
                base ^ loss_mask(qubit0),
            ));
            branches.push((
                parent,
                Mutation::Loss { qubit: qubit1 },
                p * probabilities[1] / 2.0,
                base ^ loss_mask(qubit1),
            ));
            self.entries[parent].1 *= 1.0 - probabilities[0] - probabilities[1];
        }
        if self.insert_lazy_branches(branches) {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

fn push_loss<A, I, H>(
    mixture: &mut GeneralizedTableauMixture<A, I, H>,
    branches: &mut Vec<LazyBranch>,
    parent: usize,
    qubit: usize,
    probability: f64,
) where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    branches.push((
        parent,
        Mutation::Loss { qubit },
        mixture.entries[parent].1 * probability,
        mixture.fingerprints[parent] ^ loss_mask(qubit),
    ));
    mixture.entries[parent].1 *= 1.0 - probability;
}

impl<A, I, H> ResetLossChannel for GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    fn reset_loss_channel(&mut self, qubit: usize) {
        self.rebuild_buckets();
        let mut indices: Vec<_> = (0..self.entries.len())
            .filter(|&i| self.entries[i].0.is_lost[qubit])
            .collect();
        indices.reverse();
        let mut branches: Vec<Branch<A, I, H>> = Vec::with_capacity(indices.len());
        for index in indices {
            let (mut tab, probability) = self.entries.swap_remove(index);
            let fp = self.fingerprints.swap_remove(index) ^ loss_mask(qubit);
            tab.is_lost[qubit] = false;
            branches.push((tab, probability, fp));
        }
        self.mark_dirty();
        let _ = self.insert_branches(branches);
    }
}
