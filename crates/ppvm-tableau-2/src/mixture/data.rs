// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashMap;
use std::hash::BuildHasher;

use fxhash::FxBuildHasher;
use ppvm_pauli_word_2::DefaultStorage;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use super::equality::{Mutation, apply_mutation, structurally_equal, structurally_equal_mutated};
use super::fingerprint::fingerprint;
use crate::{Bitstring, GeneralizedTableau, RowStorage};

pub(crate) type Branch<A, I, H> = (GeneralizedTableau<A, I, H>, f64, u64);
pub(crate) type LazyBranch = (usize, Mutation, f64, u64);

/// A classical probability distribution over complete generalized-tableau states.
///
/// Fingerprints only select collision buckets. Entries merge only after a full
/// frame/loss comparison and coefficient-wise approximate comparison.
#[derive(Clone)]
pub struct GeneralizedTableauMixture<
    A: RowStorage = DefaultStorage,
    I: Bitstring = usize,
    H = FxBuildHasher,
> {
    pub n_qubits: usize,
    pub entries: Vec<(GeneralizedTableau<A, I, H>, f64)>,
    pub(crate) rng: SmallRng,
    pub(crate) sum_cutoff: f64,
    pub(crate) fingerprints: Vec<u64>,
    pub(crate) buckets: HashMap<u64, Vec<usize>, H>,
    pub(crate) dirty: bool,
}

impl<A, I, H> GeneralizedTableauMixture<A, I, H>
where
    A: RowStorage,
    I: Bitstring,
    H: BuildHasher + Clone + Default,
{
    pub fn new(n_qubits: usize, coefficient_threshold: f64, sum_cutoff: f64) -> Self {
        let mut rng: SmallRng = rand::make_rng();
        let seed = rng.random();
        Self::from_seeds(n_qubits, coefficient_threshold, sum_cutoff, rng, seed)
    }

    pub fn new_with_seed(
        n_qubits: usize,
        coefficient_threshold: f64,
        sum_cutoff: f64,
        seed: u64,
    ) -> Self {
        let mut rng = SmallRng::seed_from_u64(seed);
        let tab_seed = rng.random();
        Self::from_seeds(n_qubits, coefficient_threshold, sum_cutoff, rng, tab_seed)
    }

    fn from_seeds(
        n_qubits: usize,
        coefficient_threshold: f64,
        sum_cutoff: f64,
        rng: SmallRng,
        tab_seed: u64,
    ) -> Self {
        let tab = GeneralizedTableau::new_with_seed(n_qubits, coefficient_threshold, tab_seed);
        let fp = fingerprint(&tab);
        let mut mixture = Self {
            n_qubits,
            entries: Vec::new(),
            rng,
            sum_cutoff,
            fingerprints: Vec::new(),
            buckets: HashMap::with_capacity_and_hasher(1, H::default()),
            dirty: false,
        };
        // Construction obeys the same strict `probability > sum_cutoff` rule as
        // every later branch insertion. Old therefore starts empty at
        // `sum_cutoff >= 1.0`; bypassing this door silently changed that boundary.
        let _ = mixture.insert_branches(vec![(tab, 1.0, fp)]);
        mixture
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&GeneralizedTableau<A, I, H>, &f64)> {
        self.entries
            .iter()
            .map(|(tab, probability)| (tab, probability))
    }

    pub fn normalize_probabilities(&mut self) {
        let norm: f64 = self.entries.iter().map(|entry| entry.1).sum();
        for (_, probability) in &mut self.entries {
            *probability /= norm;
        }
    }

    pub fn truncate(&mut self) {
        let before = self.entries.len();
        self.entries
            .retain(|(_, probability)| *probability > self.sum_cutoff);
        if self.entries.len() != before {
            self.normalize_probabilities();
            self.dirty = true;
        }
    }

    pub(crate) fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    pub(crate) fn rebuild_buckets(&mut self) {
        if !self.dirty {
            return;
        }
        self.fingerprints.clear();
        self.buckets.clear();
        for (index, (tab, _)) in self.entries.iter().enumerate() {
            let fp = fingerprint(tab);
            self.fingerprints.push(fp);
            self.buckets.entry(fp).or_default().push(index);
        }
        self.dirty = false;
    }

    pub(crate) fn insert_branches(&mut self, branches: Vec<Branch<A, I, H>>) -> bool {
        self.rebuild_buckets();
        let mut dropped = false;
        for (tab, probability, fp) in branches {
            let found = self.buckets.get(&fp).and_then(|candidates| {
                candidates
                    .iter()
                    .copied()
                    .find(|&i| structurally_equal(&self.entries[i].0, &tab))
            });
            if let Some(index) = found {
                self.entries[index].1 += probability;
            } else if probability > self.sum_cutoff {
                let index = self.entries.len();
                self.entries.push((tab, probability));
                self.fingerprints.push(fp);
                self.buckets.entry(fp).or_default().push(index);
            } else {
                dropped = true;
            }
        }
        dropped
    }

    pub(crate) fn insert_lazy_branches(&mut self, branches: Vec<LazyBranch>) -> bool {
        self.rebuild_buckets();
        let parent_count = self.entries.len();
        let mut dropped = false;
        for (parent, mutation, probability, fp) in branches {
            debug_assert!(parent < parent_count);
            let found = self.buckets.get(&fp).and_then(|candidates| {
                candidates.iter().copied().find(|&i| {
                    structurally_equal_mutated(
                        &self.entries[i].0,
                        &self.entries[parent].0,
                        mutation,
                    )
                })
            });
            if let Some(index) = found {
                self.entries[index].1 += probability;
            } else if probability > self.sum_cutoff {
                let mut tab = self.entries[parent].0.clone();
                apply_mutation(&mut tab, mutation);
                let index = self.entries.len();
                self.entries.push((tab, probability));
                self.fingerprints.push(fp);
                self.buckets.entry(fp).or_default().push(index);
            } else {
                dropped = true;
            }
        }
        dropped
    }
}

/// Compatibility spelling retained for old Rust and future Python adapters.
pub type GeneralizedTableauSum<A = DefaultStorage, I = usize, H = FxBuildHasher> =
    GeneralizedTableauMixture<A, I, H>;
