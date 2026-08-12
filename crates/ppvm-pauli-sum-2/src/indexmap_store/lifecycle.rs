// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use indexmap::IndexMap;
use ppvm_traits_2::{Coefficient, IdentityBuildHasher, Indexable, TermProducer};

use super::IndexMapStore;
use crate::store::{AddTerm, ApplyProducer, InsertTerm, StoreAlloc};

impl<K, C> StoreAlloc for IndexMapStore<K, C> {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            primary: IndexMap::with_capacity_and_hasher(capacity, IdentityBuildHasher),
            aux: IndexMap::with_capacity_and_hasher(capacity, IdentityBuildHasher),
            scratch: Vec::with_capacity(capacity),
            batch: ppvm_traits_2::TermBatch::new(),
        }
    }

    fn reset(&mut self) {
        self.primary.clear();
    }
}

impl<K: Indexable, C: Coefficient> AddTerm<K, C> for IndexMapStore<K, C> {
    fn add_term(&mut self, key: K, coeff: C) {
        self.primary
            .entry(key)
            .and_modify(|value| *value += coeff.clone())
            .or_insert(coeff);
    }
}

impl<K: Indexable, C: Coefficient> InsertTerm<K, C> for IndexMapStore<K, C> {
    fn insert_term(&mut self, key: K, coeff: C) {
        // IndexMap replacement retains the original insertion position.
        self.primary.insert(key, coeff);
    }
}

impl<K: Indexable, C: Coefficient> ApplyProducer<K, C> for IndexMapStore<K, C> {
    fn apply_producer<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<K, C>,
    {
        self.batch.clear();
        for (key, coeff) in &self.primary {
            producer.produce(key, coeff, &mut self.batch);
        }
        self.primary.clear();
        for (key, coeff) in self.batch.iter() {
            self.primary
                .entry(key.clone())
                .and_modify(|value| *value += coeff.clone())
                .or_insert_with(|| coeff.clone());
        }
        self.batch.clear();
    }
}
