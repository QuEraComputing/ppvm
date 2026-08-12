// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{
    Accumulate, Coefficient, Conjugate, ImaginaryUnit, Indexable, KeyBatch, KeyProduct, Multiply,
    Pair, Retain, Scale, Support, TermBatch,
};

use super::IndexMapStore;
use crate::store::{AddTerm, product_capacity_hint};

impl<K: Indexable, C: Coefficient> Support for IndexMapStore<K, C> {
    type Key = K;
    type Coeff = C;

    fn len(&self) -> usize {
        self.primary.len()
    }

    fn get(&self, key: &K) -> Option<C> {
        self.primary.get(key).cloned()
    }

    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        self.primary.iter().map(|(k, c)| (k.clone(), c.clone()))
    }

    fn for_each_ref(&self, mut f: impl FnMut(&K, &C)) {
        for (k, c) in &self.primary {
            f(k, c);
        }
    }
}

impl<K: Indexable, C: Coefficient> Accumulate for IndexMapStore<K, C> {
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        for (k, c) in terms.iter() {
            self.primary
                .entry(k.clone())
                .and_modify(|v| *v += c.clone())
                .or_insert_with(|| c.clone());
        }
    }

    fn reduce(&mut self) {
        self.primary.retain(|_, c| !c.is_zero());
    }
}

impl<K: Indexable, C: Coefficient> Scale for IndexMapStore<K, C> {
    fn scale(&mut self, factor: &C) {
        for c in self.primary.values_mut() {
            *c *= factor.clone();
        }
    }
}

impl<K: Indexable, C: Coefficient> Pair for IndexMapStore<K, C> {
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        debug_assert!(out.len() >= keys.keys().len());
        for (slot, key) in out.iter_mut().zip(keys.keys()) {
            *slot = self.primary.get(key).cloned();
        }
    }

    fn overlap(&self, other: &Self) -> C {
        if self.len() <= other.len() {
            self.primary
                .iter()
                .filter_map(|(k, a)| other.primary.get(k).map(|b| a.clone() * b.clone()))
                .sum()
        } else {
            other
                .primary
                .iter()
                .filter_map(|(k, b)| self.primary.get(k).map(|a| a.clone() * b.clone()))
                .sum()
        }
    }

    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        if self.len() <= other.len() {
            self.primary
                .iter()
                .filter_map(|(k, a)| other.primary.get(k).map(|b| a.conj() * b.clone()))
                .sum()
        } else {
            other
                .primary
                .iter()
                .filter_map(|(k, b)| self.primary.get(k).map(|a| a.conj() * b.clone()))
                .sum()
        }
    }
}

impl<K: Indexable, C: Coefficient> Retain<K, C> for IndexMapStore<K, C> {
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        self.primary.retain(|k, c| keep(k, c));
    }
}

impl<K, C> Multiply for IndexMapStore<K, C>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    fn multiply_into(&self, other: &Self, acc: &mut Self) {
        acc.primary
            .reserve(product_capacity_hint(self.len(), other.len()));
        for (p, a) in &self.primary {
            for (q, b) in &other.primary {
                let (key, phase) = p.key_mul(q);
                let coeff = phase.apply(&(a.clone() * b.clone()));
                AddTerm::add_term(acc, key, coeff);
            }
        }
    }
}
