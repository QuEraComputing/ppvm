// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Coefficient, ImaginaryUnit, Indexable, KeyProduct};

use super::IndexMapStore;
use crate::store::{BranchInPlace, MultiplyInPlace, product_capacity_hint};

impl<K: Indexable, C: Coefficient> BranchInPlace<K, C> for IndexMapStore<K, C> {
    fn branch_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C, &mut Vec<(K, C)>) + Send + Sync,
    {
        self.scratch.clear();
        for (key, coeff) in &mut self.primary {
            f(key, coeff, &mut self.scratch);
        }
        for (key, _) in &self.scratch {
            let _ = key.key_hash();
        }

        // Build the deduplicated branch map first. Legacy `consume` compares map
        // cardinalities, not the raw fan-out count; duplicate branch keys must
        // not force branch-first ordering.
        self.aux.clear();
        self.aux.reserve(self.scratch.len());
        for (key, coeff) in self.scratch.drain(..) {
            self.aux
                .entry(key)
                .and_modify(|value| *value += coeff.clone())
                .or_insert(coeff);
        }

        // Merge the smaller map into the larger. The chosen destination defines
        // observable insertion order.
        if self.aux.len() > self.primary.len() {
            self.aux.reserve(self.primary.len());
            for (key, coeff) in self.primary.drain(..) {
                self.aux
                    .entry(key)
                    .and_modify(|value| *value += coeff.clone())
                    .or_insert(coeff);
            }
            std::mem::swap(&mut self.primary, &mut self.aux);
            self.aux.clear();
        } else {
            for (key, coeff) in self.aux.drain(..) {
                self.primary
                    .entry(key)
                    .and_modify(|value| *value += coeff.clone())
                    .or_insert(coeff);
            }
        }
    }
}

impl<K, C> MultiplyInPlace<K, C> for IndexMapStore<K, C>
where
    K: Indexable + KeyProduct,
    C: ImaginaryUnit,
{
    fn multiply_in_place(&mut self, other: &Self) {
        self.aux.clear();
        self.aux.reserve(product_capacity_hint(
            self.primary.len(),
            other.primary.len(),
        ));
        for (p, a) in &self.primary {
            for (q, b) in &other.primary {
                let (key, phase) = p.key_mul(q);
                let coeff = phase.apply(&(a.clone() * b.clone()));
                self.aux
                    .entry(key)
                    .and_modify(|value| *value += coeff.clone())
                    .or_insert(coeff);
            }
        }
        std::mem::swap(&mut self.primary, &mut self.aux);
        self.aux.clear();
    }
}
