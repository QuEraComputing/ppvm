// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Coefficient, Indexable};

use super::IndexMapStore;
use crate::store::{RekeyBijective, RotateInPlace, ScaleByKey, SignFlipByKey};

impl<K: Indexable, C: Coefficient> ScaleByKey<K, C> for IndexMapStore<K, C> {
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C) + Send + Sync,
    {
        for (key, coeff) in &mut self.primary {
            f(key, coeff);
        }
    }
}

impl<K: Indexable, C: Coefficient> SignFlipByKey<K, C> for IndexMapStore<K, C> {
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8 + Send + Sync,
    {
        for (key, coeff) in &mut self.primary {
            let sign = f(key);
            if sign != 1 {
                *coeff = coeff.mul_sign(sign);
            }
        }
    }
}

impl<K: Indexable, C: Coefficient> RekeyBijective<K, C> for IndexMapStore<K, C> {
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C) + Send + Sync,
    {
        self.aux.clear();
        self.aux.reserve(self.primary.len());
        for (key, coeff) in self.primary.drain(..) {
            let (new_key, new_coeff) = f(key, coeff);
            let displaced = self.aux.insert(new_key, new_coeff);
            debug_assert!(
                displaced.is_none(),
                "RekeyBijective requires an injective re-key"
            );
        }
        std::mem::swap(&mut self.primary, &mut self.aux);
    }
}

impl<K: Indexable, C: Coefficient> RotateInPlace<K, C> for IndexMapStore<K, C> {
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)> + Send + Sync,
    {
        self.scratch.clear();
        for (key, coeff) in &mut self.primary {
            if let Some(term) = f(key, coeff) {
                let _ = term.0.key_hash();
                self.scratch.push(term);
            }
        }
        for (key, coeff) in self.scratch.drain(..) {
            self.primary
                .entry(key)
                .and_modify(|value| *value += coeff.clone())
                .or_insert(coeff);
        }
    }
}
