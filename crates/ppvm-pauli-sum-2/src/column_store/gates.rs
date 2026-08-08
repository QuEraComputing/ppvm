// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::ColumnStore;
use super::*;

impl<K, C> ScaleByKey<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K, &mut C) + Send + Sync,
    {
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.keys.get(i);
            f(&key, &mut self.primary.coeffs[i]);
        }
    }
}

impl<K, C> SignFlipByKey<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8 + Send + Sync,
    {
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.keys.get(i);
            let sign = f(&key);
            if sign != 1 {
                let flipped = self.primary.coeffs[i].mul_sign(sign);
                self.primary.coeffs[i] = flipped;
            }
        }
    }
}

impl<K, C> RekeyBijective<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C) + Send + Sync,
    {
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.keys.get(i);
            let coeff = std::mem::replace(&mut self.primary.coeffs[i], C::zero());
            let (new_key, new_coeff) = f(key, coeff);
            self.primary.hashes[i] = new_key.key_hash();
            self.primary.keys.set(i, new_key);
            self.primary.coeffs[i] = new_coeff;
        }
        self.primary.reindex();
    }
}

impl<K, C> BranchInPlace<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    #[inline]
    fn branch_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C, &mut Vec<(K, C)>) + Send + Sync,
    {
        self.scratch.clear();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.keys.get(i);
            f(&key, &mut self.primary.coeffs[i], &mut self.scratch);
        }
        for (key, _) in self.scratch.iter() {
            let _ = key.key_hash();
        }
        self.primary
            .reserve_for_live_len(self.primary.len() + self.scratch.len());
        for (key, coeff) in self.scratch.drain(..) {
            self.primary.add(key, coeff);
        }
    }
}

impl<K, C> MultiplyInPlace<K, C> for ColumnStore<K, C>
where
    K: Columnar + KeyProduct,
    C: ImaginaryUnit,
{
    fn multiply_in_place(&mut self, other: &Self) {
        self.aux.clear();
        self.aux.reserve_for_live_len(product_capacity_hint(
            self.primary.len(),
            other.primary.len(),
        ));
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let p = self.primary.key(i);
            for j in 0..other.primary.rows() {
                if !other.primary.is_live(j) {
                    continue;
                }
                let q = other.primary.key(j);
                let (k, phase) = p.key_mul(&q);
                let c = phase
                    .apply(&(self.primary.coeffs[i].clone() * other.primary.coeffs[j].clone()));
                self.aux.add(k, c);
            }
        }
        std::mem::swap(&mut self.primary, &mut self.aux);
        self.aux.clear();
    }
}
