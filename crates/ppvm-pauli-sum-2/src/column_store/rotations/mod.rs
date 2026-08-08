// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::ColumnStore;
use super::*;

mod rx;

impl<K, C> RotateInPlace<K, C> for ColumnStore<K, C>
where
    K: Columnar,
    C: Coefficient,
{
    /// Scale all diagonals before merging any branch.
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)> + Send + Sync,
    {
        self.scratch.clear();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let key = self.primary.keys.get(i);
            if let Some(term) = f(&key, &mut self.primary.coeffs[i]) {
                let _ = term.0.key_hash();
                self.scratch.push(term);
            }
        }
        self.primary
            .reserve_for_live_len(self.primary.len() + self.scratch.len());
        for (key, coeff) in self.scratch.drain(..) {
            self.primary.add_likely_present(key, coeff);
        }
    }

    #[inline(always)]
    fn rotate_x(&mut self, qubit: usize, sin: C, cos: C)
    where
        K: PauliBits,
    {
        if self.primary.is_dense() {
            self.rotate_x_kernel::<false>(qubit, sin, cos);
        } else {
            self.primary.ensure_sparse_cache();
            self.rotate_x_kernel::<true>(qubit, sin, cos);
        }
    }

    #[inline]
    fn rotate_zz(&mut self, a: usize, b: usize, sin: C, cos: C)
    where
        K: PauliBits + Clone,
    {
        self.scratch.clear();
        for i in 0..self.primary.rows() {
            if !self.primary.is_live(i) {
                continue;
            }
            let lost_a = self.primary.keys.is_lost(i, a);
            let lost_b = self.primary.keys.is_lost(i, b);
            let (key, sign) = if lost_a || lost_b {
                let site = if lost_a { b } else { a };
                if (lost_a && lost_b)
                    || self.primary.keys.is_lost(i, site)
                    || !self.primary.keys.x_bit(i, site)
                {
                    continue;
                }
                (
                    self.primary.keys.toggled_bits(i, site, false, true),
                    if self.primary.keys.z_bit(i, site) {
                        1
                    } else {
                        -1
                    },
                )
            } else {
                let xa = self.primary.keys.x_bit(i, a);
                let xb = self.primary.keys.x_bit(i, b);
                if xa == xb {
                    continue;
                }
                (
                    self.primary
                        .keys
                        .toggled_bits2(i, a, false, true, b, false, true),
                    if if xa {
                        self.primary.keys.z_bit(i, a)
                    } else {
                        self.primary.keys.z_bit(i, b)
                    } {
                        1
                    } else {
                        -1
                    },
                )
            };
            let branch = self.primary.coeffs[i].clone() * sin.mul_sign(sign);
            self.primary.coeffs[i] *= cos.clone();
            let _ = key.key_hash();
            self.scratch.push((key, branch));
        }
        self.primary
            .reserve_for_live_len(self.primary.len() + self.scratch.len());
        for (key, coeff) in self.scratch.drain(..) {
            self.primary.add_likely_present(key, coeff);
        }
    }
}
