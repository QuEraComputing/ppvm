// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;

impl<K, C> ColumnStore<K, C>
where
    K: Columnar + PauliBits,
    C: Coefficient,
{
    #[inline(always)]
    pub(super) fn rotate_x_kernel<const SPARSE: bool>(&mut self, qubit: usize, sin: C, cos: C) {
        self.scratch.clear();
        let original_rows = self.primary.rows();
        let scan_len = if SPARSE {
            self.primary.live_len
        } else {
            original_rows
        };
        if original_rows <= 512 {
            if SPARSE {
                for run in 0..self.primary.live_runs.len() {
                    let (start, end) = self.primary.live_runs[run];
                    for i in start as usize..end as usize {
                        self.rotate_x_small_row::<true>(i, qubit, &sin, &cos);
                    }
                }
            } else {
                for i in 0..original_rows {
                    self.rotate_x_small_row::<false>(i, qubit, &sin, &cos);
                }
            }
            self.merge_rotation_scratch::<SPARSE>();
            return;
        }

        let closed_support = (0..scan_len)
            .map(|position| {
                if SPARSE {
                    self.primary.sparse_rows[position] as usize
                } else {
                    position
                }
            })
            .find(|&i| !self.primary.keys.is_lost(i, qubit) && self.primary.keys.z_bit(i, qubit))
            .is_none_or(|i| {
                let key = self.primary.keys.toggled_bits(i, qubit, true, false);
                if SPARSE {
                    self.primary.find(&key, key.key_hash()).is_some()
                } else {
                    self.primary.find_any(&key, key.key_hash()).is_some()
                }
            });

        if !closed_support {
            for position in 0..scan_len {
                let i = if SPARSE {
                    self.primary.sparse_rows[position] as usize
                } else {
                    position
                };
                if self.primary.keys.is_lost(i, qubit) || !self.primary.keys.z_bit(i, qubit) {
                    continue;
                }
                let sign = if self.primary.keys.x_bit(i, qubit) {
                    -1
                } else {
                    1
                };
                let branch = self.primary.coeffs[i].clone() * sin.mul_sign(sign);
                self.primary.coeffs[i] *= cos.clone();
                let key = self.primary.keys.toggled_bits(i, qubit, true, false);
                let _ = key.key_hash();
                self.scratch.push((key, branch));
            }
            self.primary
                .reserve_for_live_len(self.primary.len() + self.scratch.len());
            if SPARSE {
                for (key, coeff) in self.scratch.drain(..) {
                    self.primary.add(key, coeff);
                }
            } else {
                for (key, coeff) in self.scratch.drain(..) {
                    self.primary.add_dense(key, coeff);
                }
            }
            return;
        }

        self.visited.clear();
        self.visited.resize(original_rows, false);
        for position in 0..scan_len {
            let i = if SPARSE {
                self.primary.sparse_rows[position] as usize
            } else {
                position
            };
            if self.visited[i]
                || self.primary.keys.is_lost(i, qubit)
                || !self.primary.keys.z_bit(i, qubit)
            {
                continue;
            }
            let key = self.primary.keys.toggled_bits(i, qubit, true, false);
            let hash = key.key_hash();
            let partner = if SPARSE {
                self.primary.find(&key, hash)
            } else {
                self.primary.find_any(&key, hash)
            };
            if let Some(j) = partner {
                debug_assert!(j < original_rows);
                self.visited[i] = true;
                self.visited[j] = true;
                let ci = self.primary.coeffs[i].clone();
                let cj = self.primary.coeffs[j].clone();
                let sign_i = if self.primary.keys.x_bit(i, qubit) {
                    -1
                } else {
                    1
                };
                let sign_j = if self.primary.keys.x_bit(j, qubit) {
                    -1
                } else {
                    1
                };
                self.primary.coeffs[i] =
                    ci.clone() * cos.clone() + cj.clone() * sin.mul_sign(sign_j);
                self.primary.coeffs[j] = cj * cos.clone() + ci * sin.mul_sign(sign_i);
            } else {
                let sign = if self.primary.keys.x_bit(i, qubit) {
                    -1
                } else {
                    1
                };
                let branch = self.primary.coeffs[i].clone() * sin.mul_sign(sign);
                self.primary.coeffs[i] *= cos.clone();
                self.visited[i] = true;
                self.scratch.push((key, branch));
            }
        }
        self.merge_rotation_scratch::<SPARSE>();
    }

    #[inline(always)]
    fn rotate_x_small_row<const SPARSE: bool>(&mut self, i: usize, qubit: usize, sin: &C, cos: &C) {
        if self.primary.keys.is_lost(i, qubit) || !self.primary.keys.z_bit(i, qubit) {
            return;
        }
        let key = self.primary.keys.toggled_bits(i, qubit, true, false);
        let hash = key.key_hash();
        let partner = if SPARSE {
            self.primary.find(&key, hash)
        } else {
            self.primary.find_any(&key, hash)
        };
        if let Some(j) = partner {
            if j < i {
                return;
            }
            let ci = self.primary.coeffs[i].clone();
            let cj = self.primary.coeffs[j].clone();
            let sign_i = if self.primary.keys.x_bit(i, qubit) {
                -1
            } else {
                1
            };
            let sign_j = if self.primary.keys.x_bit(j, qubit) {
                -1
            } else {
                1
            };
            self.primary.coeffs[i] = ci.clone() * cos.clone() + cj.clone() * sin.mul_sign(sign_j);
            self.primary.coeffs[j] = cj * cos.clone() + ci * sin.mul_sign(sign_i);
        } else {
            let sign = if self.primary.keys.x_bit(i, qubit) {
                -1
            } else {
                1
            };
            let branch = self.primary.coeffs[i].clone() * sin.mul_sign(sign);
            self.primary.coeffs[i] *= cos.clone();
            self.scratch.push((key, branch));
        }
    }

    #[inline(always)]
    fn merge_rotation_scratch<const SPARSE: bool>(&mut self) {
        self.primary
            .reserve_for_live_len(self.primary.len() + self.scratch.len());
        if SPARSE {
            for (key, coeff) in self.scratch.drain(..) {
                self.primary.add_likely_present(key, coeff);
            }
        } else {
            for (key, coeff) in self.scratch.drain(..) {
                self.primary.add_likely_present_dense(key, coeff);
            }
        }
    }
}
