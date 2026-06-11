// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_runtime::traits::*;
use ppvm_runtime::char::Pauli;
use ppvm_runtime::config::Config;
use crate::sum::PauliSum;

impl<T: Config> Projection for PauliSum<T>
where
    T::Coeff: std::ops::MulAssign + std::ops::Neg<Output = T::Coeff> + Clone,
    T::Map: ACMapInsert<T::Storage, T::Coeff, T::BuildHasher, T::PauliWordType> + ACMapConsume,
{
    fn p0(&mut self, pos: usize) {
        self.map_insert(|k, v| {
            let half = v.half();
            match k.get(pos) {
                Pauli::I => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::Z);
                    Some((nk, v.clone()))
                }
                Pauli::Z => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::I);
                    Some((nk, v.clone()))
                }
                _ => None,
            }
        });
    }

    fn p1(&mut self, pos: usize) {
        self.map_insert(|k, v| {
            let half = v.half();
            match k.get(pos) {
                Pauli::I => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::Z);
                    Some((nk, -v.clone()))
                }
                Pauli::Z => {
                    *v *= half;
                    let nk = k.set_new(pos, Pauli::I);
                    Some((nk, -v.clone()))
                }
                _ => None,
            }
        });
    }
}
