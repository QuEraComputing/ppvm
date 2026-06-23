// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use super::data::LossyPauliWord;
use ppvm_traits::traits::{HashFinalize, PauliStorage, PauliWordTrait, Trace};

impl<'a, A, H> Trace<'a, LossyPauliWord<A, H>> for LossyPauliWord<A, H>
where
    A: PauliStorage + 'a,
    H: Default + BuildHasher + Clone + HashFinalize + 'a,
{
    type Output = bool;
    fn trace(&'a self, value: &'a LossyPauliWord<A, H>) -> Self::Output {
        debug_assert_eq!(
            self.n_qubits(),
            value.n_qubits(),
            "#qubits mismatch, got {} and {}",
            self.n_qubits(),
            value.n_qubits()
        );
        self == value
    }
}
