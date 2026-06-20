// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::PrimInt;
use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::RotationTwo;
use ppvm_tableau::sparsevec::SparseVector;
use ppvm_tableau::tableau_index::TableauIndex;

use crate::data::GeneralizedTableauSum;
use crate::storage::EntryStore;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S> RotationTwo<T>
    for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    T::Coeff: Zero + One + Send + Sync + num::Num + PartialOrd,
    I: TableauIndex + Send + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    S: EntryStore<T, I, C>,
{
    fn rotate_2(
        &mut self,
        axis_a: [u8; 2],
        axis_b: [u8; 2],
        a: usize,
        b: usize,
        theta: <T as Config>::Coeff,
    ) {
        self.entries.for_each_mut(|tab, _p| {
            tab.rotate_2(axis_a, axis_b, a, b, theta.clone());
        });
        self.entries.mark_keys_dirty();
    }
}
