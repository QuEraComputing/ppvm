// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::{config::Config, traits::U3Gate};
use ppvm_tableau::{sparsevec::SparseVector, tableau_index::TableauIndex};

use crate::{data::GeneralizedTableauSum, storage::EntryStore};

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S> U3Gate<T>
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
    fn u3(&mut self, addr0: usize, theta: T::Coeff, phi: T::Coeff, lambda: T::Coeff) {
        self.entries.for_each_mut(|tab, _p| {
            tab.u3(addr0, theta.clone(), phi.clone(), lambda.clone());
        });
        self.entries.mark_keys_dirty();
    }
}
