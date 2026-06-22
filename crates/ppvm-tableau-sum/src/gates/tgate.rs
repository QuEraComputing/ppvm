// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::ops::{BitAnd, Shl};

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_tableau::{sparsevec::SparseVector, tableau_index::TableauIndex};
use ppvm_traits::{config::Config, traits::TGate};

use super::impl_generalized_tableau_sum_gate;
use crate::data::GeneralizedTableauSum;
use crate::storage::EntryStore;

impl<T, I, C, S> TGate<T> for GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    impl_generalized_tableau_sum_gate!(t, index);
    impl_generalized_tableau_sum_gate!(t_adj, index);
}
