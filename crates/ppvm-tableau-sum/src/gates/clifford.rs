// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::Complex;
use num::PrimInt;
use ppvm_tableau::sparsevec::SparseVector;
use ppvm_traits::config::Config;
use ppvm_traits::traits::{Clifford, CliffordExtensions};

use super::impl_generalized_tableau_sum_gate;
use crate::data::GeneralizedTableauSum;
use crate::storage::EntryStore;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> Clifford
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    impl_generalized_tableau_sum_gate!(x);
    impl_generalized_tableau_sum_gate!(y);
    impl_generalized_tableau_sum_gate!(z);
    impl_generalized_tableau_sum_gate!(h);
    impl_generalized_tableau_sum_gate!(s);
    impl_generalized_tableau_sum_gate!(cnot);
    impl_generalized_tableau_sum_gate!(cz);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> CliffordExtensions
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    impl_generalized_tableau_sum_gate!(s_dag);
    impl_generalized_tableau_sum_gate!(sqrt_x);
    impl_generalized_tableau_sum_gate!(sqrt_x_dag);
    impl_generalized_tableau_sum_gate!(sqrt_y);
    impl_generalized_tableau_sum_gate!(sqrt_y_dag);
    impl_generalized_tableau_sum_gate!(cy);
}
