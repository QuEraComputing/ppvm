// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use num::Complex;
use ppvm_runtime::{
    config::Config,
    traits::{Clifford, CliffordExtensions},
};
use ppvm_tableau::sparsevec::SparseVector;

use super::impl_generalized_tableau_sum_gate;
use crate::data::GeneralizedTableauSum;
use crate::storage::EntryStore;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> Clifford
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_sum_gate!(x, index);
    impl_generalized_tableau_sum_gate!(y, index);
    impl_generalized_tableau_sum_gate!(z, index);
    impl_generalized_tableau_sum_gate!(h, index);
    impl_generalized_tableau_sum_gate!(s, index);
    impl_generalized_tableau_sum_gate!(cnot, control, target);
    impl_generalized_tableau_sum_gate!(cz, control, target);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>, S: EntryStore<T, I, C>> CliffordExtensions
    for GeneralizedTableauSum<T, I, C, S>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_sum_gate!(s_adj, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_x, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_x_adj, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_y, addr0);
    impl_generalized_tableau_sum_gate!(sqrt_y_adj, addr0);
    impl_generalized_tableau_sum_gate!(cy, addr0, addr1);
}
