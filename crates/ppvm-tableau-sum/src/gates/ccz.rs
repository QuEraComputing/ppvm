// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_tableau::{sparsevec::SparseVector, tableau_index::TableauIndex};
use ppvm_traits::{config::Config, traits::CCZGate};

use crate::{data::GeneralizedTableauSum, storage::EntryStore};

impl<T, I, C, S> CCZGate<T> for GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    <T::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
{
    fn ccz(&mut self, a: usize, b: usize, c: usize) {
        assert!(a != b && a != c && b != c, "CCZ requires distinct qubits");
        assert!(a < self.n_qubits && b < self.n_qubits && c < self.n_qubits);
        self.entries.for_each_mut(|tab, _| tab.ccz(a, b, c));
        self.entries.mark_keys_dirty();
    }
}
