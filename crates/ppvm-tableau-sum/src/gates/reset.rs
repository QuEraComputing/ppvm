// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;

use crate::data::GeneralizedTableauSum;
use crate::storage::entry_store::EntryStore;
use bitvec::view::BitView;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, PrimInt, ToPrimitive, Zero};
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::{Clifford, Reset};
use ppvm_tableau::sparsevec::SparseVector;
use ppvm_tableau::tableau_index::TableauIndex;

impl<T, I, C, S> Reset for GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Debug + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + PartialOrd
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    S: EntryStore<T, I, C>,
{
    fn reset(&mut self, addr0: usize) {
        self.for_each_z_branch(addr0, |tab, outcome, _p| {
            // Flip back to |0⟩ on the outcome-1 branch; the helper takes
            // care of fingerprint hygiene and the branch / merge plumbing.
            if outcome == Some(true) {
                tab.x(addr0);
                true
            } else {
                false
            }
        });
    }
}
