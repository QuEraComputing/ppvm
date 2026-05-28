use ppvm_runtime::config::Config;
use ppvm_runtime::traits::LossyMeasure;
use ppvm_tableau::sparsevec::SparseVector;
use ppvm_tableau::tableau_index::TableauIndex;

use crate::prelude::*;
use crate::storage::EntryStore;
use bitvec::view::BitView;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use std::fmt::Debug;

impl<T, I, C, S> LossyMeasure for GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
    S: EntryStore<T, I, C>,
{
    fn measure(&mut self, _addr0: usize) -> Option<bool> {
        todo!(
            "Measure needs to branch and needs a different return value (e.g. probabilities for 0, 1, lost)"
        )
    }
}
