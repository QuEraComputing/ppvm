use std::fmt::Debug;

use crate::data::GeneralizedTableauSum;
use crate::storage::entry_store::EntryStore;
use bitvec::view::BitView;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, PrimInt, ToPrimitive, Zero};
use ppvm_runtime::config::Config;
use ppvm_runtime::traits::Reset;
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
    fn reset(&mut self, _addr0: usize) {
        todo!("This needs to branch to account for resets where the qubit is entangled");
    }
    // impl_generalized_tableau_sum_gate!(reset, index);
}
