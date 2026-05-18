use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::{config::Config, traits::LossyMeasure};

use crate::{data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex};

pub trait LossyMeasureAll {
    fn measure_all(&mut self) -> Vec<Option<bool>>;
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> LossyMeasureAll
    for GeneralizedTableau<T, I, C>
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
{
    fn measure_all(&mut self) -> Vec<Option<bool>> {
        (0..self.n_qubits()).map(|idx| self.measure(idx)).collect()
    }
}
