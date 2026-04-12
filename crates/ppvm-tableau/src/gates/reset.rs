use std::fmt::Debug;

use crate::prelude::*;
use bitvec::view::BitView;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, PrimInt, ToPrimitive, Zero};

impl<T: Config> Reset for Tableau<T>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    fn reset(&mut self, addr0: usize) {
        let m = self.measure(addr0);
        if m {
            self.x(addr0);
        }
    }
}

impl<T, I, C> Reset for GeneralizedTableau<T, I, C>
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
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
{
    fn reset(&mut self, addr0: usize) {
        let m = self.measure(addr0);

        if let Some(true) = m {
            self.x(addr0);
        }
    }
}
