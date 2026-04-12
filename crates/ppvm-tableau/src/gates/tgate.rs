use crate::prelude::*;
use bitvec::view::BitView;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, Zero};
use std::ops::{BitAnd, Shl};

const COS_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.8535533905932737,
    im: 0.3535533905932738,
}; // exp(im * pi / 8) * cos(pi/8)
const ISIN_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.14644660940672624,
    im: -0.3535533905932738,
}; // -im * exp(im * pi / 8) * sin(pi/8)

impl<T, I, C> TGate<T> for GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I>,
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn t(&mut self, index: usize) {
        if self.is_lost[index] {
            return;
        }

        let complex_cos: Complex<T::Coeff> = COS_PI_OVER_8_TIMES_EXPIPI8.into();
        let complex_sin: Complex<T::Coeff> = ISIN_PI_OVER_8_TIMES_EXPIPI8.into();
        self.branch_with_coefficients(index, Pauli::Z, complex_cos, complex_sin);
    }

    fn t_adj(&mut self, index: usize) {
        if self.is_lost[index] {
            return;
        }

        let complex_cos: Complex<T::Coeff> = COS_PI_OVER_8_TIMES_EXPIPI8.conj().into();
        let complex_sin: Complex<T::Coeff> = ISIN_PI_OVER_8_TIMES_EXPIPI8.conj().into();
        self.branch_with_coefficients(index, Pauli::Z, complex_cos, complex_sin);
    }
}
