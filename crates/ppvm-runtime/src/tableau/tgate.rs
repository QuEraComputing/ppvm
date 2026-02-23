use super::sparsevec::SparseVector;
use super::traits::TGate;
use crate::config::Config;
use crate::tableau::GeneralizedTableau;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, Zero};
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

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
    C: SparseVector<Complex<T::Coeff>, I>,
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize>
        + BitOrAssign<<I as Shl<usize>>::Output>
        + BitAnd<<I as Shl<usize>>::Output, Output = I>
        + BitXor<Output = I>,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn t(&mut self, index: usize) {
        if self.is_lost[index] {
            return;
        }

        let complex_cos: Complex<T::Coeff> = COS_PI_OVER_8_TIMES_EXPIPI8.into();
        let complex_sin: Complex<T::Coeff> = ISIN_PI_OVER_8_TIMES_EXPIPI8.into();
        self.branch_with_coefficients(index, crate::char::Pauli::Z, complex_cos, complex_sin);
    }

    fn t_adj(&mut self, index: usize) {
        if self.is_lost[index] {
            return;
        }

        let complex_cos: Complex<T::Coeff> = COS_PI_OVER_8_TIMES_EXPIPI8.conj().into();
        let complex_sin: Complex<T::Coeff> = ISIN_PI_OVER_8_TIMES_EXPIPI8.conj().into();
        self.branch_with_coefficients(index, crate::char::Pauli::Z, complex_cos, complex_sin);
    }
}
