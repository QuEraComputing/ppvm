use num::{
    Complex, One, Zero,
    complex::{Complex64, ComplexFloat},
};
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

use crate::char::Pauli;
use crate::config::Config;
use crate::tableau::GeneralizedTableau;
use crate::tableau::sparsevec::SparseVector;
use crate::traits::*;

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> RotationOne<T>
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: Zero + One,
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize>
        + BitOrAssign<<I as Shl<usize>>::Output>
        + BitAnd<<I as Shl<usize>>::Output, Output = I>
        + BitXor<Output = I>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
{
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: <T as Config>::Coeff) {
        let (cos, sin) = theta.sin_cos();

        let complex_cos: Complex<T::Coeff> = Complex {
            re: cos,
            im: T::Coeff::zero(),
        };

        let i_complex_sin: Complex<T::Coeff> = Complex {
            re: T::Coeff::zero(),
            im: -sin,
        };

        self.branch_with_coefficients(addr0, axis, complex_cos, i_complex_sin);
    }
}
