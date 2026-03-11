use num::Complex;
use num::complex::{Complex64, ComplexFloat};
use num::{One, ToPrimitive, Zero};

use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use crate::tableau::{GeneralizedTableau, Tableau};
use crate::traits::Clifford;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

pub trait TGate<T: Config> {
    fn t(&mut self, addr0: usize);
    fn t_adj(&mut self, addr0: usize);
}

pub trait Measure {
    fn measure(&mut self, addr0: usize) -> bool;
}

pub trait CliffordExtensions: Clifford {
    fn s_adj(&mut self, addr0: usize);
    fn sqrt_x(&mut self, addr0: usize) {
        self.h(addr0);
        self.s(addr0);
        self.h(addr0);
    }

    fn sqrt_x_adj(&mut self, addr0: usize) {
        self.h(addr0);
        self.s_adj(addr0);
        self.h(addr0);
    }

    fn sqrt_y(&mut self, addr0: usize) {
        // NOTE: in matmul (RHS applied first)
        // SqrtY == S * SqrtX * S'
        self.s_adj(addr0);
        self.sqrt_x(addr0);
        self.s(addr0);
    }

    fn sqrt_y_adj(&mut self, addr0: usize) {
        self.s_adj(addr0);
        self.sqrt_x_adj(addr0);
        self.s(addr0);
    }
}

pub trait TableauIndex:
    PartialEq
    + Eq
    + Hash
    + Copy
    + From<u8>
    + Shl<usize, Output = Self>
    + BitOrAssign<Self>
    + BitAnd<Self, Output = Self>
    + BitXor<Output = Self>
{
}

impl<I> TableauIndex for I where
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize, Output = Self>
        + BitOrAssign<Self>
        + BitAnd<I, Output = I>
        + BitXor<Output = I>
{
}

pub trait Reset: Clifford + Measure {
    fn reset(&mut self, addr0: usize) {
        let m = self.measure(addr0);
        if m {
            self.x(addr0);
        }
    }
}

impl<T: Config> Reset for Tableau<T> {}

impl<T, I, C> Reset for GeneralizedTableau<T, I, C>
where
    T: Config,
    I: TableauIndex + Debug,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
{
}
