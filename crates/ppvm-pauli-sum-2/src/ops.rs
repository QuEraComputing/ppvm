// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Arithmetic accumulates coefficients without reducing or truncating terms.
//! Results retain the left operand's policy and preserved keys.

use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use ppvm_traits_2::{Accumulate, Scale, Word};

use crate::Sum;

impl<S: Accumulate, P, K: Into<S::Key>> AddAssign<(K, S::Coeff)> for Sum<S, P>
where
    S::Key: Word,
{
    fn add_assign(&mut self, (key, coeff): (K, S::Coeff)) {
        self.add_term(key.into(), coeff);
    }
}

impl<S: Accumulate, P> AddAssign<&Self> for Sum<S, P> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.n_sites(),
            rhs.n_sites(),
            "addition requires equal widths"
        );
        for (key, coeff) in rhs.iter() {
            self.storage.accumulate(key, coeff);
        }
    }
}

impl<S: Accumulate, P> AddAssign<Self> for Sum<S, P> {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl<S: Accumulate, P> Add<&Self> for Sum<S, P> {
    type Output = Self;

    fn add(mut self, rhs: &Self) -> Self {
        self += rhs;
        self
    }
}

impl<S: Accumulate, P> Add<Self> for Sum<S, P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        self + &rhs
    }
}

impl<S: Accumulate + Clone, P: Clone> Add<&Sum<S, P>> for &Sum<S, P> {
    type Output = Sum<S, P>;

    fn add(self, rhs: &Sum<S, P>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<S: Accumulate + Clone, P: Clone> Add<Sum<S, P>> for &Sum<S, P> {
    type Output = Sum<S, P>;

    fn add(self, rhs: Sum<S, P>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<S: Accumulate, P> SubAssign<&Self> for Sum<S, P> {
    fn sub_assign(&mut self, rhs: &Self) {
        assert_eq!(
            self.n_sites(),
            rhs.n_sites(),
            "subtraction requires equal widths"
        );
        for (key, coeff) in rhs.iter() {
            self.storage.accumulate(key, -coeff);
        }
    }
}

impl<S: Accumulate, P> SubAssign<Self> for Sum<S, P> {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl<S: Accumulate, P> Sub<&Self> for Sum<S, P> {
    type Output = Self;

    fn sub(mut self, rhs: &Self) -> Self {
        self -= rhs;
        self
    }
}

impl<S: Accumulate, P> Sub<Self> for Sum<S, P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self - &rhs
    }
}

impl<S: Accumulate + Clone, P: Clone> Sub<&Sum<S, P>> for &Sum<S, P> {
    type Output = Sum<S, P>;

    fn sub(self, rhs: &Sum<S, P>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<S: Accumulate + Clone, P: Clone> Sub<Sum<S, P>> for &Sum<S, P> {
    type Output = Sum<S, P>;

    fn sub(self, rhs: Sum<S, P>) -> Self::Output {
        self.clone() - rhs
    }
}

// Each scalar type supports assignment, an owned sum, and a borrowed sum.
macro_rules! scalar_ops {
    ($coeff:ty, $scalar:ty) => {
        impl<S: Scale<Coeff = $coeff>, P> MulAssign<$scalar> for Sum<S, P> {
            fn mul_assign(&mut self, rhs: $scalar) {
                self.scale(&rhs);
            }
        }

        impl<S: Scale<Coeff = $coeff>, P> Mul<$scalar> for Sum<S, P> {
            type Output = Self;

            fn mul(mut self, rhs: $scalar) -> Self {
                self *= rhs;
                self
            }
        }

        impl<S: Scale<Coeff = $coeff> + Clone, P: Clone> Mul<$scalar> for &Sum<S, P> {
            type Output = Sum<S, P>;

            fn mul(self, rhs: $scalar) -> Self::Output {
                self.clone() * rhs
            }
        }
    };
}

scalar_ops!(f64, f64);
scalar_ops!(f64, &f64);
scalar_ops!(num::Complex<f64>, num::Complex<f64>);
scalar_ops!(num::Complex<f64>, &num::Complex<f64>);
