use ppvm_runtime::traits::Coefficient;

use crate::{Term, term::Inner};

impl Coefficient for Term {
    fn half(&self) -> Self {
        Term::from_f64(0.5)
    }

    fn mul_sign(&self, sign: i8) -> Self {
        let mut ret = self.clone();
        ret *= sign as f64;
        ret
    }

    fn sin_cos(&self) -> (Self, Self) {
        (self.clone().sin(), self.clone().cos())
    }
}

impl std::iter::Sum for Term {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut total = Term::from_f64(0.0);
        for t in iter {
            total += t;
        }
        total
    }
}

impl std::ops::Neg for Term {
    type Output = Term;

    fn neg(self) -> Self::Output {
        let mut ret = self;
        ret *= -1.0;
        ret
    }
}

impl From<f32> for Term {
    fn from(value: f32) -> Self {
        Term::from_f64(value as f64)
    }
}

impl From<f64> for Term {
    fn from(value: f64) -> Self {
        Term::from_f64(value)
    }
}

impl From<i32> for Term {
    fn from(value: i32) -> Self {
        Term::from_f64(value as f64)
    }
}

impl From<i64> for Term {
    fn from(value: i64) -> Self {
        Term::from_f64(value as f64)
    }
}

impl num::Zero for Term {
    fn zero() -> Self {
        Term::from_f64(0.0)
    }

    fn is_zero(&self) -> bool {
        if let Inner::Const(c) = self.inner {
            c.abs() < self.min_eps
        } else {
            false
        }
    }
}
