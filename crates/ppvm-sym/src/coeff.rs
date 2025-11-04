use ppvm_runtime::traits::{Coefficient, ComplexCoefficient};

use crate::{Prod, Sum, Term, term::Inner};

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

    fn cutoff(&self, threshold: f64) -> bool {
        if let Inner::Const(c) = self.inner {
            c.abs() < threshold
        } else {
            false
        }
    }
}

impl ComplexCoefficient for Term {
    fn mul_phase(&self, phase: u8) -> Self {
        match self.inner {
            Inner::Sum(ref s) => {
                let mut ret = Sum::new();
                for (p, c) in s.terms.iter() {
                    let mut new_p = p.clone();
                    new_p.add_phase(phase);
                    ret.add_term(new_p, *c, self.max_sin, self.min_eps);
                }
                let mut c0 = Prod::new();
                c0.add_phase(phase);
                ret.add_term(c0, s.c0, self.max_sin, self.min_eps);
                Term {
                    inner: Inner::Sum(ret),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Const(f) => {
                let mut ret = Prod::new();
                ret.add_phase(phase);
                Term {
                    inner: Inner::One(ret, f),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::One(ref p, c) => {
                let mut ret = p.clone();
                ret.add_phase(phase);
                Term {
                    inner: Inner::One(ret, c),
                    max_sin: self.max_sin,
                    min_eps: self.min_eps,
                }
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                );
            }
        }
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
