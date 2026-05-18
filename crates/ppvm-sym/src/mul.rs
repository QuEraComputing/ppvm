// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

use core::panic;

use ppvm_runtime::prelude::*;

use crate::{
    Term,
    term::{Inner, Prod, Sum},
};

impl_op_mul_assign_coefficient!(Term);

impl Prod {
    /// Multiply this product by an additional `sin(x_u)`.
    pub fn mul_sin(&mut self, u: u32) {
        *self.sin.entry(u).or_insert(0) += 1;
        self.sin_pow += 1;
    }

    /// Multiply this product by an additional `cos(x_u)`.
    pub fn mul_cos(&mut self, u: u32) {
        *self.cos.entry(u).or_insert(0) += 1;
        self.cos_pow += 1;
    }
}

impl Sum {
    /// Multiply this sum in place by `coeff · p`, respecting the same
    /// `max` / `min_eps` truncation bounds used elsewhere.
    pub fn mul_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max || coeff.abs() < min_eps {
            self.terms.clear();
            self.c0 = 0.0;
            return;
        }

        if p.pow() == 0 {
            *self *= coeff;
            return;
        }

        let mut old_terms = std::mem::take(&mut self.terms);
        self.add_term(p.clone(), self.c0 * coeff, max, min_eps);
        for (k, v) in old_terms.drain() {
            self.add_term(k * p.clone(), v * coeff, max, min_eps);
        }
        self.c0 = 0.0;
    }
}

impl std::ops::Mul<Prod> for Prod {
    type Output = Prod;

    fn mul(self, rhs: Prod) -> Self::Output {
        let mut new = self;
        new *= rhs;
        new
    }
}

impl std::ops::MulAssign<Prod> for Prod {
    fn mul_assign(&mut self, rhs: Prod) {
        for (k, v) in rhs.sin {
            *self.sin.entry(k).or_insert(0) += v;
        }
        self.sin_pow += rhs.sin_pow;
        for (k, v) in rhs.cos {
            *self.cos.entry(k).or_insert(0) += v;
        }
        self.cos_pow += rhs.cos_pow;
    }
}

impl std::ops::MulAssign<f64> for Sum {
    fn mul_assign(&mut self, rhs: f64) {
        self.c0 *= rhs;
        for v in self.terms.values_mut() {
            *v *= rhs;
        }
    }
}

impl std::ops::MulAssign<f64> for Term {
    fn mul_assign(&mut self, rhs: f64) {
        match self.inner {
            Inner::Sum(ref mut s) => {
                *s *= rhs;
            }
            Inner::One(_, ref mut c) => {
                *c *= rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
            Inner::Const(ref mut c) => {
                *c *= rhs;
            }
        }
    }
}

impl std::ops::MulAssign<Term> for Term {
    fn mul_assign(&mut self, rhs: Term) {
        match self.inner {
            Inner::Sum(ref mut s1) => match rhs.inner {
                Inner::Sum(ref s2) => {
                    let mut new_sum = Sum::new();
                    new_sum.c0 = s1.c0 * s2.c0;
                    if s2.c0 > self.min_eps {
                        for (p1, c1) in s1.terms.iter() {
                            let p = p1.clone();
                            new_sum.add_term(p, c1 * s2.c0, self.max_sin, self.min_eps);
                        }
                    }

                    if s1.c0 > self.min_eps {
                        for (p2, c2) in s2.terms.iter() {
                            let p = p2.clone();
                            new_sum.add_term(p, c2 * s1.c0, self.max_sin, self.min_eps);
                        }
                    }

                    for (p1, c1) in s1.terms.iter() {
                        for (p2, c2) in s2.terms.iter() {
                            new_sum.add_term(
                                p1.clone() * p2.clone(),
                                c1 * c2,
                                self.max_sin,
                                self.min_eps,
                            );
                        }
                    }
                    self.inner = Inner::Sum(new_sum);
                }
                Inner::Const(c) => {
                    *s1 *= c;
                }
                Inner::One(p2, c2) => {
                    s1.mul_term(p2, c2, self.max_sin, self.min_eps);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Const(c) => match rhs.inner {
                Inner::Sum(s) => {
                    self.inner = Inner::Sum(s);
                    *self *= c;
                }
                Inner::Const(c2) => {
                    self.inner = Inner::Const(c * c2);
                }
                Inner::One(p2, c2) => {
                    self.inner = Inner::One(p2, c * c2);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::One(ref p, c) => match rhs.inner {
                Inner::Sum(s) => {
                    let mut new_sum = s;
                    new_sum.mul_term(p.clone(), c, self.max_sin, self.min_eps);
                }
                Inner::Const(c2) => {
                    self.inner = Inner::One(p.clone(), c * c2);
                }
                Inner::One(p2, c2) => {
                    self.inner = Inner::One(p.clone() * p2, c * c2);
                }
                Inner::Var(_) => {
                    panic!(
                        "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions, bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::Mul<f64> for Term {
    type Output = Term;

    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl std::ops::Mul<Term> for f64 {
    type Output = Term;

    fn mul(self, rhs: Term) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl std::ops::Mul<Term> for Term {
    type Output = Term;

    fn mul(self, rhs: Term) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}
