// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Addition: the additive group on [`Prod`]-keyed [`Sum`]s and on [`Term`].
//!
//! Ported from `ppvm-sym/src/add.rs`, arm for arm, including its two
//! deliberately *untruncated* operator forms (`Sum += f64` and `Sum += Prod`,
//! behavioural contract 11) and the `Const + Sum` **operand adoption** that moves
//! the right-hand table into the accumulator (perf feature 6). The one
//! divergence is the `One` receiver of `Term += f64`, documented on that impl.

use crate::term::{Inner, Prod, Sum, Term};

/// Bare `c₀ += rhs`, with **no `min_eps` check** — asymmetric with the method
/// form [`Sum::add_const`], which does drop `|c| < min_eps`. Old behaviour,
/// reproduced deliberately (behavioural contract 11).
impl std::ops::AddAssign<f64> for Sum {
    fn add_assign(&mut self, rhs: f64) {
        self.c0 += rhs;
    }
}

/// Bare table insert, with **no `max_sin`/`min_eps` check at all** — so
/// `sum += prod` can insert a monomial [`Sum::add_term`] would have rejected,
/// including a `pow() == 0` monomial. Old behaviour, reproduced deliberately
/// (behavioural contract 11).
impl std::ops::AddAssign<Prod> for Sum {
    fn add_assign(&mut self, rhs: Prod) {
        *self.terms.entry(rhs).or_insert(0.0) += 1.0;
    }
}

/// # Divergence from old (`oldSuspectedBugs` #1)
///
/// Old's `Inner::One` arm (`ppvm-sym/src/add.rs:24-28`) built the promoted `Sum`
/// into a local and **never assigned it back**, so `t += 2.0` was a silent no-op
/// whenever `t` was a single weighted monomial — the most common non-constant
/// form during propagation. It is a copy-paste omission: the structurally
/// identical `AddAssign<Prod> for Term` arm two impls below *does* assign, and
/// old's two `add` unit tests only exercise `Const` and `Sum` receivers.
///
/// The coefficient ring must be an additive monoid for the graded-map laws
/// (`lean/PPVM/Algebra/GradedMap.lean` `accumulate_comm`/`accumulate_assoc`),
/// which forbids `x + c == x` for `c != 0`, so this impl is **correct** here and
/// the divergence is pinned by a differential test rather than ported. Nothing in
/// the propagation path reaches it (the engine only ever adds `Term`s, and
/// `AddAssign<Term>`'s own `One` receiver arms already assign), so the golden
/// masters are unaffected.
impl std::ops::AddAssign<f64> for Term {
    fn add_assign(&mut self, rhs: f64) {
        match self.inner {
            Inner::Const(ref mut c) => {
                *c += rhs;
            }
            Inner::One(ref p, c) => {
                let mut sum = Sum::new();
                sum.c0 = rhs;
                sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                self.inner = Inner::Sum(sum);
            }
            Inner::Sum(ref mut s) => {
                *s += rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::AddAssign<Prod> for Term {
    fn add_assign(&mut self, rhs: Prod) {
        match self.inner {
            Inner::Const(c) => {
                let mut sum = Sum::new();
                sum.c0 = c;
                sum.add_term(rhs, 1.0, self.max_sin, self.min_eps);
                self.inner = Inner::Sum(sum);
            }
            Inner::One(ref p, c) => {
                let mut sum = Sum::new();
                sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                sum.add_term(rhs, 1.0, self.max_sin, self.min_eps);
                self.inner = Inner::Sum(sum);
            }
            Inner::Sum(ref mut s) => {
                *s += rhs;
            }
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

/// The truncation parameters are inherited from **`self` only**; `rhs`'s
/// `max_sin`/`min_eps` are silently ignored (behavioural contract 1).
impl std::ops::AddAssign<Term> for Term {
    fn add_assign(&mut self, rhs: Term) {
        match self.inner {
            Inner::Sum(ref mut s1) => match rhs.inner {
                Inner::Sum(s2) => {
                    for (p, c) in s2.terms {
                        s1.add_term(p, c, self.max_sin, self.min_eps);
                    }
                    s1.c0 += s2.c0;
                }
                Inner::One(p, c) => {
                    s1.add_term(p, c, self.max_sin, self.min_eps);
                }
                Inner::Const(c) => {
                    s1.add_const(c, self.min_eps);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::One(ref p, c) => match rhs.inner {
                // Operand adoption: the right-hand table is MOVED in.
                Inner::Sum(s) => {
                    let mut new_sum = s;
                    new_sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(new_sum);
                }
                Inner::Const(c2) => {
                    let mut sum = Sum::new();
                    sum.c0 = c2;
                    sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                Inner::One(p2, c2) => {
                    let mut sum = Sum::new();
                    sum.add_term(p.clone(), c, self.max_sin, self.min_eps);
                    sum.add_term(p2, c2, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Const(c) => match rhs.inner {
                // Operand adoption: `std::iter::Sum`'s fold starts here, so the
                // first summand's table becomes the accumulator (behavioural
                // contract 10).
                Inner::Sum(s) => {
                    self.inner = Inner::Sum(s);
                    *self += c;
                }
                Inner::One(p2, c2) => {
                    let mut sum = Sum::new();
                    sum.c0 = c;
                    sum.add_term(p2, c2, self.max_sin, self.min_eps);
                    self.inner = Inner::Sum(sum);
                }
                Inner::Const(c2) => {
                    self.inner = Inner::Const(c + c2);
                }
                _ => {
                    panic!(
                        "variable is not used in sin/cos expressions,\
                            bare variable is not allowed in expression"
                    );
                }
            },
            Inner::Var(_) => {
                panic!(
                    "variable is not used in sin/cos expressions,\
                    bare variable is not allowed in expression"
                );
            }
        }
    }
}

impl std::ops::Add for Term {
    type Output = Term;

    fn add(mut self, rhs: Term) -> Self::Output {
        self += rhs;
        self
    }
}

impl std::ops::Add<f64> for Term {
    type Output = Term;

    fn add(mut self, rhs: f64) -> Self::Output {
        self += rhs;
        self
    }
}

impl std::ops::Add<Term> for f64 {
    type Output = Term;

    fn add(self, mut rhs: Term) -> Self::Output {
        rhs += self;
        rhs
    }
}

impl std::ops::Sub<Term> for Term {
    type Output = Term;

    fn sub(mut self, rhs: Term) -> Self::Output {
        self += -rhs;
        self
    }
}

impl std::ops::Sub<f64> for Term {
    type Output = Term;

    fn sub(mut self, rhs: f64) -> Self::Output {
        self += -rhs;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_add_sum() {
        let mut sum = Sum::new();
        sum += 1.0;
        sum += Prod::sin(0);
        assert_eq!(sum.eval(&[1.1]).unwrap(), 1.0 + 1.1f64.sin());
    }

    #[test]
    fn test_add_term() {
        let mut t = Term::from_f64(2.0);
        t += Prod::sin(0);
        t += Prod::cos(1);
        assert_eq!(
            t.eval(&[1.1, 2.2]).unwrap(),
            2.0 + 1.1f64.sin() + 2.2f64.cos()
        );

        let mut t = Term::var(0).sin();
        t += Term::from_f64(3.0);
        t += Term::var(1).cos();
        assert_eq!(
            t.eval(&[1.1, 2.2]).unwrap(),
            3.0 + 1.1f64.sin() + 2.2f64.cos()
        );
    }

    #[test]
    fn operator_add_bypasses_truncation() {
        // Behavioural contract 11.
        let mut sum = Sum::new();
        sum += 1e-9;
        assert_eq!(sum.c0(), 1e-9);

        let mut sum = Sum::new();
        sum += Prod::new();
        assert_eq!(sum.len(), 1, "`+= Prod` inserts even an empty monomial");
    }

    #[test]
    fn add_const_drops_below_min_eps() {
        let mut sum = Sum::new();
        sum.add_const(1e-9, 1e-3);
        assert_eq!(sum.c0(), 0.0);
    }

    #[test]
    fn add_f64_to_one_promotes() {
        // Divergence from old (`oldSuspectedBugs` #1): old silently dropped the
        // addend here.
        let mut t = Term::var(0).sin();
        t += 2.0;
        assert_eq!(t.eval(&[0.5]).unwrap(), 2.0 + 0.5f64.sin());
    }

    #[test]
    fn exact_cancellation_keeps_a_zero_monomial() {
        // Behavioural contract 4.
        let mut t = Term::var(0).sin();
        t += Term::var(0).sin() * -1.0;
        assert_eq!(t.n_monomials(), 1);
        assert!(t.to_string().contains("0.000"));
    }

    #[test]
    fn truncation_parameters_come_from_the_left() {
        // Behavioural contract 1: LHS wins, so swapping operands changes the
        // result. (A map-backed receiver is needed to observe it: the
        // `Const × One` / `One × One` fast arms never consult `max_sin` at all,
        // in old exactly as here.)
        let mut a = Term::from_f64(1.0) + Term::var(1).cos();
        a.set_max_sin(0);
        let b = Term::var(0).sin(); // max_sin = usize::MAX
        let ab = a.clone() * b.clone();
        let ba = b * a;
        assert_eq!(ab.n_monomials(), 0, "LHS max_sin = 0 truncated everything");
        assert!(ba.n_monomials() > 0, "LHS max_sin = MAX kept the monomials");
    }
}
