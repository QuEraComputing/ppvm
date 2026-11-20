// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Addition and subtraction on [`Term`], truncating with the left operand's
//! settings. A [`Sum`] carries no settings, so it is only changed through
//! [`Sum::add_term`]/[`Sum::mul_term`], which take them explicitly.

use crate::term::{Inner, Sum, Term};

impl Sum {
    /// Add a non-`Sum` representation, with the truncation of [`Sum::add_term`].
    #[inline]
    fn add_small(&mut self, rhs: Inner, max: usize, min_eps: f64) {
        match rhs {
            Inner::One(p, c) => self.add_term(p, c, max, min_eps),
            Inner::Const(c) => self.c0 += c,
            Inner::Sum(_) => unreachable!(),
        }
    }
}

/// A `Sum` operand's table is reused (moved in), whichever side it is on; a
/// right-hand table is pruned to `self`'s settings when its own are looser.
impl std::ops::AddAssign<Term> for Term {
    #[inline]
    fn add_assign(&mut self, rhs: Term) {
        let (max, eps) = (self.max_sin, self.min_eps);
        let looser = self.is_looser(&rhs);
        self.inner = match (std::mem::take(&mut self.inner), rhs.inner) {
            (Inner::Const(a), Inner::Const(b)) => Inner::Const(a + b),
            (Inner::Sum(mut s1), Inner::Sum(s2)) => {
                for (p, c) in s2.maps.into_iter().flat_map(|m| m.terms) {
                    s1.add_term(p, c, max, eps);
                }
                s1.c0 += s2.c0;
                Inner::Sum(s1)
            }
            (Inner::Sum(mut s), small) => {
                s.add_small(small, max, eps);
                Inner::Sum(s)
            }
            (small, Inner::Sum(mut s)) => {
                if looser {
                    s.prune(max, eps);
                }
                s.add_small(small, max, eps);
                Inner::Sum(s)
            }
            (a, b) => {
                let mut s = Sum::new();
                s.add_small(a, max, eps);
                s.add_small(b, max, eps);
                Inner::Sum(s)
            }
        };
    }
}

impl Term {
    /// `self += sign · rhs` without cloning a `Sum` right operand's table.
    fn add_scaled_ref(&mut self, rhs: &Term, sign: f64) {
        if let (Inner::Sum(lhs), Inner::Sum(rhs)) = (&mut self.inner, &rhs.inner) {
            for (p, c) in rhs.iter() {
                lhs.add_term(p.clone(), sign * c, self.max_sin, self.min_eps);
            }
            lhs.c0 += sign * rhs.c0;
            return;
        }
        let mut rhs = rhs.clone();
        if sign < 0.0 {
            rhs.negate();
        }
        *self += rhs;
    }
}

impl std::ops::AddAssign<&Term> for Term {
    #[inline]
    fn add_assign(&mut self, rhs: &Term) {
        self.add_scaled_ref(rhs, 1.0);
    }
}

impl std::ops::SubAssign<&Term> for Term {
    #[inline]
    fn sub_assign(&mut self, rhs: &Term) {
        self.add_scaled_ref(rhs, -1.0);
    }
}

impl std::ops::SubAssign<Term> for Term {
    #[inline]
    #[allow(clippy::suspicious_op_assign_impl)] // subtraction adds the negation
    fn sub_assign(&mut self, mut rhs: Term) {
        rhs.negate();
        *self += rhs;
    }
}

impl std::ops::AddAssign<f64> for Term {
    #[inline]
    fn add_assign(&mut self, rhs: f64) {
        *self += Term::from_f64(rhs);
    }
}

impl std::ops::SubAssign<f64> for Term {
    #[inline]
    fn sub_assign(&mut self, rhs: f64) {
        *self += Term::from_f64(-rhs);
    }
}

binary_ops!(Add::add from AddAssign::add_assign for Term, &Term, f64);
binary_ops!(Sub::sub from SubAssign::sub_assign for Term, &Term, f64);

impl std::ops::Add<Term> for f64 {
    type Output = Term;

    #[inline]
    fn add(self, rhs: Term) -> Term {
        rhs + self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::{Prod, Var};
    use ppvm_traits_2::ImaginaryUnit;

    #[test]
    fn term_operators_truncate_with_the_left_settings() {
        let mut sum = Sum::new();
        // The identity folds into c₀, which no threshold truncates.
        sum.add_term(Prod::new(), 1.0, 0, 1e-3);
        sum.add_term(Prod::new(), 1e-9, 0, 1e-3);
        assert_eq!((sum.len(), sum.c0()), (0, 1.0 + 1e-9));

        let mut t = Term::from(1.0) + Var(1).cos();
        t.set_max_sin(0);
        t.set_min_eps(1e-3);
        t += Var(0).sin(); // above `max_sin`: dropped
        t += Var(1).cos() * 1e-6; // below `min_eps`: dropped
        t += 1e-9; // a constant: kept
        assert_eq!(t.inner(), (Term::from(1.0 + 1e-9) + Var(1).cos()).inner());
    }

    #[test]
    fn cancellation_removes_the_monomial() {
        let mut t = Var(0).sin();
        t += Var(0).sin() * -1.0;
        assert_eq!(t.n_monomials(), 0);
        assert_eq!(t.to_string(), "[]");
        // A merge that lands below `min_eps` (not exactly 0) is removed too.
        let mut t = Var(0).sin() + Var(1).cos();
        t.set_min_eps(1e-3);
        t += Var(0).sin() * -0.9995;
        assert_eq!(t.n_monomials(), 1);
    }

    #[test]
    fn scaling_down_prunes_but_a_sign_flip_does_not_walk_twice() {
        let mut t = Var(0).sin() + Var(1).cos() * 1e-3;
        t.set_min_eps(1e-6);
        t *= -1.0;
        assert_eq!(t.n_monomials(), 2);
        t *= 1e-4;
        assert_eq!(t.n_monomials(), 1, "1e-7 fell below min_eps");
        t *= 0.0;
        assert_eq!(t.n_monomials(), 0);
    }

    #[test]
    fn truncation_parameters_come_from_the_left() {
        // The left operand's `max_sin` wins.
        let mut a = Term::from_f64(1.0) + Var(1).cos();
        a.set_max_sin(0);
        let b = Var(0).sin(); // max_sin = usize::MAX
        let ab = a.clone() * b.clone();
        let ba = b * a;
        assert_eq!(ab.n_monomials(), 0, "LHS max_sin = 0 truncated everything");
        assert!(ba.n_monomials() > 0, "LHS max_sin = MAX kept the monomials");
    }

    #[test]
    fn identity_monomial_folds_into_c0_through_term_ops() {
        // `Term` addition folds `i⁰ · 1` into `c₀`; only `Sum += Prod` does not.
        let mut t = Term::from(1.0) + Var(0).sin();
        t += Term::imaginary_unit() * Term::imaginary_unit();
        assert_eq!(t.n_monomials(), 1);
        assert!(matches!(t.inner(), Inner::Sum(s) if s.c0() == 0.0));
    }

    #[test]
    fn adding_a_constant_to_a_single_monomial_keeps_both() {
        let mut t = Var(0).sin();
        t += 2.0;
        assert_eq!(t.eval(&[1.1]).unwrap(), 2.0 + 1.1f64.sin());
    }

    #[test]
    fn subtraction_matches_adding_the_negation() {
        let a = Term::from(1.0) + Var(0).sin() * 2.0;
        let b = Var(0).sin() + Var(1).cos();
        let want = a.clone() + -b.clone();
        assert_eq!(a.clone() - b.clone(), want);
        assert_eq!(a.clone() - &b, want);
        let mut c = a.clone();
        c -= &b;
        assert_eq!(c, want);
        assert_eq!(a.clone() - 1.5, a + -1.5);
    }

    #[test]
    fn adding_a_term_ignores_the_right_operands_thresholds() {
        let mut big = Term::from(0.0) + Var(0).sin();
        big.set_min_eps(0.5);
        let mut small = Term::from(1.0) + Var(1).sin() * 1e-6;
        small.set_min_eps(1e-12);
        let mut l = small.clone();
        l += big.clone();
        assert_eq!(l.n_monomials(), 2, "left `min_eps` = 1e-12 keeps both");
        big += small;
        assert_eq!(
            big.n_monomials(),
            1,
            "left `min_eps` = 0.5 drops the 1e-6 monomial"
        );
    }
}
