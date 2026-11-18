// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Term`] as a `ppvm-traits-2` coefficient: [`Coefficient`], [`Angle`],
//! [`ImaginaryUnit`], [`Conjugate`] and `Halvable`.

use ppvm_traits_2::{Angle, Coefficient, Conjugate, ImaginaryUnit};

use crate::term::{Inner, ONE, Prod, Sum, Term, Var};

impl Coefficient for Term {
    /// Flips the stored signs in place (or does nothing for `+1`); never a
    /// `Term` product. Panics on a bare variable.
    #[inline]
    fn mul_sign(&self, sign: i8) -> Self {
        let mut ret = self.clone();
        ret.mul_sign_assign(sign);
        ret
    }

    #[inline]
    fn mul_sign_assign(&mut self, sign: i8) {
        if sign < 0 {
            self.negate();
        }
    }

    /// `|c|` for a constant (including a `Sum` with no monomials), `+∞` otherwise.
    ///
    /// So a `CoefficientThreshold` policy never drops a symbolic coefficient;
    /// symbolic terms truncate through `max_sin`/`min_eps` instead. This breaks the
    /// law `N(x) == 0 iff x == 0`, but no absolute value exists on this ring, and an
    /// `ℓ¹` norm (a usable substitute) would start dropping symbolic coefficients.
    fn magnitude(&self) -> f64 {
        self.as_const().map_or(f64::INFINITY, f64::abs)
    }
}

/// `x_u ↦ (sin x_u, cos x_u)`. The ring is the free polynomial ring, so
/// `sin² + cos² = 1` holds only after evaluation.
impl Angle<Term> for Var {
    #[inline]
    fn sin_cos(&self) -> (Term, Term) {
        (self.sin(), self.cos())
    }
}

/// An `f64` angle on a symbolic sum, folded to two constants.
impl Angle<Term> for f64 {
    #[inline]
    fn sin_cos(&self) -> (Term, Term) {
        let (s, c) = f64::sin_cos(*self);
        (Term::from_f64(s), Term::from_f64(c))
    }
}

impl Term {
    /// Multiply by `i^phase`: `eval_complex(x.mul_phase(k)) == iᵏ · eval_complex(x)`.
    /// Even phases scale in place; odd ones flip every key's `i` bit and rebuild.
    #[inline]
    pub fn mul_phase(&self, phase: u8) -> Self {
        match (&self.inner, phase & 3) {
            (_, 0) => self.clone(),
            (_, 2) => -self.clone(),
            (_, k) => self.rotate_phase(k),
        }
    }

    /// [`Term::mul_phase`] for odd `k`.
    fn rotate_phase(&self, k: u8) -> Self {
        let relabel = |p: &Prod, c: f64| {
            let mut q = p.clone();
            let sign = q.mul_phase(k);
            (q, c * sign)
        };
        match self.inner {
            Inner::Sum(ref s) => {
                let mut ret = Sum::new();
                for (p, c) in s.iter().chain([(&ONE, s.c0)]) {
                    let (q, c) = relabel(p, c);
                    ret.add_term(q, c, self.max_sin, self.min_eps);
                }
                self.with_inner(Inner::Sum(ret))
            }
            Inner::Const(c) => {
                let (q, c) = relabel(&ONE, c);
                self.with_inner(Inner::truncated(q, c, self.max_sin, self.min_eps))
            }
            Inner::One(ref p, c) => {
                let (q, c) = relabel(p, c);
                self.with_inner(Inner::truncated(q, c, self.max_sin, self.min_eps))
            }
        }
    }
}

/// `i` is `One(i, 1.0)`; `i·i` normalises back to `Const(−1.0)`.
impl ImaginaryUnit for Term {
    fn imaginary_unit() -> Self {
        Term::from_f64(1.0).mul_phase(1)
    }

    /// Folds `i` into the monomials rather than running a `Term` product.
    fn mul_i(&self) -> Self {
        self.mul_phase(1)
    }

    /// Folds `iᵏ` into the monomials.
    fn mul_i_pow(&self, k: u8) -> Self {
        self.mul_phase(k & 3)
    }
}

/// `i ↦ −i`: negates every imaginary monomial in place, since no key changes.
impl Conjugate for Term {
    fn conj(&self) -> Self {
        let mut ret = self.clone();
        match ret.inner {
            Inner::Sum(ref mut s) => {
                if let Some(maps) = &mut s.maps {
                    for (_, c) in maps.terms.iter_mut().filter(|(p, _)| p.imag) {
                        *c = -*c;
                    }
                }
            }
            Inner::One(ref p, ref mut c) if p.imag => *c = -*c,
            Inner::One(..) | Inner::Const(_) => {}
        }
        ret
    }
}

/// Folds from the first summand, so its truncation settings apply, as for `a + b
/// + …`; an empty iterator gives `Const(0.0)` with the default settings.
impl std::iter::Sum for Term {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let mut total = iter.next().unwrap_or_else(|| Term::from_f64(0.0));
        for t in iter {
            total += t;
        }
        total
    }
}

/// Exact in `f64`: scales the coefficients by `0.5` without pruning, so
/// `x.half() + x.half() == x`.
impl ppvm_traits_2::Halvable for Term {
    #[inline]
    fn half(&self) -> Self {
        let mut ret = self.clone();
        match ret.inner {
            Inner::Sum(ref mut s) => s.scale(0.5),
            Inner::One(_, ref mut c) | Inner::Const(ref mut c) => *c *= 0.5,
        }
        ret
    }
}

impl std::ops::Neg for Term {
    type Output = Term;

    #[inline]
    fn neg(mut self) -> Self::Output {
        self.negate();
        self
    }
}

macro_rules! from_number {
    ($($t:ty),*) => {$(
        impl From<$t> for Term {
            #[inline]
            fn from(value: $t) -> Self {
                Term::from_f64(value as f64)
            }
        }
    )*};
}
from_number!(f32, f64, i32, i64);

/// `is_zero()` holds for a constant (`Const`, or a `Sum` with no monomials) below `min_eps`.
impl num::Zero for Term {
    fn zero() -> Self {
        Term::from_f64(0.0)
    }

    fn is_zero(&self) -> bool {
        self.as_const().is_some_and(|c| c.abs() < self.min_eps)
    }
}

/// `one()` is `Const(1.0)`.
impl num::One for Term {
    fn one() -> Self {
        Term::from_f64(1.0)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::term::Var;
    use num::{Complex, Zero};

    /// Constants in any form are zero below `min_eps` and have `|c|` as their
    /// magnitude; symbolic forms never are zero and never get dropped by a threshold.
    #[test]
    fn constants_in_any_form_are_zero_and_have_a_magnitude() {
        assert!(Term::from_f64(0.0).is_zero());
        assert_eq!(Term::from_f64(-1e-30).magnitude(), 1e-30);
        let mut t = Var(0).sin();
        t += Var(0).sin() * -1.0;
        assert!(matches!(t.inner(), Inner::Sum(_)));
        assert!(t.is_zero(), "a cancelled `Sum` is zero");
        assert_eq!(t.magnitude(), 0.0);
        t += -2.5;
        assert!(!t.is_zero());
        assert_eq!(t.magnitude(), 2.5);
        let mut tiny = Var(0).sin();
        tiny.set_min_eps(0.0);
        tiny *= 1e-30;
        assert!(!tiny.is_zero());
        assert_eq!(tiny.magnitude(), f64::INFINITY);
        // With the default `min_eps` the product truncates to an exact zero.
        assert!((Var(0).sin() * 1e-30).is_zero());
    }

    #[test]
    fn equality_looks_through_the_representation() {
        let with = |inner| Term {
            inner,
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        };
        let mut c0_only = Sum::new();
        c0_only.c0 = 1.0;
        assert_eq!(with(Inner::Sum(c0_only.clone())), Term::from(1.0));
        assert_eq!(Term::from(1.0), with(Inner::Sum(c0_only)));
        assert_eq!(with(Inner::Sum(Sum::new())), Term::from(0.0));

        let x = Var(0).sin() * 2.0;
        let mut one_entry = Term::from(0.0) + x.clone();
        assert!(matches!(one_entry.inner(), Inner::Sum(_)));
        assert_eq!(one_entry, x);
        assert_eq!(x, one_entry);
        assert_ne!(one_entry, x.mul_phase(1), "the `i` bit is part of the key");
        one_entry += 1.0;
        assert_ne!(one_entry, x, "an extra c₀ differs");
        assert_ne!(Var(0).sin() * 3.0, x, "coefficients differ");
        assert_ne!(Var(1).sin() * 2.0, x, "keys differ");

        // `x − x` is the zero constant in every form.
        let mut t = Var(0).sin() + Var(1).cos();
        t += (Var(0).sin() + Var(1).cos()) * -1.0;
        assert_eq!(t, Term::from(0.0));
        assert_eq!(
            Term::imaginary_unit() * Term::imaginary_unit(),
            -Term::from(1.0)
        );

        assert_ne!(
            Term::from(1.0),
            Term::from(1.0) + x.clone(),
            "extra monomials differ"
        );
        assert_ne!(Term::from(1.0) + x.clone(), Term::from(1.0));
        let mut tight = Term::from(1.0);
        tight.set_max_sin(2);
        assert_ne!(tight, Term::from(1.0), "truncation settings are compared");
    }

    #[test]
    fn half_is_exact_and_sum_keeps_the_first_settings() {
        use ppvm_traits_2::Halvable;
        for a in forms() {
            assert_eq!(a.half() + a.half(), a, "{a}");
        }
        let mut first = Term::from(1.0) + Var(0).sin();
        first.set_max_sin(1);
        let total: Term = [first, Var(1).sin(), Var(0).sin() * Var(1).sin()]
            .into_iter()
            .sum();
        assert_eq!((total.max_sin(), total.n_monomials()), (1, 2));
    }

    #[test]
    fn conversions_and_neg() {
        for t in [Term::from(2i32), Term::from(2i64), Term::from(2.0f32)] {
            assert_eq!(t, Term::from(2.0f64));
        }
        assert_eq!(-Term::from(2.0), Term::from(-2.0));
    }

    /// One of every representation and phase: `Const`, `One`, map-backed sums
    /// with real, imaginary and negative parts, and a variable past the eval cache.
    pub(crate) fn forms() -> Vec<Term> {
        let x = Var(0).sin();
        let s = Term::from(0.5) + Var(1).cos();
        let mut c0_only = Term::from(0.0) + Var(1).cos();
        c0_only *= 0.0;
        c0_only += 0.75;
        vec![
            Term::from(2.0),
            Term::from(0.0),
            x.clone(),
            x.clone() * x.clone(),
            s.clone(),
            s.clone() + x.clone(),
            Term::imaginary_unit(),
            x.mul_phase(1),
            x.mul_phase(3),
            s.mul_phase(3),
            s.mul_phase(2), // c₀ = −0.5
            (s.clone() + x.clone()).mul_phase(1) + x.clone(),
            c0_only,
            Var(40).cos() + Var(0).cos() * 2.0,
        ]
    }

    pub(crate) const PT: [f64; 41] = {
        let mut v = [0.0; 41];
        v[0] = 0.3;
        v[1] = 1.1;
        v[40] = -0.7;
        v
    };

    pub(crate) fn ev(t: &Term) -> Complex<f64> {
        t.eval_complex(&PT).unwrap()
    }

    pub(crate) fn close(a: Complex<f64>, b: Complex<f64>) -> bool {
        (a - b).norm() < 1e-12 * (1.0 + a.norm().max(b.norm()))
    }

    /// `eval_complex` is a ring homomorphism to `ℂ` commuting with `conj`,
    /// `mul_phase` and `mul_sign`; `eval` is its real part, and borrowed ops
    /// match owned ones.
    #[test]
    fn evaluation_is_a_ring_homomorphism() {
        let i = Term::imaginary_unit();
        let forms = forms();
        for a in &forms {
            let ea = ev(a);
            assert!((a.eval(&PT).unwrap() - ea.re).abs() < 1e-12, "eval of {a}");
            assert!(close(ev(&a.conj()), ea.conj()), "conj of {a}");
            assert!(close(ev(&(a.clone() + 1.5)), ea + 1.5), "{a} + 1.5");
            assert!(close(ev(&(1.5 + a.clone())), ea + 1.5), "1.5 + {a}");
            assert!(close(ev(&(a.clone() - 1.5)), ea - 1.5), "{a} - 1.5");
            assert!(close(ev(&(a.clone() * -3.0)), ea * -3.0), "{a} · -3");
            assert!(close(ev(&(-3.0 * a.clone())), ea * -3.0), "-3 · {a}");
            assert!(close(ev(&-a.clone()), -ea), "-{a}");
            assert_eq!(a.mul_sign(-1), -a.clone(), "{a}");
            assert_eq!(a.mul_sign(1), *a, "{a}");
            let mut b = a.clone();
            b.mul_sign_assign(-1);
            assert_eq!(b, a.mul_sign(-1));
            for k in 0..4u8 {
                assert!(close(ev(&a.mul_phase(k)), ea * Complex::i().powu(k as u32)));
            }
            assert_eq!(a.clone() * i.clone(), a.mul_phase(1), "{a} · i");
            for b in &forms {
                let eb = ev(b);
                assert!(close(ev(&(a.clone() + b.clone())), ea + eb), "{a} + {b}");
                assert!(close(ev(&(a.clone() - b.clone())), ea - eb), "{a} - {b}");
                assert!(close(ev(&(a.clone() * b.clone())), ea * eb), "{a} · {b}");

                let mut by_ref = a.clone();
                by_ref += b;
                assert_eq!(by_ref, a.clone() + b.clone(), "{a} += &{b}");
                let mut by_ref = a.clone();
                by_ref -= b;
                assert_eq!(by_ref, a.clone() - b.clone(), "{a} -= &{b}");
                let mut by_ref = a.clone();
                by_ref *= b;
                assert_eq!(by_ref, a.clone() * b.clone(), "{a} *= &{b}");
            }
        }
        assert_eq!(i.conj(), -i);
    }

    /// Associativity and distributivity hold in value for every triple of forms.
    #[test]
    fn ring_axioms_hold_after_evaluation() {
        let forms = forms();
        for a in &forms {
            for b in &forms {
                for c in &forms {
                    let (a, b, c) = (a.clone(), b.clone(), c.clone());
                    let abc = ev(&((a.clone() * b.clone()) * c.clone()));
                    assert!(close(abc, ev(&(a.clone() * (b.clone() * c.clone())))));
                    let lhs = ev(&(a.clone() * (b.clone() + c.clone())));
                    assert!(close(
                        lhs,
                        ev(&(a.clone() * b.clone() + a.clone() * c.clone()))
                    ));
                    let lhs = ev(&((a.clone() + b.clone()) + c.clone()));
                    assert!(close(lhs, ev(&(a + (b + c)))));
                }
            }
        }
    }

    #[test]
    fn conj_is_an_involution_and_fixes_real_forms() {
        for a in forms() {
            assert_eq!(a.conj().conj(), a, "{a}");
        }
        let x = Var(0).sin() + Var(1).cos().mul_phase(1);
        assert_eq!(x.conj().conj(), x);
        let real = Term::from(0.5) + Var(1).cos();
        assert_eq!(real.conj(), real);
        assert_eq!(Term::from(2.0).conj(), Term::from(2.0));
    }

    #[test]
    fn mul_phase_reduces_any_exponent_mod_4() {
        for a in [Var(0).sin(), Var(0).sin().mul_phase(1), Term::from(2.0)] {
            for k in 0..=u8::MAX {
                assert_eq!(a.mul_phase(k), a.mul_phase(k % 4), "{a}, k = {k}");
                assert_eq!(a.mul_i_pow(k), a.mul_phase(k % 4), "{a}, k = {k}");
            }
        }
        let mut p = Prod::sin(0);
        assert_eq!(p.mul_phase(1), 1.0);
        assert_eq!(p.mul_phase(255), 1.0); // i · i³ = 1
        assert!(!p.is_imag());
    }

    #[test]
    fn four_factors_of_i_are_the_identity() {
        let x = Var(0).sin() + Var(1).cos();
        let i = Term::imaginary_unit();
        let mut y = x.clone();
        for _ in 0..4 {
            y *= &i;
        }
        assert_eq!(y, x);
        assert_eq!(x.mul_phase(2).mul_phase(2), x);
    }

    #[test]
    fn iter_sum_matches_folded_addition() {
        let forms = forms();
        let total: Term = forms.iter().cloned().sum();
        let expected: Complex<f64> = forms.iter().map(ev).sum();
        assert!(close(ev(&total), expected));
        let empty: Term = std::iter::empty().sum();
        assert_eq!(empty, Term::from(0.0));
    }

    /// The identity monomial never appears as a key or inside `One`.
    #[test]
    fn identity_monomial_is_always_normalised() {
        let i = Term::imaginary_unit();
        assert!(matches!(i.inner(), Inner::One(p, c) if p.is_imag() && *c == 1.0));
        assert_eq!(i.clone() * i.clone(), -Term::from(1.0));
        assert_eq!(i.mul_phase(3), Term::from(1.0));
        assert_eq!(Term::from(2.0).mul_i_pow(0), Term::from(2.0));
        assert_eq!(Term::from(2.0).mul_i_pow(2), Term::from(-2.0));
        let no_identity_key = |t: &Term| match t.inner() {
            Inner::Sum(s) => s.iter().all(|(p, _)| !p.is_identity()),
            Inner::One(p, _) => !p.is_identity(),
            _ => true,
        };
        for a in forms() {
            for b in forms() {
                let ab = a.clone() * b.clone();
                assert!(no_identity_key(&ab), "{a} · {b}");
                assert!(no_identity_key(&(a.clone() + b.clone())), "{a} + {b}");
                assert!(no_identity_key(&ab.mul_phase(1)), "i · {a} · {b}");
            }
        }
    }

    #[test]
    fn opposite_phases_share_a_key_and_cancel() {
        let x = Var(0).sin();
        let t = x.mul_phase(1) + x.mul_phase(3); // i·x − i·x
        assert_eq!(t.n_monomials(), 0);
        assert_eq!(t.eval_complex(&[0.4]).unwrap(), Complex::new(0.0, 0.0));
        // i² is a sign on the real key, not a separate key.
        assert_eq!(x.mul_phase(2), -x);
    }

    #[test]
    fn sin_cos_is_the_angle_domain() {
        let (s, c) = Angle::<Term>::sin_cos(&Var(0));
        assert_eq!(s, Var(0).sin());
        assert_eq!(c, Var(0).cos());
        let (s, c) = Angle::<Term>::sin_cos(&0.3f64);
        assert_eq!(s, Term::from_f64(0.3f64.sin()));
        assert_eq!(c, Term::from_f64(0.3f64.cos()));
    }
}
