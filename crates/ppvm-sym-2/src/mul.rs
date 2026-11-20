// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The ring product on [`Prod`], [`Sum`] and [`Term`]. Operators take their
//! right operand by value so a `Sum` operand can be moved in instead of copied.

use crate::term::{Factor, Inner, Prod, Sum, Term, merge_factors};

impl Prod {
    /// Multiply this product by an additional `sin(x_u)`.
    #[inline(always)]
    pub(crate) fn mul_sin(&mut self, u: u32) {
        self.merge_factor(Factor::new(u, 1, 0));
        self.sin_pow += 1;
        self.debug_check();
    }

    /// Multiply this product by an additional `cos(x_u)`.
    #[inline(always)]
    pub(crate) fn mul_cos(&mut self, u: u32) {
        self.merge_factor(Factor::new(u, 0, 1));
        self.debug_check();
    }
}

impl Sum {
    /// Multiply in place by `coeff · p`, with the truncation of [`Sum::add_term`].
    ///
    /// If `p.sin_pow() > max` every product exceeds `max`, so the table is
    /// cleared without walking it (leaving an empty `Sum`, not `Const(0.0)`).
    #[inline]
    pub(crate) fn mul_term(&mut self, p: &Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max {
            if let Some(maps) = &mut self.maps {
                maps.terms.clear();
            }
            self.c0 = 0.0;
            return;
        }

        // `i` is not a scalar even though it has no factors.
        if p.is_identity() {
            self.scale(coeff);
            self.prune(max, min_eps);
            return;
        }

        // Ping-pong: `aux` is empty on entry, so after the swap `terms` holds
        // `aux`'s retained capacity and `aux` holds the table to rebuild from.
        let mut maps = self.maps.take().unwrap_or_default();
        debug_assert!(maps.aux.is_empty());
        std::mem::swap(&mut maps.terms, &mut maps.aux);
        let mut old_terms = std::mem::take(&mut maps.aux);
        maps.terms.reserve(old_terms.len());
        self.maps = Some(maps);

        // Reset `c₀` *before* the rebuild: `i · i` lands back on the identity
        // monomial and folds into `c₀`, which a reset afterwards would wipe.
        let c0 = std::mem::take(&mut self.c0) * coeff;
        if c0.abs() > min_eps {
            self.add_term(p.clone(), c0, max, min_eps);
        }
        for (mut k, v) in old_terms.drain() {
            let sign = k.mul_prod(p);
            self.add_term(k, v * coeff * sign, max, min_eps);
        }
        // Hand the (now empty, still-allocated) buffer back for the next multiply.
        self.maps.as_mut().unwrap().aux = old_terms;
    }

    /// Multiply every coefficient, `c₀` included, by `r`, without truncating.
    #[inline]
    pub(crate) fn scale(&mut self, r: f64) {
        self.c0 *= r;
        for v in self.maps.iter_mut().flat_map(|m| m.terms.values_mut()) {
            *v *= r;
        }
    }
}

impl Prod {
    /// Multiply in place by `rhs`, returning the sign (`±1.0`) the coefficient
    /// picks up: `−1` when both factors carry an `i`.
    #[inline]
    pub(crate) fn mul_prod(&mut self, rhs: &Prod) -> f64 {
        self.sin_pow += rhs.sin_pow;
        match (self.factors.len(), rhs.factors.as_slice()) {
            (_, []) => {}
            (0, _) => self.factors = rhs.factors.clone(),
            // The common case on the propagation hot path: the right operand is
            // a single `sin(x_u)` or `cos(x_u)` atom.
            (_, &[f]) => self.merge_factor(f),
            (_, fs) => self.factors = merge_factors(&self.factors, fs),
        }
        self.debug_check();
        self.compose_imag(rhs)
    }

    /// `self · rhs` as a new monomial plus the coefficient sign; see
    /// [`Prod::mul_prod`]. Allocates exactly one factor vector.
    #[inline]
    pub(crate) fn times(&self, rhs: &Prod) -> (Prod, f64) {
        let mut out = Prod {
            factors: merge_factors(&self.factors, &rhs.factors),
            sin_pow: self.sin_pow + rhs.sin_pow,
            imag: self.imag,
        };
        out.debug_check();
        let sign = out.compose_imag(rhs);
        (out, sign)
    }

    #[inline]
    fn compose_imag(&mut self, rhs: &Prod) -> f64 {
        self.set_phase_exponent(self.imag as u8 + rhs.imag as u8)
    }
}

impl Term {
    /// Flip the sign of every coefficient in place (exact; never truncates).
    #[inline]
    pub(crate) fn negate(&mut self) {
        match self.inner {
            Inner::Sum(ref mut s) => s.scale(-1.0),
            Inner::One(_, ref mut c) | Inner::Const(ref mut c) => *c = -*c,
        }
    }

    /// Whether `rhs` truncates less than `self`, so a representation taken over
    /// from it must be truncated again.
    #[inline]
    pub(crate) fn is_looser(&self, rhs: &Term) -> bool {
        rhs.max_sin > self.max_sin || rhs.min_eps < self.min_eps
    }

    /// Apply `self`'s settings to its current representation.
    pub(crate) fn retruncate(&mut self) {
        let (max, eps) = (self.max_sin, self.min_eps);
        match &mut self.inner {
            Inner::Sum(s) => s.prune(max, eps),
            Inner::One(..) => {
                let Inner::One(p, c) = std::mem::take(&mut self.inner) else {
                    unreachable!()
                };
                self.inner = Inner::truncated(p, c, max, eps);
            }
            Inner::Const(_) => {}
        }
    }
}

/// Scales every coefficient; shrinking (`|rhs| < 1`) prunes monomials at or below
/// `min_eps`, so a sign flip stays a single pass.
impl std::ops::MulAssign<f64> for Term {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        match self.inner {
            Inner::Sum(ref mut s) => {
                s.scale(rhs);
                if rhs.abs() < 1.0 {
                    s.prune(self.max_sin, self.min_eps);
                }
            }
            Inner::One(_, ref mut c) => {
                *c *= rhs;
                if c.abs() <= self.min_eps {
                    self.inner = Inner::Const(0.0);
                }
            }
            Inner::Const(ref mut c) => *c *= rhs,
        }
    }
}

/// Uses `self`'s truncation settings; `rhs`'s are ignored. Only an owned `Sum`
/// operand is moved in; everything else goes through the borrowed product.
impl std::ops::MulAssign<Term> for Term {
    #[inline]
    fn mul_assign(&mut self, rhs: Term) {
        let (max, eps) = (self.max_sin, self.min_eps);
        let looser = self.is_looser(&rhs);
        let rhs_settings = (rhs.max_sin, rhs.min_eps);
        match (&mut self.inner, rhs.inner) {
            (Inner::Const(c), Inner::Sum(s)) => {
                let c = *c;
                self.inner = Inner::Sum(s);
                if looser {
                    self.retruncate();
                }
                *self *= c;
            }
            (Inner::One(p, c), Inner::Sum(mut s)) => {
                s.mul_term(p, *c, max, eps);
                self.inner = Inner::Sum(s);
            }
            (_, inner) => {
                let mut rhs = Term::default_with(inner);
                (rhs.max_sin, rhs.min_eps) = rhs_settings;
                *self *= &rhs;
            }
        }
    }
}

/// Never clones `rhs` unless the result must own its table (`Const`/`One` × `Sum`).
impl std::ops::MulAssign<&Term> for Term {
    #[inline]
    fn mul_assign(&mut self, rhs: &Term) {
        let (max, eps) = (self.max_sin, self.min_eps);
        let looser = self.is_looser(rhs);
        let inner = match (&mut self.inner, &rhs.inner) {
            (_, Inner::Const(c)) => return *self *= *c,
            (Inner::Const(c), rhs) => {
                let c = *c;
                self.inner = rhs.clone();
                if looser {
                    self.retruncate();
                }
                return *self *= c;
            }
            (Inner::Sum(s), Inner::One(p, c)) => return s.mul_term(p, *c, max, eps),
            (Inner::One(p, c), Inner::Sum(s)) => {
                let mut s = s.clone();
                s.mul_term(p, *c, max, eps);
                Inner::Sum(s)
            }
            (Inner::One(p, c), Inner::One(p2, c2)) => {
                let (p, sign) = p.times(p2);
                Inner::truncated(p, *c * c2 * sign, max, eps)
            }
            // `c₀` takes part as the identity monomial, which `add_term` folds back.
            (Inner::Sum(s1), Inner::Sum(s2)) => {
                let mut out = Sum::new();
                let rhs: Vec<_> = s2.parts().collect();
                for (p1, c1) in s1.parts() {
                    for &(p2, c2) in &rhs {
                        let (p, sign) = p1.times(p2);
                        out.add_term(p, c1 * c2 * sign, max, eps);
                    }
                }
                Inner::Sum(out)
            }
        };
        self.inner = inner;
    }
}

binary_ops!(Mul::mul from MulAssign::mul_assign for Term, &Term, f64);

impl std::ops::Mul<Term> for f64 {
    type Output = Term;

    #[inline]
    fn mul(self, rhs: Term) -> Term {
        rhs * self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::Var;
    use crate::term::{DEFAULT_MAX_SIN, DEFAULT_MIN_EPS};

    #[test]
    fn monomial_canonical_under_multiplication_order() {
        // Same monomial built by different multiplication orders must be the
        // SAME key (canonicality is what stops the term count exploding).
        let mut a = Prod::sin(3);
        a.mul_prod(&Prod::cos(1));
        a.mul_prod(&Prod::sin(1));

        let mut b = Prod::sin(1);
        b.mul_prod(&Prod::sin(3));
        b.mul_prod(&Prod::cos(1));
        assert_eq!(a, b.times(&Prod::new()).0);

        assert_eq!(a, b);
        assert_eq!(hash_of(&a), hash_of(&b));
        assert_eq!(a.sin_pow(), 2);
        assert_eq!(a.cos_pow(), 1);
    }

    fn hash_of(p: &Prod) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = fxhash::FxHasher64::default();
        p.hash(&mut h);
        h.finish()
    }

    #[test]
    fn i_times_i_is_a_sign_on_the_real_key() {
        let mut a = Prod::sin(0);
        assert_eq!(a.mul_phase(1), 1.0);
        let mut b = Prod::cos(0);
        assert_eq!(b.mul_phase(3), -1.0); // i³ = −i: key `i`, sign −1
        let (c, sign) = a.times(&b); // the keys multiply as i · i = −1
        assert!(!c.is_imag());
        assert_eq!(sign, -1.0);
    }

    #[test]
    fn mul_term_zero_shortcut_leaves_empty_sum() {
        // The shortcut leaves an EMPTY `Sum`, not `Const(0.0)`.
        let mut s = Sum::new();
        s.add_term(Prod::cos(0), 1.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        let mut p = Prod::sin(0);
        p.mul_sin(1);
        s.mul_term(&p, 1.0, 1, f64::EPSILON);
        assert!(s.is_empty());
        assert_eq!(s.c0(), 0.0);
    }

    #[test]
    fn mul_term_reuses_the_aux_allocation() {
        // After a rebuild both buffers are still allocated.
        let mut s = Sum::new();
        for i in 0..16u32 {
            s.add_term(Prod::cos(i), 1.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        }
        s.mul_term(&Prod::sin(0), 1.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        let maps = s.maps.as_ref().expect("first monomial creates the map box");
        assert!(maps.aux.is_empty());
        assert!(maps.aux.capacity() > 0, "aux allocation was freed");
        assert_eq!(s.len(), 16);
    }

    /// `One × Sum` keeps the product (owned and borrowed), and a negative `c₀`
    /// takes part in `Sum × Sum`.
    #[test]
    fn products_with_sums_keep_every_cross_term() {
        let x = [1.1f64, 2.2];
        let (s0, s1, c1) = (x[0].sin(), x[1].sin(), x[1].cos());
        let a = Term::from_f64(-1.0) + Var(0).sin();
        let b = Var(1).sin() + Var(1).cos();
        let close = |t: Term, want: f64| (t.eval(&x).unwrap() - want).abs() < 1e-12;
        assert!(close(Var(0).sin() * b.clone(), s0 * (s1 + c1)));
        assert!(close(Var(0).sin() * &b, s0 * (s1 + c1)));
        assert!(close(a * b, (s0 - 1.0) * (s1 + c1)));
    }

    #[test]
    fn products_of_small_forms_stay_out_of_the_map() {
        let mut capped = Var(0).sin();
        capped.set_max_sin(1);
        assert_eq!((capped * Var(1).sin()).inner(), &Inner::Const(0.0));
        let t = Var(0).sin() * Var(1).cos();
        assert!(matches!(t.inner(), Inner::One(..)));
        let t = Term::from_f64(2.0) * Var(0).sin();
        assert!(matches!(t.inner(), Inner::One(_, c) if *c == 2.0));
    }

    use crate::coeff::tests::forms;
    use ppvm_traits_2::{Conjugate, ImaginaryUnit};
    use std::collections::HashMap;

    fn sum_term(max_sin: usize) -> Term {
        Term {
            inner: Inner::Sum(Sum::new()),
            max_sin,
            min_eps: f64::EPSILON,
        }
    }

    /// The nonzero monomials of a term, `c₀` under the identity key.
    fn table(t: &Term) -> HashMap<Prod, f64> {
        let mut s = sum_term(usize::MAX);
        s += t.clone();
        let Inner::Sum(s) = s.inner else {
            unreachable!()
        };
        s.parts()
            .filter(|(_, c)| *c != 0.0)
            .map(|(p, c)| (p.clone(), c))
            .collect()
    }

    type BinOp = fn(Term, &Term) -> Term;

    /// With `max_sin = k` on the left operand, in any representation, a product,
    /// sum or difference is exactly the untruncated result with every monomial of
    /// sine degree `> k` removed; every other operation's result obeys `k` too.
    #[test]
    fn max_sin_truncation_is_exact_in_every_form() {
        let forms = forms();
        let obeys = |t: &Term, k| table(t).keys().all(|p| p.sin_pow() <= k);
        for k in 0..=4 {
            for a in &forms {
                let mut lhs = a.clone();
                lhs.set_max_sin(k);
                assert!(obeys(&lhs, k), "k = {k}: {a}");
                for t in [
                    lhs.mul_phase(1),
                    lhs.conj(),
                    lhs.clone() * 0.5,
                    -lhs.clone(),
                ] {
                    assert!(obeys(&t, k), "k = {k}: {a} → {t}");
                }
                for b in &forms {
                    let ops: [(BinOp, &str); 5] = [
                        (|x, y| x * y.clone(), "*"),
                        (|x, y| x * y, "*&"),
                        (|x, y| x + y.clone(), "+"),
                        (|x, y| x + y, "+&"),
                        (|x, y| x - y, "-&"),
                    ];
                    for (op, name) in ops {
                        let full = table(&op(a.clone(), b));
                        let got = table(&op(lhs.clone(), b));
                        let want: HashMap<_, _> =
                            full.into_iter().filter(|(p, _)| p.sin_pow() <= k).collect();
                        assert_eq!(got.len(), want.len(), "k = {k}: {a} {name} {b}");
                        for (p, c) in &want {
                            assert!((got[p] - c).abs() < 1e-12, "k = {k}: {a} {name} {b} at {p}");
                        }
                    }
                }
            }
        }
    }

    /// `min_eps` applies to monomials in every form; the constant is never dropped.
    #[test]
    fn min_eps_applies_to_monomials_only() {
        let mut one = Var(0).sin();
        one.set_min_eps(1e-3);
        assert_eq!((one.clone() * 1e-4).inner(), &Inner::Const(0.0));
        assert_eq!(
            (one.clone() * (Var(1).cos() * 1e-4)).inner(),
            &Inner::Const(0.0)
        );
        let mut c = Term::from(1.0);
        c.set_min_eps(1e-3);
        assert_eq!((c.clone() * 1e-9).as_const(), Some(1e-9));
        // A table adopted from a looser right operand is pruned to the left's settings.
        let loose = Term::from(1.0) + Var(0).sin() * 1e-6;
        assert_eq!((c.clone() * loose.clone()).n_monomials(), 0);
        assert_eq!((c.clone() * &loose).n_monomials(), 0);
        assert_eq!((c + loose).n_monomials(), 0);
    }

    /// `min_eps` drops `|c| <= min_eps`, so `min_eps = 0` drops exactly the zeros.
    #[test]
    fn min_eps_drops_coefficients_at_the_threshold() {
        let eps = 1e-3;
        let mut s = Sum::new();
        s.add_term(Prod::sin(0), eps * 4.0, usize::MAX, eps);
        s.add_term(Prod::sin(1), eps, usize::MAX, eps);
        assert_eq!(s.len(), 1);
        s.mul_term(&Prod::cos(0), 0.5, usize::MAX, eps);
        assert_eq!(s.len(), 1, "the product stays above `min_eps`");
        s.mul_term(&Prod::cos(0), 0.5, usize::MAX, eps);
        assert!(s.is_empty(), "the product reached `min_eps`");

        let x = Var(0).sin() + Var(1).sin();
        let mut t = x.clone();
        t.set_min_eps(0.0);
        assert_eq!(
            (t.clone() - &x).n_monomials(),
            0,
            "cancellation removes the entry"
        );
        assert_eq!((t.clone() * 0.0).n_monomials(), 0);
        assert!(num::Zero::is_zero(&(t * 0.0)));
    }

    #[test]
    fn mul_term_on_an_unallocated_sum_moves_c0_into_the_table() {
        let mut s = Sum::new();
        s.c0 = 2.0;
        s.mul_term(&Prod::sin(3), 0.5, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        assert_eq!(s.c0(), 0.0);
        assert_eq!(s.iter().collect::<Vec<_>>(), vec![(&Prod::sin(3), 1.0)]);
    }

    #[test]
    fn product_of_constant_only_sums_stays_constant() {
        let mut a = sum_term(usize::MAX);
        a += 2.0;
        let mut b = sum_term(usize::MAX);
        b += -3.0;
        let ab = a * b;
        assert!(matches!(ab.inner(), Inner::Sum(s) if s.c0() == -6.0 && s.is_empty()));
    }

    /// A long chain of in-place multiplies, cloning along the way, tracks the
    /// running product in `ℂ`: the `aux` ping-pong never leaks stale monomials.
    #[test]
    fn long_multiply_chains_track_the_running_product() {
        let mut vals = vec![0.0; 41];
        vals[0] = 0.4;
        vals[1] = -1.3;
        vals[40] = 2.2;
        let ev = |t: &Term| t.eval_complex(&vals).unwrap();
        let factors = [
            Var(0).sin(),
            Term::from(0.5) + Var(1).cos(),
            Term::imaginary_unit(),
            Var(40).cos() * 1.5,
            Var(0).cos().mul_phase(3) + Var(1).sin(),
            Term::from(-2.0),
        ];
        let mut t = Term::from(0.25) + Var(1).sin();
        let mut want = ev(&t);
        for step in 0..60 {
            let f = &factors[step % factors.len()];
            want *= ev(f);
            t *= f;
            if step % 7 == 0 {
                t = t.clone();
            }
            let got = ev(&t);
            assert!(
                (got - want).norm() <= 1e-9 * want.norm().max(1.0),
                "step {step}"
            );
            if let Inner::Sum(s) = t.inner() {
                assert!(s.maps.as_ref().is_none_or(|m| m.aux.is_empty()));
            }
        }
        assert!(want.norm() > 0.0);
    }
}
