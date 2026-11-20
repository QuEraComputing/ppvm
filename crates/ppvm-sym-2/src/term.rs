// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The user-facing [`Term`], its representation [`Inner`] and the angle [`Var`].

use crate::prod::Prod;
use crate::sum::Sum;

/// The `max_sin` of a new [`Term`]: no sine-degree truncation.
pub const DEFAULT_MAX_SIN: usize = usize::MAX;
/// The `min_eps` of a new [`Term`].
pub const DEFAULT_MIN_EPS: f64 = f64::EPSILON;

/// The representation of a [`Term`]. `Sum` alone could hold every value; `One`
/// and `Const` exist only for performance, keeping the common single-monomial
/// and scalar cases out of a hash map. `One` never holds the identity monomial.
///
/// `PartialEq` compares `c₀` and the monomials, not the variant, so `Const(c)`
/// equals a `Sum` holding only `c₀ = c`, and `One(p, c)` a `Sum` holding only `c · p`.
#[derive(Debug, Clone)]
pub(crate) enum Inner {
    /// A general sum of products.
    Sum(Sum),
    /// A single weighted product.
    One(Prod, f64),
    /// A numeric constant.
    Const(f64),
}

/// `Const(0.0)`, the placeholder while an operator rebuilds the representation.
impl Default for Inner {
    fn default() -> Self {
        Inner::Const(0.0)
    }
}

impl PartialEq for Inner {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Inner::Sum(a), Inner::Sum(b)) => a == b,
            // At most one side is a `Sum`; put it on the right, so the left
            // has at most one monomial to look up.
            _ => {
                let (a, b) = match self {
                    Inner::Sum(_) => (other, self),
                    _ => (self, other),
                };
                a.c0() == b.c0()
                    && a.n_monomials() == b.n_monomials()
                    && match a {
                        Inner::One(p, c) => b.coeff_of(p) == Some(*c),
                        _ => true,
                    }
            }
        }
    }
}

impl Inner {
    fn c0(&self) -> f64 {
        match self {
            Inner::Const(c) => *c,
            Inner::Sum(s) => s.c0,
            Inner::One(..) => 0.0,
        }
    }

    fn n_monomials(&self) -> usize {
        match self {
            Inner::Sum(s) => s.len(),
            Inner::One(..) => 1,
            Inner::Const(_) => 0,
        }
    }

    fn coeff_of(&self, p: &Prod) -> Option<f64> {
        match self {
            Inner::Sum(s) => s.maps.as_ref()?.terms.get(p).copied(),
            Inner::One(q, c) => (q == p).then_some(*c),
            Inner::Const(_) => None,
        }
    }

    /// `c · p` with the truncation of [`Sum::add_term`]: `Const(c)` for the
    /// identity, `Const(0.0)` when `p` fails a threshold, otherwise `One`.
    #[inline]
    pub(crate) fn truncated(p: Prod, c: f64, max: usize, min_eps: f64) -> Self {
        if p.is_identity() {
            Inner::Const(c)
        } else if p.sin_pow() > max || c.abs() <= min_eps {
            Inner::Const(0.0)
        } else {
            Inner::One(p, c)
        }
    }
}

/// A symbolic polynomial in `sin(x_i)` and `cos(x_i)`.
///
/// Carries its truncation settings, `max_sin` (maximum sine degree) and
/// `min_eps` (largest dropped coefficient magnitude), which every operation applies to its
/// result with the left operand's settings.
///
/// # Examples
///
/// ```
/// use ppvm_sym_2::Var;
///
/// // sin²(x0) at x0 = π/2 equals 1.
/// let expr = Var(0).sin() * Var(0).sin();
/// let v = expr.eval(&[std::f64::consts::FRAC_PI_2]).unwrap();
/// assert!((v - 1.0).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Term {
    pub(crate) inner: Inner,
    pub(crate) max_sin: usize, // max sin power
    pub(crate) min_eps: f64,   // coefficients with |c| <= min_eps are dropped
}

impl Term {
    /// Set the maximum sine power, dropping the monomials above it.
    #[inline]
    pub fn set_max_sin(&mut self, max: usize) {
        self.max_sin = max;
        self.retruncate();
    }

    /// Set the coefficient cutoff, dropping the monomials with `|c| <= eps`.
    #[inline]
    pub fn set_min_eps(&mut self, eps: f64) {
        self.min_eps = eps;
        self.retruncate();
    }

    /// The maximum sine power retained during arithmetic ([`DEFAULT_MAX_SIN`] by default).
    #[inline]
    pub fn max_sin(&self) -> usize {
        self.max_sin
    }

    /// The coefficient cutoff used during arithmetic ([`DEFAULT_MIN_EPS`] by default).
    #[inline]
    pub fn min_eps(&self) -> f64 {
        self.min_eps
    }

    /// The value as an `f64` if it has no monomials (`Const`, or a `Sum` with only `c₀`).
    #[inline]
    pub fn as_const(&self) -> Option<f64> {
        match &self.inner {
            Inner::Const(c) => Some(*c),
            Inner::Sum(s) if s.is_empty() => Some(s.c0),
            _ => None,
        }
    }

    /// The internal representation, read-only, for tests and diagnostics.
    #[cfg(test)]
    #[inline]
    pub(crate) fn inner(&self) -> &Inner {
        &self.inner
    }

    /// The angle variable `x_u`; shorthand for [`Var`]`(u)`.
    pub fn var(u: u32) -> Var {
        Var(u)
    }

    /// Build a constant term from `c`, with the default truncation settings.
    pub fn from_f64(c: f64) -> Self {
        Self::default_with(Inner::Const(c))
    }

    pub(crate) fn default_with(inner: Inner) -> Self {
        Self {
            inner,
            max_sin: DEFAULT_MAX_SIN,
            min_eps: DEFAULT_MIN_EPS,
        }
    }

    /// Monomials in this term: `0` for a constant, `1` for a single product,
    /// otherwise the table size.
    pub fn n_monomials(&self) -> usize {
        self.inner.n_monomials()
    }

    /// A term with this term's truncation parameters and a new representation.
    #[inline]
    pub(crate) fn with_inner(&self, inner: Inner) -> Self {
        Term {
            inner,
            max_sin: self.max_sin,
            min_eps: self.min_eps,
        }
    }
}

/// A symbolic angle `x_u`: the argument of `sin`/`cos`, never a [`Term`] itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Var(pub u32);

impl Var {
    /// `sin(x_u)`, with the default truncation settings.
    pub fn sin(self) -> Term {
        Term::default_with(Inner::One(Prod::sin(self.0), 1.0))
    }

    /// `cos(x_u)`, with the default truncation settings.
    pub fn cos(self) -> Term {
        Term::default_with(Inner::One(Prod::cos(self.0), 1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn setting_thresholds_prunes_but_keeps_the_constant() {
        let mut t = Term::from_f64(1e-9) + Var(0).sin() * Var(1).sin() + Var(1).cos() * 1e-6;
        t.set_max_sin(0);
        assert_eq!(t.n_monomials(), 1);
        t.set_min_eps(1e-3);
        assert_eq!((t.n_monomials(), t.as_const()), (0, Some(1e-9)));
        let mut one = Var(0).sin();
        one.set_max_sin(0);
        assert_eq!(one.inner(), &Inner::Const(0.0));
    }
}
