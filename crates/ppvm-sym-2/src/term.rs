// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The monomial [`Prod`], the monomial table [`Sum`] and the user-facing [`Term`].

use std::cmp::Ordering;
use std::collections::hash_map::Entry;

use fxhash::FxHashMap;

/// The `max_sin` of a new [`Term`]: no sine-degree truncation.
pub const DEFAULT_MAX_SIN: usize = usize::MAX;
/// The `min_eps` of a new [`Term`].
pub const DEFAULT_MIN_EPS: f64 = f64::EPSILON;

/// One variable's contribution to a monomial: `sin(x_var)^sin · cos(x_var)^cos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) struct Factor {
    /// The variable id.
    pub(crate) var: u32,
    /// The exponent of `sin(x_var)` (may be `0` if only `cos` is present).
    pub(crate) sin: u32,
    /// The exponent of `cos(x_var)` (may be `0` if only `sin` is present).
    pub(crate) cos: u32,
}

/// `[i] sin^m cos^n`: a monomial over symbolic variables.
///
/// Canonical, so the derived `Hash`/`Eq` identify the monomial: factors sorted
/// by variable, one per variable, none all-zero. A phase `iᵏ` keeps only its
/// `i` bit here; the sign `(−1)^⌊k/2⌋` goes into the coefficient, which is why
/// every phase operation returns a `±1.0`.
#[derive(Debug, Clone, Default)]
pub(crate) struct Prod {
    /// Ascending by `var`, one entry per variable, no all-zero entry.
    pub(crate) factors: Vec<Factor>,
    /// Total sine degree, cached for the `max_sin` check on every insert.
    pub(crate) sin_pow: usize,
    /// Whether the monomial carries a factor `i`.
    pub(crate) imag: bool,
}

/// `sin_pow` is derived from `factors`, so identity is `(factors, imag)`.
impl PartialEq for Prod {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.imag == other.imag && self.factors == other.factors
    }
}

impl Eq for Prod {}

impl std::hash::Hash for Prod {
    #[inline]
    fn hash<H: std::hash::Hasher>(&self, h: &mut H) {
        self.factors.hash(h);
        self.imag.hash(h);
    }
}

/// The real identity monomial `1`.
pub(crate) static ONE: Prod = Prod::new();

impl Factor {
    #[inline]
    pub(crate) const fn new(var: u32, sin: u32, cos: u32) -> Self {
        Self { var, sin, cos }
    }
}

impl Prod {
    /// Construct the empty product (the multiplicative identity, value `1`).
    pub(crate) const fn new() -> Self {
        Self {
            factors: Vec::new(),
            sin_pow: 0,
            imag: false,
        }
    }

    /// Multiply by `iᵏ`, returning the sign (`±1.0`) the coefficient must be
    /// multiplied by to keep the phase canonical.
    #[inline]
    pub(crate) fn mul_phase(&mut self, k: u8) -> f64 {
        self.set_phase_exponent(self.imag as u8 + (k & 3))
    }

    /// Whether the monomial carries a factor `i`.
    #[cfg(test)]
    #[inline]
    pub(crate) fn is_imag(&self) -> bool {
        self.imag
    }

    /// Whether this is the real identity monomial `1`, which is never a table key.
    #[inline]
    pub(crate) fn is_identity(&self) -> bool {
        self.factors.is_empty() && !self.imag
    }

    /// Set the phase to `iᵉ`: keep `e`'s odd part in the `imag` bit and return
    /// its sign, `(−1)^⌊e/2⌋` (mod 4).
    #[inline]
    pub(crate) fn set_phase_exponent(&mut self, e: u8) -> f64 {
        self.imag = e & 1 == 1;
        if e & 2 == 0 { 1.0 } else { -1.0 }
    }

    /// Build the singleton product `sin(x_id)`.
    pub(crate) fn sin(id: u32) -> Self {
        let mut p = Self::new();
        p.mul_sin(id);
        p
    }

    /// Build the singleton product `cos(x_id)`.
    pub(crate) fn cos(id: u32) -> Self {
        let mut p = Self::new();
        p.mul_cos(id);
        p
    }

    /// Total power of the sine factors (cached; `O(1)`).
    #[inline]
    pub(crate) fn sin_pow(&self) -> usize {
        self.sin_pow
    }

    /// Total power of the cosine factors (`O(factors)`; only `Display` needs it).
    #[inline]
    pub(crate) fn cos_pow(&self) -> usize {
        self.factors.iter().map(|f| f.cos as usize).sum()
    }

    /// Merge `f` into the sorted factor vector (exponents add).
    #[inline]
    pub(crate) fn merge_factor(&mut self, f: Factor) {
        match self.factors.binary_search_by_key(&f.var, |g| g.var) {
            Ok(i) => {
                self.factors[i].sin += f.sin;
                self.factors[i].cos += f.cos;
            }
            Err(i) => self.factors.insert(i, f),
        }
    }

    /// Debug-only check of canonicality and the cached sine degree.
    #[cfg(debug_assertions)]
    pub(crate) fn debug_check(&self) {
        for (i, f) in self.factors.iter().enumerate() {
            assert!(f.sin > 0 || f.cos > 0, "all-zero factor in Prod");
            if i > 0 {
                assert!(self.factors[i - 1].var < f.var, "Prod factors not sorted");
            }
        }
        let s: usize = self.factors.iter().map(|f| f.sin as usize).sum();
        assert_eq!(s, self.sin_pow, "sin_pow drifted");
    }

    #[cfg(not(debug_assertions))]
    #[inline]
    pub(crate) fn debug_check(&self) {}
}

/// `c₀ + Σᵢ cᵢ · pᵢ` over distinct monomials `pᵢ`.
///
/// The table is allocated on first insert. `aux` is a second, empty-between-uses
/// table that [`Sum::mul_term`] rebuilds into, so neither allocation is freed.
#[derive(Debug, Default)]
pub(crate) struct Sum {
    pub(crate) c0: f64,
    pub(crate) maps: Option<Box<SumMaps>>,
}

#[derive(Debug, Default)]
pub(crate) struct SumMaps {
    pub(crate) terms: FxHashMap<Prod, f64>,
    /// Persistent rebuild buffer for [`Sum::mul_term`], empty between operations.
    pub(crate) aux: FxHashMap<Prod, f64>,
}

impl Clone for Sum {
    /// Clone the value only. `aux` is workspace, empty between operations, so a
    /// clone starts it fresh rather than copying a buffer.
    #[inline]
    fn clone(&self) -> Self {
        Self {
            c0: self.c0,
            maps: self.maps.as_ref().map(|maps| {
                Box::new(SumMaps {
                    terms: maps.terms.clone(),
                    aux: FxHashMap::default(),
                })
            }),
        }
    }
}

/// Representational equality: exact `c0` and table contents. `aux` is not compared.
impl PartialEq for Sum {
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0
            && match (&self.maps, &other.maps) {
                (Some(lhs), Some(rhs)) => lhs.terms == rhs.terms,
                (Some(maps), None) | (None, Some(maps)) => maps.terms.is_empty(),
                (None, None) => true,
            }
    }
}

impl Sum {
    /// Construct an empty sum (value `0`).
    #[inline]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// The monomials, then a nonzero `c₀` as the identity monomial.
    pub(crate) fn parts(&self) -> impl Iterator<Item = (&Prod, f64)> + '_ {
        let c0 = (self.c0 != 0.0).then_some((&ONE, self.c0));
        self.iter().chain(c0)
    }

    /// The constant part `c₀`.
    #[cfg(test)]
    #[inline]
    pub(crate) fn c0(&self) -> f64 {
        self.c0
    }

    /// The number of monomials in the table (excluding `c₀`).
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.maps.as_ref().map_or(0, |maps| maps.terms.len())
    }

    /// Whether the monomial table is empty (`c₀` is not considered).
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.maps.as_ref().is_none_or(|maps| maps.terms.is_empty())
    }

    /// Iterate the monomial table in the backend's (deterministic, seed-free)
    /// order.
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&Prod, f64)> + '_ {
        self.maps
            .iter()
            .flat_map(|maps| maps.terms.iter())
            .map(|(p, c)| (p, *c))
    }

    /// Add `coeff · p`, dropping it if `p.sin_pow() > max` or `|coeff| <= min_eps`,
    /// and removing the monomial if the merged coefficient reaches `min_eps`.
    ///
    /// The identity monomial folds into `c₀`, which is never truncated. Dropping at
    /// insert is exact for `max` but not for `min_eps`: two sub-threshold
    /// contributions to one monomial are both dropped, where a post-pass would
    /// keep their sum.
    #[inline]
    pub(crate) fn add_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.is_identity() {
            self.c0 += coeff;
            return;
        }
        if p.sin_pow() > max || coeff.abs() <= min_eps {
            return;
        }
        let maps = self.maps.get_or_insert_with(Default::default);
        match maps.terms.entry(p) {
            Entry::Occupied(mut e) => {
                *e.get_mut() += coeff;
                if e.get().abs() <= min_eps {
                    e.remove();
                }
            }
            Entry::Vacant(e) => {
                e.insert(coeff);
            }
        }
    }

    /// Remove every monomial with sine degree above `max` or `|c| <= min_eps`.
    pub(crate) fn prune(&mut self, max: usize, min_eps: f64) {
        if let Some(maps) = &mut self.maps {
            maps.terms
                .retain(|p, c| p.sin_pow <= max && c.abs() > min_eps);
        }
    }
}

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

/// Merge two sorted factor vectors (exponents add), preserving canonicality.
pub(crate) fn merge_factors(lhs: &[Factor], rhs: &[Factor]) -> Vec<Factor> {
    let mut out = Vec::with_capacity(lhs.len() + rhs.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < lhs.len() && j < rhs.len() {
        let (a, b) = (lhs[i], rhs[j]);
        match a.var.cmp(&b.var) {
            Ordering::Less => {
                out.push(a);
                i += 1;
            }
            Ordering::Greater => {
                out.push(b);
                j += 1;
            }
            Ordering::Equal => {
                out.push(Factor::new(a.var, a.sin + b.sin, a.cos + b.cos));
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&lhs[i..]);
    out.extend_from_slice(&rhs[j..]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::{DEFAULT_MAX_SIN, DEFAULT_MIN_EPS};

    #[test]
    fn merge_factors_is_a_sorted_union_adding_exponents() {
        let f = Factor::new;
        let a = [f(0, 1, 0), f(2, 0, 1), f(5, 2, 0)];
        let b = [f(1, 1, 0), f(2, 1, 1), f(9, 0, 3)];
        let m = merge_factors(&a, &b);
        assert_eq!(
            m,
            vec![f(0, 1, 0), f(1, 1, 0), f(2, 1, 2), f(5, 2, 0), f(9, 0, 3)]
        );
        assert_eq!(merge_factors(&a, &[]), a.to_vec());
        assert_eq!(merge_factors(&[], &b), b.to_vec());
    }

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

    #[test]
    fn empty_and_unallocated_sums_compare_equal() {
        let mut s = Sum::new();
        s.add_term(Prod::sin(0), 1.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        s.maps.as_mut().unwrap().terms.clear();
        assert_eq!(s, Sum::new());
        assert_eq!(Sum::new(), s);
    }

    #[test]
    fn identity_monomial_folds_into_c0_but_i_does_not() {
        let mut s = Sum::new();
        s.add_term(Prod::new(), 2.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        assert_eq!((s.c0(), s.len()), (2.0, 0));
        let mut i = Prod::new();
        i.mul_phase(1);
        s.add_term(i, 3.0, DEFAULT_MAX_SIN, DEFAULT_MIN_EPS);
        assert_eq!((s.c0(), s.len()), (2.0, 1));
    }
}
