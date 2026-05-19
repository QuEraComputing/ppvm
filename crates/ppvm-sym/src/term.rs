// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{collections::BTreeMap, hash::Hash};

use fxhash::FxHashMap;

/// `<coeff> sin^m cos^n`
///
/// A single product of trigonometric atoms over symbolic variables.
/// Sines and cosines are grouped by variable id; powers are tracked as
/// totals so we can quickly compute and bound them. The order of
/// variables is kept canonical (ascending) so two `Prod` values
/// representing the same monomial compare equal.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Prod {
    pub(crate) sin: BTreeMap<u32, u32>,
    pub(crate) cos: BTreeMap<u32, u32>,
    pub(crate) sin_pow: usize,
    pub(crate) cos_pow: usize,
    /// phase factor mod 4, encoded as:
    /// |  | sign | imag |
    /// |--|------|------|
    /// |+1|    0 |    0 |
    /// |+i|    0 |    1 |
    /// |-1|    1 |    0 |
    /// |-i|    1 |    1 |
    pub(crate) phase: u8,
}

impl Default for Prod {
    fn default() -> Self {
        Self::new()
    }
}

impl Prod {
    /// Construct the empty product (the multiplicative identity, value `1`).
    pub fn new() -> Self {
        Self {
            sin: BTreeMap::new(),
            cos: BTreeMap::new(),
            sin_pow: 0,
            cos_pow: 0,
            phase: 0,
        }
    }

    /// Multiply the phase by `i^phase` (modulo `4`).
    pub fn add_phase(&mut self, phase: u8) {
        self.phase = (self.phase + phase) % 4;
    }

    /// Build the singleton product `sin(x_id)`.
    pub fn sin(id: u32) -> Self {
        let mut p = Self::new();
        p.sin.insert(id, 1);
        p.sin_pow = 1;
        p
    }

    /// Total power of all sine and cosine factors.
    pub fn pow(&self) -> usize {
        self.sin_pow + self.cos_pow
    }

    /// Total power of the sine factors.
    pub fn sin_pow(&self) -> usize {
        self.sin_pow
    }

    /// Total power of the cosine factors.
    pub fn cos_pow(&self) -> usize {
        self.cos_pow
    }

    /// Build the singleton product `cos(x_id)`.
    pub fn cos(id: u32) -> Self {
        let mut p = Self::new();
        p.cos.insert(id, 1);
        p.cos_pow = 1;
        p
    }
}

/// A formal sum `c₀ + Σᵢ cᵢ · pᵢ`, where each `pᵢ` is a [`Prod`] and
/// `cᵢ` is an `f64` coefficient.
#[derive(Debug, Clone, PartialEq)]
pub struct Sum {
    pub(crate) c0: f64,
    pub(crate) terms: FxHashMap<Prod, f64>,
}

impl Default for Sum {
    fn default() -> Self {
        Self::new()
    }
}

impl Sum {
    /// Construct an empty sum (value `0`).
    pub fn new() -> Self {
        Self {
            c0: 0.0,
            terms: FxHashMap::default(),
        }
    }

    /// Add the constant `c` into the sum's `c₀`, dropping it if
    /// `|c| < min_eps`.
    pub fn add_const(&mut self, c: f64, min_eps: f64) {
        if c.abs() < min_eps {
            return;
        }
        self.c0 += c;
    }

    /// Add `coeff · p` into the sum, subject to the same trunation
    /// constraints used elsewhere (`max` caps the sine power, `min_eps`
    /// drops near-zero coefficients).
    pub fn add_term(&mut self, p: Prod, coeff: f64, max: usize, min_eps: f64) {
        if p.sin_pow() > max || coeff.abs() < min_eps {
            return;
        }

        if p.pow() == 0 {
            self.c0 += coeff;
            return;
        }
        *self.terms.entry(p).or_insert(0.0) += coeff;
    }
}

/// Internal representation of a [`Term`].
///
/// `Sum` holds a full formal sum; `One` is the optimisation for a
/// single weighted product; `Var` is a bare symbolic variable used
/// before it has been wrapped in a `sin` or `cos`; `Const` is a
/// numeric scalar.
#[derive(Debug, Clone, PartialEq)]
pub enum Inner {
    /// A general sum of products.
    Sum(Sum),
    /// A single weighted product.
    One(Prod, f64),
    /// A bare symbolic variable (only valid as the argument of `sin`
    /// or `cos`).
    Var(u32),
    /// A numeric constant.
    Const(f64),
}

/// A symbolic polynomial in `sin(x_i)` and `cos(x_i)`.
///
/// `Term` is the public-facing wrapper around its [`Inner`] enum. It
/// also carries two truncation parameters, applied during multiplication
/// and addition:
///
/// # Examples
///
/// ```
/// use ppvm_sym::Term;
///
/// // sin²(x0) at x0 = π/2 equals 1.
/// let expr = Term::var(0).sin() * Term::var(0).sin();
/// let v = expr.eval(&[std::f64::consts::FRAC_PI_2]).unwrap();
/// assert!((v - 1.0).abs() < 1e-12);
/// ```
///
///
/// * `max_sin` — drop terms whose total sine power exceeds this bound.
/// * `min_eps` — drop terms whose coefficient magnitude falls below
///   this threshold.
#[derive(Debug, Clone, PartialEq)]
pub struct Term {
    pub(crate) inner: Inner,
    pub(crate) max_sin: usize, // max sin power
    pub(crate) min_eps: f64,   // min coefficient to keep
}

impl Term {
    /// Set the maximum sine power retained during arithmetic.
    pub fn set_max_sin(&mut self, max: usize) {
        self.max_sin = max;
    }

    /// Set the coefficient cutoff used during arithmetic.
    pub fn set_min_eps(&mut self, eps: f64) {
        self.min_eps = eps;
    }

    /// Construct a bare symbolic variable.
    pub fn var(u: u32) -> Self {
        Self {
            inner: Inner::Var(u),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }

    /// Apply `sin(·)` to the term. Only valid on variables and
    /// constants.
    pub fn sin(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::sin(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).sin();
            }
            _ => {
                panic!("only variable or constant can be input of sin");
            }
        }
        self
    }

    /// Apply `cos(·)` to the term. Only valid on variables and
    /// constants.
    pub fn cos(mut self) -> Self {
        match &mut self.inner {
            Inner::Var(u) => {
                self.inner = Inner::One(Prod::cos(*u), 1.0);
            }
            Inner::Const(c) => {
                *c = (*c).cos();
            }
            _ => {
                panic!("only variable or constant can be input of cos");
            }
        }
        self
    }

    /// Build a constant term from `c`.
    pub fn from_f64(c: f64) -> Self {
        Self {
            inner: Inner::Const(c),
            max_sin: usize::MAX,
            min_eps: f64::EPSILON,
        }
    }
}
