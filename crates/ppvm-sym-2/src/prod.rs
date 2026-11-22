// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The monomial [`Prod`] and its per-variable [`Factor`].

use std::cmp::Ordering;

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
}
