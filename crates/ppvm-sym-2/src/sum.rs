// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The monomial table [`Sum`].

use std::collections::hash_map::Entry;

use fxhash::FxHashMap;

use crate::prod::{ONE, Prod};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::term::{DEFAULT_MAX_SIN, DEFAULT_MIN_EPS};

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
