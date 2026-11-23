// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Readout by pattern: [`Sum::trace`] adds the coefficients of keys that a
//! per-site [`Pattern`] accepts. With `{I, Z}` on every site this is
//! `⟨0…0|O|0…0⟩`; conjugate by `X` first to read other basis states. `I`
//! denotes a live site, so lost sites only match patterns that accept `L`.

use std::ops::BitOr;

use num::Zero;

use ppvm_traits_2::{PauliBits, Support, loss::LossInfo};

use crate::Sum;

/// Accepted values for one site, from `I`, `X`, `Y`, `Z`, and lost (`L`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SiteSet(u8);

impl SiteSet {
    // Bit `1 << code` for packed codes `0=I, 1=X, 2=Z, 3=Y`; bit 4 is lost.
    pub const I: Self = Self(1);
    pub const X: Self = Self(1 << 1);
    pub const Z: Self = Self(1 << 2);
    pub const Y: Self = Self(1 << 3);
    pub const L: Self = Self(1 << 4);
    pub const ANY: Self = Self(0b1_1111);

    #[inline]
    fn accepts<W: PauliBits + LossInfo>(self, word: &W, site: usize) -> bool {
        let value = if word.is_lost(site) {
            Self::L.0
        } else {
            1 << word.pauli_code(site)
        };
        self.0 & value != 0
    }
}

impl BitOr for SiteSet {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// One accepted [`SiteSet`] per site; a word matches when every site does.
/// Build one from a `Vec<SiteSet>` in site order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pattern {
    sites: Vec<SiteSet>,
}

impl Pattern {
    /// `{I, Z}` on every site: readout in the all-zero state.
    pub fn zero_state(n_sites: usize) -> Self {
        vec![SiteSet::I | SiteSet::Z; n_sites].into()
    }

    fn matches<W: PauliBits + LossInfo>(&self, word: &W) -> bool {
        self.sites
            .iter()
            .enumerate()
            .all(|(site, set)| set.accepts(word, site))
    }
}

/// One set per site, in site order.
impl From<Vec<SiteSet>> for Pattern {
    fn from(sites: Vec<SiteSet>) -> Self {
        Self { sites }
    }
}

impl<S: Support, P> Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
{
    /// Sum of the coefficients whose keys `pattern` accepts, or zero if none do.
    /// Panics if the pattern's width differs from the sum's.
    pub fn trace(&self, pattern: &Pattern) -> S::Coeff {
        assert_eq!(
            pattern.sites.len(),
            self.n_sites(),
            "pattern width must match sum width"
        );
        let mut total = S::Coeff::zero();
        self.for_each_ref(|key, coeff| {
            if pattern.matches(key) {
                total += coeff;
            }
        });
        total
    }
}
