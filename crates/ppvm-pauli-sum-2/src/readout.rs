// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Readout by pattern: [`Sum::trace`] adds the coefficients of keys that a
//! per-site [`Pattern`] accepts. With `{I, Z}` on every site this is
//! `⟨0…0|O|0…0⟩`; conjugate by `X` first to read other basis states. `I`
//! denotes a live site, so lost sites only match patterns that accept `L`.

use std::{fmt, ops::BitOr, str::FromStr};

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
/// Parses from one token per site: `I`, `X`, `Y`, `Z`, `L`, `_` (any), or a
/// bracketed union such as `[IZ]`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Pattern {
    sites: Vec<SiteSet>,
}

impl Pattern {
    /// The same set on every site.
    pub fn uniform(n_sites: usize, set: SiteSet) -> Self {
        Self {
            sites: vec![set; n_sites],
        }
    }

    /// `{I, Z}` on every site: readout in the all-zero state.
    pub fn zero_state(n_sites: usize) -> Self {
        Self::uniform(n_sites, SiteSet::I | SiteSet::Z)
    }

    pub fn n_sites(&self) -> usize {
        self.sites.len()
    }

    /// Whether every site of `word` is accepted. Panics if the widths differ.
    pub fn matches<W: PauliBits + LossInfo>(&self, word: &W) -> bool {
        assert_eq!(
            word.n_sites(),
            self.n_sites(),
            "pattern width must match word width"
        );
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

/// An invalid pattern string, with the reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsePatternError(String);

impl fmt::Display for ParsePatternError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid pattern: {}", self.0)
    }
}

impl std::error::Error for ParsePatternError {}

fn letter(c: char) -> Option<SiteSet> {
    Some(match c {
        'I' => SiteSet::I,
        'X' => SiteSet::X,
        'Y' => SiteSet::Y,
        'Z' => SiteSet::Z,
        'L' => SiteSet::L,
        _ => return None,
    })
}

impl FromStr for Pattern {
    type Err = ParsePatternError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let error = |reason: String| Err(ParsePatternError(reason));
        let mut sites = Vec::new();
        let mut chars = s.chars();
        while let Some(c) = chars.next() {
            let set = match c {
                '_' => SiteSet::ANY,
                '[' => {
                    let mut set = SiteSet(0);
                    loop {
                        match chars.next() {
                            Some(']') if set != SiteSet(0) => break set,
                            Some(']') => return error("empty `[]`".into()),
                            Some(c) => match letter(c) {
                                Some(value) => set = set | value,
                                None => return error(format!("unexpected `{c}` in `[...]`")),
                            },
                            None => return error("unclosed `[`".into()),
                        }
                    }
                }
                c => match letter(c) {
                    Some(set) => set,
                    None => return error(format!("unexpected `{c}`")),
                },
            };
            sites.push(set);
        }
        Ok(Self { sites })
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
            pattern.n_sites(),
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
