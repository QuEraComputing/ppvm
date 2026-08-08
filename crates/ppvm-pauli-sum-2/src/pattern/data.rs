// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{LossySite, Pauli};

/// The set of Paulis accepted at one site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SiteSet(pub(crate) u8);

impl SiteSet {
    pub const I: Self = Self(1);
    pub const X: Self = Self(1 << 1);
    pub const Z: Self = Self(1 << 2);
    pub const Y: Self = Self(1 << 3);
    pub const ANY: Self = Self(0b1111);
    pub const NON_IDENTITY: Self = Self(0b1110);

    #[inline]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[inline]
    pub const fn accepts(self, pauli: Pauli) -> bool {
        self.0 & pauli_bit(pauli) != 0
    }

    /// Accepted Paulis in the `-2` tower's conventional `I, X, Y, Z` order.
    pub fn iter(self) -> impl Iterator<Item = Pauli> {
        [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z]
            .into_iter()
            .filter(move |&pauli| self.accepts(pauli))
    }
}

const fn pauli_bit(pauli: Pauli) -> u8 {
    match pauli {
        Pauli::I => 1,
        Pauli::X => 1 << 1,
        Pauli::Z => 1 << 2,
        Pauli::Y => 1 << 3,
    }
}

/// Convert a word site to a Pauli.
///
/// Lost sites are represented by `None`. They fail ordinary atoms but retain
/// old's special behavior where `_` / `[XYZ]?` accepts every alphabet symbol.
pub trait PatternSite {
    fn to_pauli(self) -> Option<Pauli>;
}

impl PatternSite for Pauli {
    #[inline]
    fn to_pauli(self) -> Option<Pauli> {
        Some(self)
    }
}

impl PatternSite for LossySite<Pauli> {
    #[inline]
    fn to_pauli(self) -> Option<Pauli> {
        match self {
            LossySite::Present(pauli) => Some(pauli),
            LossySite::Lost => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OpPattern {
    Identity,
    Single(Pauli),
    Double(Pauli, Pauli),
    AnyNonIdentity,
    SingleOrIdentity(Pauli),
    DoubleOrIdentity(Pauli, Pauli),
    AnyPauliOrIdentity,
}

impl OpPattern {
    pub(crate) fn optional(self) -> Self {
        match self {
            Self::Single(pauli) => Self::SingleOrIdentity(pauli),
            Self::Double(left, right) => Self::DoubleOrIdentity(left, right),
            Self::AnyNonIdentity => Self::AnyPauliOrIdentity,
            other => other,
        }
    }

    pub(crate) fn accepts(self, pauli: Pauli) -> bool {
        self.site_set().accepts(pauli)
    }

    pub(crate) fn accepts_site(self, pauli: Option<Pauli>) -> bool {
        match pauli {
            Some(pauli) => self.accepts(pauli),
            // Old's `_` / `[XYZ]?` arm returns `true` without inspecting the
            // symbol, so it also accepts the lossy alphabet's `L`.
            None => self == Self::AnyPauliOrIdentity,
        }
    }

    pub(crate) fn site_set(self) -> SiteSet {
        match self {
            Self::Identity => SiteSet::I,
            Self::Single(pauli) => SiteSet(pauli_bit(pauli)),
            Self::Double(a, b) => SiteSet(pauli_bit(a) | pauli_bit(b)),
            Self::AnyNonIdentity => SiteSet::NON_IDENTITY,
            Self::SingleOrIdentity(pauli) => SiteSet(pauli_bit(pauli) | 1),
            Self::DoubleOrIdentity(a, b) => SiteSet(pauli_bit(a) | pauli_bit(b) | 1),
            Self::AnyPauliOrIdentity => SiteSet::ANY,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Decorated {
    Position(OpPattern, usize),
    Star(OpPattern),
    Repeat(OpPattern, usize),
}

/// A stateful Pauli-word pattern compatible with the original grammar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PauliPattern(pub(crate) Vec<Decorated>);
