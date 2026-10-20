// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Phased words with explicit prefactors, parsing, display, and equality.

use std::fmt;
use std::hash::BuildHasher;

use crate::{HashFinalize, PauliStorage, PauliWord};
use ppvm_traits_2::{Phase, Word};

/// A base word with a phase in `{+1, +i, -1, -i}`, representing `i^φ · g(w)`.
/// Equality and display include the phase; site inspection uses only the word.
/// Implements neither `Hash` nor `Indexable`, so it is not a map key.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Phased<W> {
    /// The bare base word (carries no phase of its own).
    pub(crate) word: W,
    /// The explicit `ℤ₄` prefactor `i^φ`.
    pub(crate) phase: Phase,
}

impl<W> Phased<W> {
    /// Wrap `word` with the trivial phase `+1`.
    #[inline]
    pub fn new(word: W) -> Self {
        Self {
            word,
            phase: Phase::one(),
        }
    }

    /// Wrap `word` with an explicit `phase`.
    #[inline]
    pub fn with_phase(word: W, phase: Phase) -> Self {
        Self { word, phase }
    }

    /// The `ℤ₄` prefactor `i^φ`.
    #[inline]
    pub fn phase(&self) -> Phase {
        self.phase
    }

    /// Borrow the bare base word (phase stripped).
    #[inline]
    pub fn word(&self) -> &W {
        &self.word
    }

    /// Consume the wrapper, returning the base word and its phase.
    #[inline]
    pub fn into_parts(self) -> (W, Phase) {
        (self.word, self.phase)
    }

    /// Multiply the stored phase by `delta` (the `ℤ₄` group product).
    #[inline]
    pub fn add_phase(&mut self, delta: Phase) {
        self.phase *= delta;
    }

    /// Whether the phase has a positive sign (`+1` or `+i`).
    #[inline]
    pub fn is_positive(&self) -> bool {
        matches!(self.phase, Phase::Pos1 | Phase::PosI)
    }
}

/// Site inspection delegates to the inner word; the scalar phase changes no sites.
impl<W: Word> Word for Phased<W> {
    type Site = W::Site;

    #[inline]
    fn n_sites(&self) -> usize {
        self.word.n_sites()
    }

    #[inline]
    fn get(&self, index: usize) -> Self::Site {
        self.word.get(index)
    }

    #[inline]
    fn weight(&self) -> usize {
        self.word.weight()
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = Self::Site> {
        self.word.iter()
    }
}

/// Renders the phase as `+`, `+i`, `-`, or `-i`, matching the legacy display.
fn phase_prefix(phase: Phase) -> &'static str {
    match phase {
        Phase::Pos1 => "+",
        Phase::PosI => "+i",
        Phase::Neg1 => "-",
        Phase::NegI => "-i",
    }
}

/// `Display` prints the phase prefix followed by the base word, e.g. `+iXYZI`.
impl<W: fmt::Display> fmt::Display for Phased<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", phase_prefix(self.phase), self.word)
    }
}

impl<A, H> From<&str> for Phased<PauliWord<A, H>>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    /// Parses a required `+`/`-`, optional `i`, then a [`PauliWord`] string.
    /// Panics if the sign is missing or the word is invalid.
    fn from(s: &str) -> Self {
        let mut chars = s.chars();
        let (phase, prefix_len) = match (chars.next(), chars.next()) {
            (Some('+'), Some('i')) => (Phase::PosI, 2),
            (Some('-'), Some('i')) => (Phase::NegI, 2),
            (Some('+'), _) => (Phase::Pos1, 1),
            (Some('-'), _) => (Phase::Neg1, 1),
            _ => panic!("invalid phase format: {s:?} (expected a leading +/-)"),
        };
        // Recognized prefixes are ASCII, so `prefix_len` is a UTF-8 boundary.
        Self {
            word: PauliWord::<A, H>::from(&s[prefix_len..]),
            phase,
        }
    }
}

impl<A, H> From<String> for Phased<PauliWord<A, H>>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn from(s: String) -> Self {
        Self::from(s.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_prefixes_roundtrip() {
        for (prefix, phase) in [
            ("+", Phase::Pos1),
            ("+i", Phase::PosI),
            ("-", Phase::Neg1),
            ("-i", Phase::NegI),
        ] {
            for body in ["", "XYZI"] {
                let text = format!("{prefix}{body}");
                let word = Phased::<PauliWord>::from(text.as_str());
                assert_eq!(word.phase(), phase);
                assert_eq!(word.word().to_string(), body);
                assert_eq!(word.to_string(), text);
            }
        }
    }

    #[test]
    #[should_panic(expected = "invalid Pauli character")]
    fn unicode_body_is_rejected_by_word_parser() {
        let _ = Phased::<PauliWord>::from("+λ");
    }
}
