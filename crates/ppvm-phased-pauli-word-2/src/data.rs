// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The generic [`Phased`] wrapper, its inherent constructors/accessors, the
//! delegating [`Word`] inspection impl, the signed-string parsers, and the
//! structural [`PartialEq`]/[`Eq`]/[`Clone`]/[`Display`] that agree on the
//! logical identity `(phase, word)`.
//!
//! `Phased<W>` pairs a base word `W` with an explicit `ℤ₄` phase (a
//! [`Phase`]), giving the standalone phased operator `i^φ · g(w)`. The phase does
//! **not** affect site inspection, so [`Word`] delegates verbatim to the inner
//! word; the phase *does* participate in equality, `Display`, and (later) the
//! Clifford sign and the twisted product.
//!
//! Design: `word-data-structures.md` §"Phased words" (the `Phased<W>` struct and
//! the `PhasedPauliWord` alias) and §"Ordering and serialization" (phase
//! participates in equality/serialization but never in map-key hashing, because
//! `Phased<W>` is not [`Indexable`](ppvm_traits_2)). Ported from
//! `ppvm-pauli-word/src/phase/data.rs`. Lean spec: the phased Pauli group `𝒫₁`
//! (`ℤ₄ × 𝔽₂²`) of `lean/PPVM/Pauli/Phase.lean`.

use std::fmt;
use std::hash::BuildHasher;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage, PauliWord};
use ppvm_traits_2::{Phase, Word};

/// A base word `W` paired with an explicit `ℤ₄` phase — the standalone phased
/// operator `i^φ · g(w)`.
///
/// The `phase` field is a [`Phase`] (`{+1, +i, −1, −i} ≅ ℤ/4`); the `word` field
/// is the bare base word carrying no phase of its own. Unlike the bare word,
/// whose `PhaseTrack` drops every conjugation sign, `Phased<W>` retains it: this
/// is the type that answers "what is the sign of `HYH`?".
///
/// The structural identity is `(phase, word)`; both participate in equality,
/// `Clone`, and `Display` (`word-data-structures.md` §"Ordering and
/// serialization"). `Phased<W>` deliberately implements **neither** `Hash` nor
/// [`Indexable`](ppvm_traits_2): the phase is part of its identity, so it is not
/// a production map key (`word-data-structures.md` §"Phased words").
///
/// # Examples
///
/// ```
/// use ppvm_phased_pauli_word_2::PhasedPauliWord;
/// use ppvm_traits_2::{Phase, Word};
///
/// let pw: PhasedPauliWord = "+iXYZI".into();
/// assert_eq!(pw.n_sites(), 4);
/// assert_eq!(pw.phase(), Phase::PosI); // +i
/// assert!(pw.is_positive());
///
/// let neg: PhasedPauliWord = "-XYZI".into();
/// assert_eq!(neg.phase(), Phase::Neg1); // -1
/// assert!(!neg.is_positive());
/// ```
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

    /// `true` when the real part of the prefactor is non-negative (`+1` or `+i`).
    #[inline]
    pub fn is_positive(&self) -> bool {
        matches!(self.phase, Phase::Pos1 | Phase::PosI)
    }
}

/// Site inspection delegates entirely to the inner word: the phase is a scalar
/// prefactor and does not change any site (`word-data-structures.md` §"Phased
/// words").
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

/// Render the `ℤ₄` phase as the `+ / +i / − / −i` prefix, matching the old
/// crate's `PhasedPauliWord` display (`ppvm-pauli-word/src/phase/data.rs`).
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
    /// Parse a `"[sign][i]<PauliString>"` literal into a phased Pauli word:
    /// `+`/`-` set the sign, an optional `i` immediately after sets the
    /// imaginary flag, and the remainder is the [`PauliWord`] string. Ported from
    /// `ppvm-pauli-word/src/phase/data.rs`. Panics on a missing/invalid sign.
    fn from(s: &str) -> Self {
        let mut chars = s.chars();
        let (phase, prefix_len) = match (chars.next(), chars.next()) {
            (Some('+'), Some('i')) => (Phase::PosI, 2),
            (Some('-'), Some('i')) => (Phase::NegI, 2),
            (Some('+'), _) => (Phase::Pos1, 1),
            (Some('-'), _) => (Phase::Neg1, 1),
            _ => panic!("invalid phase format: {s:?} (expected a leading +/-)"),
        };
        let body: String = s.chars().skip(prefix_len).collect();
        Self {
            word: PauliWord::<A, H>::from(body.as_str()),
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
    use crate::PhasedPauliWord;
    use ppvm_traits_2::Pauli;

    type Pw = Phased<PauliWord>;

    #[test]
    fn parse_sign_and_body() {
        let cases = [
            ("+XYZI", Phase::Pos1),
            ("-XYZI", Phase::Neg1),
            ("+iXYZI", Phase::PosI),
            ("-iXYZI", Phase::NegI),
        ];
        for (s, want) in cases {
            let pw: Pw = s.into();
            assert_eq!(pw.phase(), want, "{s} phase");
            assert_eq!(pw.n_sites(), 4, "{s} width");
            assert_eq!(pw.get(0), Pauli::X);
            assert_eq!(pw.get(1), Pauli::Y);
            assert_eq!(pw.get(2), Pauli::Z);
            assert_eq!(pw.get(3), Pauli::I);
        }
    }

    #[test]
    fn display_roundtrips_prefix() {
        for s in ["+XYZI", "-XYZI", "+iXYZI", "-iXYZI"] {
            let pw: Pw = s.into();
            assert_eq!(pw.to_string(), s);
        }
    }

    #[test]
    fn is_positive_reads_sign() {
        let pos: Pw = "+XY".into();
        let posi: Pw = "+iXY".into();
        let neg: Pw = "-XY".into();
        let negi: Pw = "-iXY".into();
        assert!(pos.is_positive() && posi.is_positive());
        assert!(!neg.is_positive() && !negi.is_positive());
    }

    #[test]
    fn phase_participates_in_equality() {
        let a: Pw = "+XY".into();
        let b: Pw = "-XY".into();
        let c: Pw = "+XY".into();
        assert_ne!(a, b, "different phase ⇒ distinct");
        assert_eq!(a, c);
    }

    #[test]
    fn weight_iter_delegate_to_word() {
        let pw: Pw = "-iXIYZ".into();
        assert_eq!(pw.weight(), 3);
        let via_iter: Vec<Pauli> = pw.iter().collect();
        assert_eq!(via_iter, vec![Pauli::X, Pauli::I, Pauli::Y, Pauli::Z]);
    }

    #[test]
    fn send_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PhasedPauliWord>();
    }
}
