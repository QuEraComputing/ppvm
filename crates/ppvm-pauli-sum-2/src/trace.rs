// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Expectation extraction: the [`Trace`] pairing of a [`Sum`] against a
//! **pattern** — the set of Pauli words a measurement contracts onto.
//!
//! This is the port of the old crate's `Trace` surface
//! (`ppvm-pauli-word/src/pattern/`, `ppvm-pauli-sum/src/sum/trace.rs`), which is
//! how every expectation value in the old test-suite is read out: propagate an
//! observable backwards through the circuit, then contract the resulting sum
//! against the state's stabilizer pattern. For the all-zero state that pattern is
//! `Z?*` — every site is `I` or `Z` — because `⟨0…0|P|0…0⟩` is `1` for such a `P`
//! and `0` for any `P` carrying an `X` or `Y`. So
//!
//! ```text
//! ⟨0…0| O |0…0⟩ = Σ_{k ∈ supp(O), k matches Z?*} O[k]
//! ```
//!
//! which is exactly [`Trace::trace`] here. It backs `tests/ghz.rs`,
//! `tests/cnot.rs` and `tests/noise.rs::test_depolarizing_error` in the old
//! crate (e.g. three independent `depolarize1(i, pᵢ)` on `ZZZ` contract to
//! `Πᵢ (1 − 4pᵢ/3)`).
//!
//! Design: §"Compatibility with current names" retains
//! [`Trace<'a, RHS>`](ppvm_traits_2::Trace) in `ppvm-traits-2` precisely so its
//! implementers can land with the pattern port; this is that landing. The
//! same-type sum-against-sum pairing is [`Sum::overlap`], whose bilinearity is
//! machine-checked in `lean/PPVM/Algebra/GradedMap.lean`.
//!
//! # Scope
//!
//! Old's `PauliPattern` carries a whole textual grammar (alternation, counted
//! repetition, absolute-position anchoring). Ported here is the fragment the
//! expectation path actually consumes — a per-site alphabet mask, optionally with
//! a repeating tail — which covers `Z?*` and every literal/wildcard word. The
//! parser and the anchored/counted decorations are deliberately left out until a
//! caller needs them; nothing in the sum engine does.

use ppvm_traits_2::{Accumulate, Indexable, Pauli, Trace, Word};

use crate::policy::Policy;
use crate::sum::Sum;

/// The set of single-qubit Paulis a pattern accepts at one site — a 4-bit mask
/// over `I`/`X`/`Y`/`Z`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SiteSet(u8);

impl SiteSet {
    /// Accepts only the identity (`I`).
    pub const I: Self = Self(1 << 0);
    /// Accepts only `X`.
    pub const X: Self = Self(1 << 1);
    /// Accepts only `Y`.
    pub const Y: Self = Self(1 << 2);
    /// Accepts only `Z`.
    pub const Z: Self = Self(1 << 3);
    /// Accepts any Pauli, identity included — old's `_` / `[XYZ]?`.
    pub const ANY: Self = Self(0b1111);
    /// Accepts any non-identity Pauli — old's `[XYZ]`.
    pub const NON_IDENTITY: Self = Self(0b1110);

    /// The union of two site sets — old's `?` suffix is `s.union(SiteSet::I)`.
    #[inline]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    /// Whether this set accepts `p`.
    #[inline]
    pub const fn accepts(self, p: Pauli) -> bool {
        let bit = match p {
            Pauli::I => 1 << 0,
            Pauli::X => 1 << 1,
            Pauli::Y => 1 << 2,
            Pauli::Z => 1 << 3,
        };
        self.0 & bit != 0
    }
}

/// A set of Pauli words, as a per-site alphabet mask with a repeating tail.
///
/// Sites `0..prefix.len()` must match the corresponding `prefix` entry; every
/// site beyond that must match `tail`. A pattern with an empty prefix is old's
/// bare `<pat>*` form.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PauliPattern {
    prefix: Vec<SiteSet>,
    tail: SiteSet,
}

impl PauliPattern {
    /// A pattern whose every site must match `tail` — old's `<pat>*`.
    #[inline]
    pub fn repeated(tail: SiteSet) -> Self {
        Self {
            prefix: Vec::new(),
            tail,
        }
    }

    /// A pattern anchoring the leading sites to `prefix`, with `tail` for the
    /// rest.
    #[inline]
    pub fn new(prefix: impl IntoIterator<Item = SiteSet>, tail: SiteSet) -> Self {
        Self {
            prefix: prefix.into_iter().collect(),
            tail,
        }
    }

    /// `Z?*` — every site is `I` or `Z`.
    ///
    /// The all-zero computational-basis state's stabilizer pattern: a Pauli word
    /// has non-zero `⟨0…0|·|0…0⟩` exactly when it is diagonal, and then that
    /// expectation is `+1`. Contracting a propagated observable against this
    /// pattern is the zero-state expectation value.
    #[inline]
    pub fn zero_state() -> Self {
        Self::repeated(SiteSet::Z.union(SiteSet::I))
    }

    /// Whether `word` is in the set this pattern denotes.
    pub fn matches<W>(&self, word: &W) -> bool
    where
        W: Word<Site = Pauli>,
    {
        word.iter().enumerate().all(|(i, p)| {
            let set = self.prefix.get(i).copied().unwrap_or(self.tail);
            set.accepts(p)
        })
    }
}

/// A pattern traces a *word* to a boolean membership — old's
/// `impl Trace<'a, W> for PauliPattern` (`Output = bool`), which is the shape the
/// sum-side trace below consumes.
impl<'a, W> Trace<'a, W> for PauliPattern
where
    W: Word<Site = Pauli> + 'a,
{
    type Output = bool;

    #[inline]
    fn trace(&'a self, value: &'a W) -> bool {
        self.matches(value)
    }
}

/// A sum traces against a pattern to the **sum of the coefficients of the
/// matching keys**.
///
/// Old's `impl Trace<'a, Rhs> for PauliSum`, which delegates to the map's
/// `fold(zero, |acc, (k, v)| { value.trace(k).then(|| acc += v); acc })`
/// (`ppvm-traits/src/map/hashmap.rs`). Zero-coefficient terms contribute zero, so
/// the contract is insensitive to whether the support carries them; the
/// *summation order* is the backend's iteration order on both sides, so a
/// floating-point contraction agrees only up to reassociation error.
///
/// This is stated for the concrete [`PauliPattern`] rather than for every
/// `R: Trace<'a, S::Key, Output = bool>`: the engine's read side
/// ([`Support::iter`](ppvm_traits_2::Support::iter)) yields **owned** keys, which
/// cannot satisfy the `&'a RHS` the generic form would demand. Old could be
/// generic because its map iterator borrows. `PauliPattern` is the only
/// implementer old ever had.
impl<'a, S, P> Trace<'a, PauliPattern> for Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word<Site = Pauli> + Indexable,
{
    type Output = S::Coeff;

    fn trace(&'a self, value: &'a PauliPattern) -> S::Coeff {
        self.iter()
            .filter(|(k, _)| value.matches(k))
            .map(|(_, c)| c)
            .sum()
    }
}
