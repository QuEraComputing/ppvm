// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Pauli, PauliBits};

use super::data::{Decorated, OpPattern, PatternSite, PauliPattern, SiteSet};

impl PauliPattern {
    /// Build a pattern with `prefix` at absolute sites and a repeating tail.
    pub fn new(prefix: impl IntoIterator<Item = SiteSet>, tail: SiteSet) -> Self {
        let mut patterns: Vec<_> = prefix
            .into_iter()
            .enumerate()
            .map(|(position, set)| Decorated::Position(op_from_set(set), position))
            .collect();
        patterns.push(Decorated::Star(op_from_set(tail)));
        Self(patterns)
    }

    /// Build a pattern consisting of a single repeated atom.
    pub fn repeated(tail: SiteSet) -> Self {
        Self(vec![Decorated::Star(op_from_set(tail))])
    }

    /// The zero-state contraction pattern, `Z?*`.
    pub fn zero_state() -> Self {
        Self::repeated(SiteSet::Z.union(SiteSet::I))
    }

    /// Whether `word` is accepted, using the original greedy matcher.
    pub fn matches<W>(&self, word: &W) -> bool
    where
        W: PauliBits,
        W::Site: PatternSite,
    {
        let mut position = 0;
        for decorated in &self.0 {
            let matched = match decorated {
                Decorated::Position(op, target) => {
                    match_position(word, &mut position, *op, *target)
                }
                Decorated::Star(op) => match_star(word, &mut position, *op),
                Decorated::Repeat(op, count) => match_repeat(word, &mut position, *op, *count),
            };
            if !matched {
                return false;
            }
        }
        while position < word.n_sites() {
            if !is_identity(word, position) {
                return false;
            }
            position += 1;
        }
        true
    }

    /// Match a Pauli bit-vector, specializing the ubiquitous zero-state
    /// contraction to one X-plane scan.
    #[inline]
    pub(crate) fn matches_bits<W>(&self, word: &W) -> bool
    where
        W: PauliBits,
        W::Site: PatternSite,
    {
        if matches!(
            self.0.as_slice(),
            [Decorated::Star(OpPattern::SingleOrIdentity(Pauli::Z))]
        ) {
            return (0..word.n_sites()).all(|i| !word.x_bit(i) && !word.is_lost(i));
        }
        self.matches(word)
    }
}

fn match_position<W>(word: &W, position: &mut usize, pattern: OpPattern, target: usize) -> bool
where
    W: PauliBits,
{
    if *position > target {
        return false;
    }
    while *position < target {
        if *position >= word.n_sites() || !is_identity(word, *position) {
            return false;
        }
        *position += 1;
    }
    if *position >= word.n_sites() || !accepts_at(pattern, word, *position) {
        return false;
    }
    *position += 1;
    true
}

fn match_star<W>(word: &W, position: &mut usize, pattern: OpPattern) -> bool
where
    W: PauliBits,
{
    if pattern == OpPattern::AnyPauliOrIdentity {
        *position = word.n_sites();
        return true;
    }
    while *position < word.n_sites() && accepts_at(pattern, word, *position) {
        *position += 1;
    }
    true
}

fn match_repeat<W>(word: &W, position: &mut usize, pattern: OpPattern, count: usize) -> bool
where
    W: PauliBits,
{
    if count != 0 && pattern == OpPattern::AnyPauliOrIdentity {
        let Some(end) = position.checked_add(count) else {
            return false;
        };
        if end > word.n_sites() {
            return false;
        }
        *position = end;
        return true;
    }
    let mut matched = 0;
    while *position < word.n_sites() {
        if accepts_at(pattern, word, *position) {
            *position += 1;
            matched += 1;
            if matched == count {
                return true;
            }
        } else {
            return matched >= count;
        }
    }
    matched == count
}

#[inline(always)]
fn is_identity<W: PauliBits>(word: &W, position: usize) -> bool {
    !word.is_lost(position) && word.pauli_code(position) == 0
}

#[inline(always)]
fn accepts_at<W: PauliBits>(pattern: OpPattern, word: &W, position: usize) -> bool {
    if word.is_lost(position) {
        return pattern == OpPattern::AnyPauliOrIdentity;
    }
    let pauli = match word.pauli_code(position) {
        0 => Pauli::I,
        1 => Pauli::X,
        2 => Pauli::Z,
        3 => Pauli::Y,
        _ => unreachable!("a Pauli code has two bits"),
    };
    pattern.accepts(pauli)
}

fn op_from_set(set: SiteSet) -> OpPattern {
    let i = set.accepts(Pauli::I);
    let non_identity: Vec<_> = [Pauli::X, Pauli::Z, Pauli::Y]
        .into_iter()
        .filter(|&pauli| set.accepts(pauli))
        .collect();
    match (i, non_identity.as_slice()) {
        (true, []) => OpPattern::Identity,
        (false, [pauli]) => OpPattern::Single(*pauli),
        (false, [a, b]) => OpPattern::Double(*a, *b),
        (false, [_, _, _]) => OpPattern::AnyNonIdentity,
        (true, [pauli]) => OpPattern::SingleOrIdentity(*pauli),
        (true, [a, b]) => OpPattern::DoubleOrIdentity(*a, *b),
        (true, [_, _, _]) => OpPattern::AnyPauliOrIdentity,
        (false, []) => panic!("a Pauli pattern site set cannot be empty"),
        _ => unreachable!("there are only three non-identity Paulis"),
    }
}
