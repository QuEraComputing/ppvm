// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::iter::Peekable;

use ppvm_traits_2::{Pauli, Word};

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
        W: Word,
        W::Site: PatternSite,
    {
        let mut sites = word
            .iter()
            .enumerate()
            .map(|(position, site)| (position, site.to_pauli()))
            .peekable();
        for decorated in &self.0 {
            let matched = match decorated {
                Decorated::Position(op, position) => match_position(&mut sites, *op, *position),
                Decorated::Star(op) => match_star(&mut sites, *op),
                Decorated::Repeat(op, count) => match_repeat(&mut sites, *op, *count),
            };
            if !matched {
                return false;
            }
        }
        sites.all(|(_, pauli)| pauli == Some(Pauli::I))
    }
}

fn match_position<I>(sites: &mut Peekable<I>, pattern: OpPattern, position: usize) -> bool
where
    I: Iterator<Item = (usize, Option<Pauli>)>,
{
    for (site_position, pauli) in sites.by_ref() {
        if site_position == position {
            return pattern.accepts_site(pauli);
        }
        if pauli != Some(Pauli::I) {
            return false;
        }
    }
    false
}

fn match_star<I>(sites: &mut Peekable<I>, pattern: OpPattern) -> bool
where
    I: Iterator<Item = (usize, Option<Pauli>)>,
{
    while sites
        .peek()
        .is_some_and(|(_, pauli)| pattern.accepts_site(*pauli))
    {
        sites.next();
    }
    true
}

fn match_repeat<I>(sites: &mut Peekable<I>, pattern: OpPattern, count: usize) -> bool
where
    I: Iterator<Item = (usize, Option<Pauli>)>,
{
    let mut matched = 0;
    while let Some((_, pauli)) = sites.peek() {
        if pattern.accepts_site(*pauli) {
            sites.next();
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
