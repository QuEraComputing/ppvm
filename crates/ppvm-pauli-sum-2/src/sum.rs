// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::collections::HashSet;

use ppvm_traits_2::{
    Accumulate, Conjugate, IdentityBuildHasher, Indexable, Pair, Retain, Scale, Support, Word,
};

use crate::{NoPolicy, Policy, StoreAlloc};

/// Keys protected from truncation, using their structural digests for lookup.
pub type PreserveSet<K> = HashSet<K, IdentityBuildHasher>;

/// A sparse sum of words and coefficients, all with the same site count.
/// Storage determines the key and coefficient types. Reduction and truncation
/// are explicit; adding or scaling terms never removes zero coefficients.
/// Arithmetic retains the left operand's policy and preserved keys.
/// Addition and subtraction panic on unequal widths.
/// Exact and approximate comparisons require matching widths and stored keys,
/// including explicit zeros, and compare coefficients independently of policy.
/// Display sorts terms by weight with three coefficient decimals; Debug uses eight.
#[derive(Clone)]
pub struct Sum<S: Support, P = NoPolicy> {
    pub(crate) storage: S,
    policy: P,
    n_sites: usize,
    preserved_keys: PreserveSet<S::Key>,
}

impl<S: Support, P> Sum<S, P> {
    pub(crate) fn check_qubit(&self, qubit: usize) {
        assert!(qubit < self.n_sites(), "qubit out of bounds");
    }

    pub(crate) fn check_pair(&self, qubit0: usize, qubit1: usize) {
        self.check_qubit(qubit0);
        self.check_qubit(qubit1);
        assert_ne!(qubit0, qubit1, "two-qubit gates require distinct qubits");
    }

    /// Number of sites, including when the sum is empty.
    pub fn n_sites(&self) -> usize {
        self.n_sites
    }

    /// The sum's truncation policy.
    pub fn policy(&self) -> &P {
        &self.policy
    }

    /// Keys protected from [`Self::truncate`]. Absent keys are not inserted.
    /// Gates transform stored terms, leaving these key exemptions fixed.
    pub fn preserved_keys(&self) -> &PreserveSet<S::Key> {
        &self.preserved_keys
    }

    /// Remove all truncation exemptions without changing stored terms.
    pub fn clear_preserved_keys(&mut self) {
        self.preserved_keys.clear();
    }

    /// Number of stored terms, including zero coefficients.
    pub fn len(&self) -> usize {
        self.storage.len()
    }

    /// Whether there are no stored terms.
    pub fn is_empty(&self) -> bool {
        self.storage.is_empty()
    }

    /// Coefficient at `key`, or `None` if it is absent.
    pub fn get(&self, key: &S::Key) -> Option<S::Coeff> {
        self.storage.get(key)
    }

    /// Whether this exact key and coefficient are present.
    pub fn contains(&self, key: &S::Key, coeff: &S::Coeff) -> bool {
        self.get(key).is_some_and(|value| value == *coeff)
    }

    /// Whether a key is present, even with a zero coefficient.
    pub fn contains_key(&self, key: &S::Key) -> bool {
        self.get(key).is_some()
    }

    /// Iterate over owned key/coefficient pairs in storage order.
    pub fn iter(&self) -> impl Iterator<Item = (S::Key, S::Coeff)> + '_ {
        self.storage.iter()
    }

    /// Visit terms by reference, avoiding copies when storage supports it.
    pub fn for_each_ref(&self, f: impl FnMut(&S::Key, &S::Coeff)) {
        self.storage.for_each_ref(f);
    }
}

impl<S: Support, P> Sum<S, P>
where
    S::Key: Word + Indexable,
{
    /// Protect a key from truncation. Panics if its width does not match.
    pub fn preserve_key(&mut self, key: S::Key) {
        assert_eq!(
            key.n_sites(),
            self.n_sites,
            "preserved key width must match sum width"
        );
        self.preserved_keys.insert(key);
    }

    /// Protect these keys from truncation; all must match the sum's width.
    #[must_use]
    pub fn preserving(mut self, keys: impl IntoIterator<Item = S::Key>) -> Self {
        for key in keys {
            self.preserve_key(key);
        }
        self
    }
}

impl<S: Support + StoreAlloc, P> Sum<S, P> {
    /// Construct an empty sum with an explicit policy and capacity hint.
    pub fn with_capacity(n_sites: usize, policy: P, capacity: usize) -> Self {
        Self {
            storage: S::with_capacity(capacity),
            policy,
            n_sites,
            preserved_keys: PreserveSet::default(),
        }
    }

    /// Remove all terms, retaining width, policy, preserved keys, and capacity.
    pub fn clear(&mut self) {
        self.storage.clear();
    }
}

impl<S, P> Sum<S, P>
where
    S: Support + StoreAlloc,
    P: Policy<S::Key, S::Coeff>,
{
    /// Construct an empty sum with the default policy.
    pub fn new(n_sites: usize) -> Self
    where
        P: Default,
    {
        Self::with_policy(n_sites, P::default())
    }

    /// Construct an empty sum using the policy's initial capacity hint.
    pub fn with_policy(n_sites: usize, policy: P) -> Self {
        let capacity = policy.capacity(n_sites);
        Self::with_capacity(n_sites, policy, capacity)
    }
}

impl<S, P> Sum<S, P>
where
    S: Accumulate + StoreAlloc,
    S::Key: Word,
    P: Policy<S::Key, S::Coeff>,
{
    /// Accumulate terms with the default policy. Panics on a width mismatch.
    pub fn from_terms(n_sites: usize, terms: impl IntoIterator<Item = (S::Key, S::Coeff)>) -> Self
    where
        P: Default,
    {
        Self::from_terms_with_policy(n_sites, P::default(), terms)
    }

    /// Accumulate terms with an explicit policy. Panics on a width mismatch.
    pub fn from_terms_with_policy(
        n_sites: usize,
        policy: P,
        terms: impl IntoIterator<Item = (S::Key, S::Coeff)>,
    ) -> Self {
        let mut sum = Self::with_policy(n_sites, policy);
        for (key, coeff) in terms {
            sum.add_term(key, coeff);
        }
        sum
    }
}

impl<S: Accumulate, P> Sum<S, P> {
    /// Add to the coefficient at `key`, inserting it if absent.
    /// Panics before changing the sum if the word's width does not match.
    pub fn add_term(&mut self, key: S::Key, coeff: S::Coeff)
    where
        S::Key: Word,
    {
        assert_eq!(
            key.n_sites(),
            self.n_sites,
            "term width must match sum width"
        );
        self.storage.accumulate(key, coeff);
    }

    /// Remove coefficients that are exactly zero, including preserved keys.
    pub fn reduce(&mut self) {
        self.storage.reduce();
    }
}

impl<S: Scale, P> Sum<S, P> {
    /// Multiply every coefficient by `scalar` without removing any terms.
    pub fn scale(&mut self, scalar: &S::Coeff) {
        self.storage.scale(scalar);
    }
}

impl<S: Pair, P> Sum<S, P> {
    /// Sum products of matching coefficients, without conjugation.
    /// Policies may differ. Panics if the sums have different widths.
    pub fn overlap<P2>(&self, other: &Sum<S, P2>) -> S::Coeff {
        assert_eq!(self.n_sites, other.n_sites, "overlap requires equal widths");
        self.storage.overlap(&other.storage)
    }

    /// Sum products of matching coefficients, conjugating this sum's coefficients.
    /// Policies may differ. Panics if the sums have different widths.
    pub fn hermitian_overlap<P2>(&self, other: &Sum<S, P2>) -> S::Coeff
    where
        S::Coeff: Conjugate,
    {
        assert_eq!(self.n_sites, other.n_sites, "overlap requires equal widths");
        self.storage.hermitian_overlap(&other.storage)
    }
}

impl<S, P> Sum<S, P>
where
    S: Support + Retain<S::Key, S::Coeff>,
    S::Key: Indexable,
    P: Policy<S::Key, S::Coeff>,
{
    /// Apply the policy, leaving preserved terms and their coefficients intact.
    /// [`NoPolicy`] leaves every term in place; absent preserved keys stay absent.
    pub fn truncate(&mut self) {
        self.policy.truncate(&mut ProtectedStorage {
            storage: &mut self.storage,
            preserved_keys: &self.preserved_keys,
        });
    }
}

/// Let policies filter storage without removing preserved keys.
struct ProtectedStorage<'a, S, K> {
    storage: &'a mut S,
    preserved_keys: &'a PreserveSet<K>,
}

impl<S, K, C> Retain<K, C> for ProtectedStorage<'_, S, K>
where
    S: Retain<K, C>,
    K: Indexable,
{
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        if self.preserved_keys.is_empty() {
            self.storage.retain(keep);
        } else {
            self.storage
                .retain(|key, coeff| keep(key, coeff) || self.preserved_keys.contains(key));
        }
    }
}
