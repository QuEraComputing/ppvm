// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`Sum<S, P>`] — the graded sparse-sum engine: a finitely-supported linear
//! combination `Σ cₖ·k` over the graded algebra `C[K]`, backed by an
//! [`Accumulate`] container `S` and truncated by a [`Policy`] `P`.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`" (the `Sum<S, P>` struct: storage + policy + `n_sites`, no owned
//! workspace) and §"apply". The module laws the engine relies on are
//! machine-checked in `lean/PPVM/Algebra/GradedMap.lean`.

use ppvm_traits_2::{
    Accumulate, Conjugate, Indexable, Pair, Retain, Scale, TermBatch, TermProducer, TermSink, Word,
};

use crate::policy::Policy;
use crate::store::{RekeyBijective, RotateInPlace, ScaleByKey, SignFlipByKey, StoreAlloc};

/// A sparse formal sum `Σ cₖ·k` over the free `C`-module on the keys `K`.
///
/// `Sum` owns **only** its storage, policy, and width — no auxiliary map, no
/// scratch buffer (Design: §"There is **no `SumStorage` trait, and no owned
/// workspace**"). `Clone` is therefore pure data, which matters because a
/// stabilizer mixture clones frequently.
///
/// The key and coefficient are `S::Key` / `S::Coeff` (associated types of
/// [`Accumulate`]), so there is no phantom axis and no way to pair a storage type
/// with the wrong key.
///
/// `PauliSum = Sum<HashMapStore<PauliWord, C>, P>` is the domain alias (see the
/// crate root).
#[derive(Debug, Clone)]
pub struct Sum<S, P = crate::policy::NoPolicy>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// The graded-algebra container (`Vec<(K, C)>` or
    /// `HashMap<K, C, IdentityBuildHasher>`).
    storage: S,
    /// The truncation policy.
    policy: P,
    /// Invariant: every key `k` in `storage` satisfies `k.n_sites() == n_sites`.
    /// An empty sum has no key to derive the width from, so the field is carried
    /// explicitly and checked by a `debug_assert!` on every insertion path.
    n_sites: usize,
}

// --- Inspection: needs only `Support` (a supertrait of `Accumulate`). --------
impl<S, P> Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Number of qubits/sites this sum is defined over.
    #[inline]
    pub fn n_sites(&self) -> usize {
        self.n_sites
    }

    /// Number of terms in the (reduced) support.
    #[inline]
    pub fn len(&self) -> usize {
        self.storage.len()
    }

    /// Whether the support is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.storage.is_empty()
    }

    /// The coefficient at `key`, if present.
    #[inline]
    pub fn get(&self, key: &S::Key) -> Option<S::Coeff> {
        self.storage.get(key)
    }

    /// Whether `key` is in the support.
    #[inline]
    pub fn contains(&self, key: &S::Key) -> bool {
        self.storage.get(key).is_some()
    }

    /// Read-only export of the support as `(key, coeff)` pairs.
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (S::Key, S::Coeff)> + '_ {
        self.storage.iter()
    }

    /// Borrow the truncation policy.
    #[inline]
    pub fn policy(&self) -> &P {
        &self.policy
    }
}

// --- The `C`-module action (L2): needs `Scale`. -----------------------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + Scale,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Multiply every coefficient by `s` (`∀ k. cₖ *= s`).
    #[inline]
    pub fn scale(&mut self, s: &S::Coeff) {
        self.storage.scale(s);
    }
}

// --- The trace pairings (L3): needs `Pair`. ---------------------------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + Pair,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// The symmetric bilinear Hilbert–Schmidt trace pairing
    /// `⟨self, other⟩ = Σ_k self[k]·other[k]`.
    ///
    /// Bilinear (not conjugated); machine-checked biadditive/symmetric/homogeneous
    /// in `lean/PPVM/Algebra/GradedMap.lean` (`overlap_add_left`/`_right`,
    /// `overlap_comm`, `overlap_smul_left`/`_right`). Delegates to the storage's
    /// [`Pair::overlap`].
    #[inline]
    pub fn overlap(&self, other: &Self) -> S::Coeff {
        self.storage.overlap(&other.storage)
    }

    /// The sesquilinear inner product `⟨self | other⟩ = Σ_k conj(self[k])·other[k]`.
    /// Requires the coefficient ring to carry a [`Conjugate`] involution.
    #[inline]
    pub fn hermitian_overlap(&self, other: &Self) -> S::Coeff
    where
        S::Coeff: Conjugate,
    {
        self.storage.hermitian_overlap(&other.storage)
    }
}

// --- Construction: needs `StoreAlloc` to size/allocate storage. --------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + StoreAlloc,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// An empty sum on `n_sites` sites, with the default policy and a capacity
    /// hint drawn from that policy.
    #[inline]
    pub fn new(n_sites: usize) -> Self {
        let policy = P::default();
        let storage = S::with_capacity(policy.capacity(n_sites));
        Self {
            storage,
            policy,
            n_sites,
        }
    }

    /// An empty sum on `n_sites` sites with an explicit `policy`.
    #[inline]
    pub fn with_policy(n_sites: usize, policy: P) -> Self {
        let storage = S::with_capacity(policy.capacity(n_sites));
        Self {
            storage,
            policy,
            n_sites,
        }
    }

    /// Build a sum on `n_sites` sites from an iterator of `(key, coeff)` terms,
    /// using the default policy. Colliding keys are combined and zero
    /// coefficients dropped (`reduce`).
    pub fn from_terms<I>(n_sites: usize, terms: I) -> Self
    where
        I: IntoIterator<Item = (S::Key, S::Coeff)>,
    {
        Self::from_terms_with_policy(n_sites, P::default(), terms)
    }

    /// Build a sum on `n_sites` sites from an iterator of `(key, coeff)` terms
    /// with an explicit `policy`. Colliding keys are combined and zero
    /// coefficients dropped (`reduce`); the policy's truncation is **not** run
    /// here (call [`Sum::truncate`] to apply it).
    pub fn from_terms_with_policy<I>(n_sites: usize, policy: P, terms: I) -> Self
    where
        I: IntoIterator<Item = (S::Key, S::Coeff)>,
    {
        let mut sum = Self::with_policy(n_sites, policy);
        let mut batch = TermBatch::new();
        for (k, c) in terms {
            debug_assert_eq!(
                k.n_sites(),
                n_sites,
                "term key width must match the sum's n_sites"
            );
            batch.push(k, c);
        }
        sum.storage.accumulate_batch(&batch);
        sum.storage.reduce();
        sum
    }
}

// --- The engine step: needs `StoreAlloc` (reset) and `Retain` (truncate). ----
impl<S, P> Sum<S, P>
where
    S: Accumulate + StoreAlloc + Retain<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Apply a term producer: read the current support through `&`, produce the
    /// transformed terms into a batch, then **replace** the support with the
    /// accumulated batch, canonicalize (`reduce`), and truncate.
    ///
    /// The support is reset between producing and accumulating (see
    /// [`StoreAlloc`] — the design's `apply` sketch glosses this, which would
    /// double-count under a bijective re-key). `TP` is a type parameter, never
    /// `dyn`, so the call monomorphizes and the producer's `#[inline]` `produce`
    /// folds into the loop.
    ///
    /// Design: §"apply" and §"Every gate is a producer feeding `accumulate`".
    pub fn apply<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<S::Key, S::Coeff>,
    {
        let mut batch = TermBatch::with_capacity(self.storage.len());
        for (k, c) in self.storage.iter() {
            producer.produce(&k, &c, &mut batch);
        }
        // Replace, not merge: a bijective re-key must not retain the old keys.
        self.storage.reset();
        self.storage.accumulate_batch(&batch);
        self.storage.reduce();
        self.policy.truncate(&mut self.storage);
    }

    /// Run the configured policy's truncation on the current support.
    #[inline]
    pub fn truncate(&mut self) {
        self.policy.truncate(&mut self.storage);
    }
}

// --- The Clifford fast path: needs `RekeyBijective` (move-based re-key). ------
impl<S, P> Sum<S, P>
where
    S: Accumulate + StoreAlloc + Retain<S::Key, S::Coeff> + RekeyBijective<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Re-key every term by a **bijection** `f: (k, c) ↦ (φ(k), c')` in place —
    /// moving each term through `f` with no key or coefficient clones and reusing
    /// the backing allocation — then run the policy's truncation.
    ///
    /// This is the fast path for a Clifford conjugation (a Pauli bijection). It
    /// deliberately bypasses [`apply`](Self::apply)'s batch round-trip, whose
    /// read-side `iter()` clone and merge-side `entry(k.clone())` clone cost
    /// against the old crate's in-place aux-map swap
    /// ([`RekeyBijective`](crate::store::RekeyBijective) friction note). Because
    /// `f` is injective on keys, no re-keyed terms collide and `reduce` is
    /// unnecessary — a `±1` sign never zeroes a coefficient — but truncation still
    /// runs, since a Clifford can change a key's Pauli weight (e.g. `CNOT`:
    /// `IX ↦ XX`).
    ///
    /// Design: §"apply" and §"Pauli algebra traits" (the sum "applies the one-row
    /// action pointwise and drains each term's phase delta to its coefficient").
    #[inline]
    pub(crate) fn rekey_bijective<F>(&mut self, f: F)
    where
        F: FnMut(S::Key, S::Coeff) -> (S::Key, S::Coeff),
    {
        self.storage.rekey_bijective(f);
        self.policy.truncate(&mut self.storage);
    }
}

// --- The pure-sign Clifford fast path: needs `SignFlipByKey` (in-place scale). -
impl<S, P> Sum<S, P>
where
    S: Accumulate + SignFlipByKey<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Multiply every term's coefficient in place by `f(&k) ∈ {−1, +1}`, keyed on
    /// that term's own key — **no** map rebuild, no key movement, no reallocation.
    ///
    /// This is the fast path for a *pure-sign* Clifford conjugation (`X`/`Y`/`Z`),
    /// whose Pauli word is fixed so only the coefficient sign changes. It skips the
    /// move-based [`rekey_bijective`](Self::rekey_bijective) entirely — there is
    /// nothing to re-key — restoring the old crate's in-place `scale`. A `±1` never
    /// zeroes a coefficient and never changes a key's magnitude or Pauli weight, so
    /// neither `reduce` nor the policy's truncation can act; both are skipped.
    ///
    /// Design: §"Pauli algebra traits" (the sum "applies the one-row action
    /// pointwise and drains each term's phase delta to its coefficient").
    #[inline]
    pub(crate) fn flip_sign_by_key<F>(&mut self, f: F)
    where
        F: Fn(&S::Key) -> i8,
    {
        self.storage.sign_flip_by_key(f);
    }
}

// --- The diagonal-noise fast path: needs `ScaleByKey` (in-place per-key scale). -
impl<S, P> Sum<S, P>
where
    S: Accumulate + ScaleByKey<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Multiply every term's coefficient in place by the per-key ring factor
    /// `f(&k)`, keyed on that term's own key — **no** map rebuild, no key
    /// movement, no reallocation. A `None` result leaves the term untouched.
    ///
    /// This is the fast path for a *diagonal* unital Pauli channel
    /// ([`PauliError`](ppvm_traits_2::PauliError)), whose Pauli words are fixed so
    /// only the coefficients pick up the channel's real transfer eigenvalue. Like
    /// [`flip_sign_by_key`](Self::flip_sign_by_key) it skips the move-based
    /// [`rekey_bijective`](Self::rekey_bijective) entirely — there is nothing to
    /// re-key — restoring the old crate's in-place `scale`. The channel is
    /// contractive (`|λ_P| ≤ 1`), so it can never grow a key's Pauli weight; a
    /// term the channel scales to *exactly* zero (a zero eigenvalue, e.g.
    /// `pauli_error(q, [0.0, 0.25, 0.25])` → `λ_X = 0`) is dropped inside
    /// [`ScaleByKey`](crate::store::ScaleByKey) so the reduced-canonical-form
    /// invariant holds, but no separate whole-map `reduce`/truncation pass runs.
    ///
    /// Design: §"Behavioral traits" (`PauliError`); the eigenvalue is
    /// machine-checked in `lean/PPVM/Algebra/Noise.lean`
    /// (`pauli_channel_eigenvalue`).
    #[inline]
    pub(crate) fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&S::Key) -> Option<S::Coeff>,
    {
        self.storage.scale_by_key(f);
    }
}

// --- The rotation fast path: needs `RotateInPlace` (fused branching re-key). ---
impl<S, P> Sum<S, P>
where
    S: Accumulate + RotateInPlace<S::Key, S::Coeff> + Retain<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Propagate a single-qubit rotation through the support in one fused pass:
    /// `f(&k, &mut c)` scales each diagonal coefficient **in place** (the `cosθ`
    /// factor) and returns the optional anticommuting branch term `(iGP, c·sinθ·ε)`
    /// to merge, then the policy's truncation runs.
    ///
    /// This is the branching analogue of the pure-sign
    /// [`flip_sign_by_key`](Self::flip_sign_by_key) / diagonal
    /// [`scale_by_key`](Self::scale_by_key) fast paths: it deliberately bypasses
    /// [`apply`](Self::apply)'s batch round-trip, whose read-side `iter()` key
    /// clone, `reset`, and re-accumulation of the *whole* fan-out re-hash the `N`
    /// untouched diagonals. Here each diagonal is scaled where it sits — cached
    /// hash intact, no bucket move — and only the ≤`N` branch terms are hashed and
    /// merged, restoring the old crate's single-pass `map_insert`
    /// ([`RotateInPlace`](crate::store::RotateInPlace) friction note). Exact-zero
    /// cancellations from a colliding merge are dropped inside the walk (the
    /// `reduce` `apply` would run).
    ///
    /// Design: §"apply", §"Every gate is a producer feeding `accumulate`", and
    /// §"Behavioral traits" (`RotationOne`).
    #[inline]
    pub(crate) fn rotate_in_place<F>(&mut self, f: F)
    where
        F: FnMut(&S::Key, &mut S::Coeff) -> Option<(S::Key, S::Coeff)>,
    {
        self.storage.rotate_in_place(f);
        self.policy.truncate(&mut self.storage);
    }
}
