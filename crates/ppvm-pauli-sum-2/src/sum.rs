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

use std::collections::HashSet;

use ppvm_traits_2::{
    Accumulate, Conjugate, IdentityBuildHasher, Indexable, Pair, Retain, Scale, TermBatch,
    TermProducer, TermSink, Word,
};

use crate::policy::Policy;
use crate::store::{
    AddTerm, ApplyProducer, RekeyBijective, RotateInPlace, ScaleByKey, SignFlipByKey, StoreAlloc,
};

/// The keep-set carried by a [`Sum`]: keys [`Sum::truncate`] restores if the
/// policy dropped them. Hashed through the same pass-through
/// [`IdentityBuildHasher`] as the support, so a membership test costs no hashing.
pub type PreserveSet<K> = HashSet<K, IdentityBuildHasher>;

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
#[derive(Clone)]
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
    /// The capacity hint the storage was sized from — the policy's
    /// [`Policy::capacity`] or the caller's explicit override. Reported by
    /// [`capacity`](Self::capacity), matching old's `PauliSum::capacity()`, which
    /// likewise returns the *requested* hint rather than the map's current
    /// allocation.
    capacity: usize,
    /// Keys [`truncate`](Self::truncate) restores if the policy dropped them —
    /// old's `PauliSum::preserve_strings`. **Empty by default**, and the empty
    /// case is a hard fast path (see [`truncate`](Self::truncate)): an empty
    /// `HashSet` allocates nothing, so an opt-out user pays one predictable
    /// branch per `truncate` and no memory.
    preserve: PreserveSet<S::Key>,
}

/// Hand-written rather than derived so the keep-set does not force
/// `S::Key: Debug` on the whole engine (the derive would make every `Sum` with a
/// non-`Debug` key un-printable). The set is summarized by size; the storage,
/// policy, width and capacity hint print as before.
impl<S, P> std::fmt::Debug for Sum<S, P>
where
    S: Accumulate + std::fmt::Debug,
    P: Policy<S::Key, S::Coeff> + std::fmt::Debug,
    S::Key: Word + Indexable,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Sum")
            .field("storage", &self.storage)
            .field("policy", &self.policy)
            .field("n_sites", &self.n_sites)
            .field("capacity", &self.capacity)
            .field("preserved_keys", &self.preserve.len())
            .finish()
    }
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

    /// The capacity hint the storage was sized from.
    ///
    /// Old's `PauliSum::capacity()`: the value the builder resolved (the
    /// strategy's hint, or the explicit `.capacity(..)` override), **not** the
    /// map's live allocation.
    #[inline]
    pub fn capacity(&self) -> usize {
        self.capacity
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

    /// Whether the support contains `key` **with exactly this coefficient**.
    ///
    /// This is old's `PauliSum::contains(&key, &value)`
    /// (`contains_with(k, |x| x == v)`): a `(key, value)` match, *not* a
    /// membership test. The membership test is [`contains_key`](Self::contains_key)
    /// — deliberately a different name, so a caller ported from old cannot get a
    /// silently different predicate out of the same spelling.
    #[inline]
    pub fn contains(&self, key: &S::Key, value: &S::Coeff) -> bool {
        self.storage.get(key).is_some_and(|c| c == *value)
    }

    /// Whether `key` is in the support, at any coefficient (including zero).
    #[inline]
    pub fn contains_key(&self, key: &S::Key) -> bool {
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

    /// Read-only access to the active keep-set — old's
    /// `PauliSum::preserve_strings()`.
    #[inline]
    pub fn preserved_keys(&self) -> &PreserveSet<S::Key> {
        &self.preserve
    }

    /// Install a keep-set: keys [`truncate`](Self::truncate) restores (at their
    /// pre-truncate coefficient) whenever the policy would drop them.
    ///
    /// This is old's builder step `PauliSum::builder().preserve_strings(set)`,
    /// in the by-value builder shape the rest of this constructor family uses:
    ///
    /// ```
    /// use ppvm_pauli_sum_2::{CoefficientThreshold, PauliSum, PauliWord};
    ///
    /// let sum = PauliSum::<f64, CoefficientThreshold>::with_policy(
    ///     3,
    ///     CoefficientThreshold { threshold: 0.5 },
    /// )
    /// .preserving(["ZII", "IZI", "IIZ"].map(PauliWord::from));
    /// assert_eq!(sum.preserved_keys().len(), 3);
    /// ```
    ///
    /// The mechanism composes with **any** policy, because the policy still runs
    /// unchanged in the middle of [`truncate`](Self::truncate) — see there for
    /// the exact snapshot/restore semantics.
    #[must_use]
    pub fn preserving<I>(mut self, keys: I) -> Self
    where
        I: IntoIterator<Item = S::Key>,
    {
        self.preserve.extend(keys);
        self
    }

    /// Add one key to the keep-set — the mutating form of
    /// [`preserving`](Self::preserving), for a caller assembling the set after
    /// construction.
    #[inline]
    pub fn preserve_key(&mut self, key: S::Key) {
        debug_assert_eq!(
            key.n_sites(),
            self.n_sites(),
            "preserved key width must match the sum's n_sites"
        );
        self.preserve.insert(key);
    }

    /// Empty the keep-set, restoring the bare-policy [`truncate`](Self::truncate)
    /// fast path.
    #[inline]
    pub fn clear_preserved_keys(&mut self) {
        self.preserve.clear();
    }

    /// Canonicalize the support to reduced finite-support form: drop every term
    /// whose coefficient is **exactly** zero.
    ///
    /// This is [`Accumulate::reduce`] surfaced on the engine — a *caller-driven*
    /// finalize step, deliberately **not** run by any gate, channel, re-key or by
    /// the L4 [`multiply_into`](Self::multiply_into) (Design: §"`reduce()` is
    /// first-class, and runs only at finalize"; several `accumulate_batch` calls
    /// or a whole outer product can precede one `reduce`). Nothing calls it
    /// implicitly, so the old crate's contract that a zero-coefficient term
    /// survives every operation is untouched: dropping a zero happens only when
    /// the caller asks for it here.
    ///
    /// Lean: `reduce_structural` in `lean/PPVM/Algebra/GradedMap.lean` — the
    /// reduced map agrees with the original everywhere and its support is exactly
    /// the non-zero keys.
    #[inline]
    pub fn reduce(&mut self) {
        self.storage.reduce();
    }

    /// Borrow the storage backend — **crate-internal**, for the sibling
    /// capability modules (`clifford`, `rotation`, `noise`, `multiply`) that
    /// drive a storage fast path the engine does not otherwise expose.
    ///
    /// The backing container stays private to the crate: nothing here widens the
    /// public surface, which remains the `ppvm-traits-2` trait impls plus the
    /// inherent methods above.
    #[inline]
    pub(crate) fn storage(&self) -> &S {
        &self.storage
    }

    /// Mutably borrow the storage backend — **crate-internal**; see
    /// [`storage`](Self::storage).
    #[inline]
    pub(crate) fn storage_mut(&mut self) -> &mut S {
        &mut self.storage
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
        Self::with_policy(n_sites, P::default())
    }

    /// An empty sum on `n_sites` sites with an explicit `policy`; the capacity
    /// hint is `policy.capacity(n_sites)`.
    #[inline]
    pub fn with_policy(n_sites: usize, policy: P) -> Self {
        let capacity = policy.capacity(n_sites);
        Self::with_capacity(n_sites, policy, capacity)
    }

    /// An empty sum on `n_sites` sites with an explicit `policy` **and an
    /// explicit capacity**, overriding the policy's hint.
    ///
    /// This is old's builder override — `PauliSum::builder().n_qubits(n)
    /// .strategy(s).capacity(c).build()` — which every real workload uses
    /// (`benches/trotter.rs` and `examples/trotter_qubit_sweep.rs` pass
    /// `n_qubits.pow(2)`, `benches/gate.rs` passes `1 << 20`). Sizing the store
    /// up front is what keeps a long propagation from rehashing a ~10⁴-term map
    /// mid-circuit, and the policy hint alone cannot express a workload's own
    /// knowledge of its support; without this the engine could only be run at
    /// `policy.capacity(n_sites)`.
    ///
    /// The hint reaches [`StoreAlloc::with_capacity`], which sizes the primary
    /// map, the auxiliary double-buffer and the branch scratch from it.
    #[inline]
    pub fn with_capacity(n_sites: usize, policy: P, capacity: usize) -> Self {
        let storage = S::with_capacity(capacity);
        Self {
            storage,
            policy,
            n_sites,
            capacity,
            // Old's builder default: `preserve_strings = {}`. An empty `HashSet`
            // does not allocate, so the opt-in machinery costs nothing here.
            preserve: PreserveSet::with_hasher(IdentityBuildHasher),
        }
    }

    /// Build a sum on `n_sites` sites from an iterator of `(key, coeff)` terms,
    /// using the default policy. Colliding keys are **combined**.
    pub fn from_terms<I>(n_sites: usize, terms: I) -> Self
    where
        I: IntoIterator<Item = (S::Key, S::Coeff)>,
    {
        Self::from_terms_with_policy(n_sites, P::default(), terms)
    }

    /// Build a sum on `n_sites` sites from an iterator of `(key, coeff)` terms
    /// with an explicit `policy`. Colliding keys are **combined**.
    ///
    /// Neither the policy's truncation nor [`reduce`](Self::reduce) runs here:
    /// this is exactly `n` applications of old's `sum += (word, coeff)`, which
    /// keeps a term whose coefficient is (or cancels to) exactly zero. Call
    /// [`Sum::truncate`] for the policy, [`Sum::reduce`] to canonicalize.
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
        sum
    }
}

// --- The engine step: needs the store's producer path. -----------------------
impl<S, P> Sum<S, P>
where
    S: Accumulate + ApplyProducer<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Apply a term producer: read the current support through `&`, produce the
    /// transformed terms into a batch, then **replace** the support with the
    /// accumulated batch.
    ///
    /// **Truncation is not run here.** The policy acts only when the caller
    /// invokes [`truncate`](Self::truncate) — the old crate's *deferred*
    /// truncation semantics (`ppvm-pauli-sum/tests/truncation_semantics.rs`):
    /// a gate must never drop a term on its own, because two individually
    /// sub-threshold contributions to the same key can merge into a surviving
    /// one (`|a + b| ≥ τ` while `|a|, |b| < τ`) when several gates run between
    /// truncations.
    ///
    /// **[`reduce`](Self::reduce) is not run here either.** Old has no `reduce`
    /// anywhere: a term whose coefficient is (or cancels to) exactly zero stays
    /// in the support, and old's exact-map `PartialEq` counts it
    /// (`ppvm-pauli-sum/tests/loss.rs::test_reset_channel` asserts a sum whose
    /// coefficients were all zeroed still equals the original key set).
    /// Canonicalization is caller-driven — [`Sum::reduce`].
    ///
    /// The support is reset between producing and accumulating (see
    /// [`ApplyProducer`] — the design's `apply` sketch glosses this, which would
    /// double-count under a bijective re-key). `TP` is a type parameter, never
    /// `dyn`, so the call monomorphizes and the producer's `#[inline]` `produce`
    /// folds into the loop.
    ///
    /// The produced-term batch is the **storage's**, not a fresh one per call:
    /// `Sum` owns no workspace, so the reusable batch lives beside the `aux`
    /// double-buffer and rotation `scratch` in the backend
    /// ([`ApplyProducer`]) — otherwise every producer-based gate would pay two
    /// `Vec` allocations, the per-gate churn architecture features 1–2 exist to
    /// remove.
    ///
    /// Design: §"apply" (incl. the "driver-owned reusable batch" note) and
    /// §"Every gate is a producer feeding `accumulate`".
    #[inline]
    pub fn apply<TP>(&mut self, producer: TP)
    where
        TP: TermProducer<S::Key, S::Coeff>,
    {
        self.storage.apply_producer(producer);
    }
}

// --- Truncation: needs `Retain` (the policy) and `AddTerm` (preserve restore). -
impl<S, P> Sum<S, P>
where
    S: Accumulate + Retain<S::Key, S::Coeff> + AddTerm<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Run the configured policy's truncation on the current support, then
    /// restore any **preserved** key the policy dropped.
    ///
    /// This is the **only** way terms leave the support by policy: gates never
    /// truncate on their own (see [`apply`](Self::apply)).
    ///
    /// # The keep-set post-filter
    ///
    /// Old's `PauliSum::truncate` (`ppvm-pauli-sum/src/sum/data.rs:271`), ported
    /// step for step:
    ///
    /// 1. **Empty keep-set → run the policy and nothing else.** This is a hard
    ///    fast path, not an optimization detail: `truncate()` is the most-called
    ///    method in the headline Trotter workload (~1500 calls over the circuit)
    ///    and the keep-set is empty in every benchmark configuration, so the
    ///    preserve machinery must cost one branch and no walk when unused.
    /// 2. Otherwise **snapshot** the preserved keys' *pre-truncate* coefficients.
    /// 3. Run the configured policy **verbatim** — which is what makes this
    ///    compose with any policy (magnitude floor, weight cap, combinations).
    /// 4. **Re-insert** every preserved key the policy dropped, at its
    ///    pre-truncate coefficient, guarded by a membership test so a survivor is
    ///    never double-added (old's `contains_with(&k, |_| true)` guard; the
    ///    restore goes through `add_term`, which would otherwise *sum* into a
    ///    kept entry).
    ///
    /// A preserved key that was **not** in the support before truncation is not
    /// inserted by this: the snapshot only records keys that were there, exactly
    /// as old's scan does.
    ///
    /// The one implementation difference from old is the snapshot: old walks the
    /// **whole support** with an always-true `retain` (it had no other route to
    /// `(k, v)` pairs through its map traits) and tests each key for membership;
    /// this probes the `|keep-set|` keys directly. Same snapshot, same restores,
    /// `O(|keep-set|)` instead of `O(|support|)` — and it is why the keep-set is
    /// hashed through the pass-through [`IdentityBuildHasher`], so a probe costs
    /// no hashing either.
    ///
    /// Lean: `truncate_preserve_eq_widened_retain` in
    /// `lean/PPVM/Algebra/Truncation.lean` — the three-step composite equals a
    /// single pass with the **widened** keep-rule `keep(k, c) ∨ k ∈ P`, which is
    /// what pins both otherwise-invisible guards (without the membership test the
    /// accumulating `add_term` would *double* a survivor; with the snapshot taken
    /// after the policy it would restore a post-truncate residue). Corollaries:
    /// `truncatePreserve_apply_of_mem` (a preserved key keeps *exactly* its
    /// pre-truncate coefficient — old's `Σᵢ Zᵢ` conservation test) and
    /// `truncatePreserve_empty` (the empty-keep-set fast path is exact).
    pub fn truncate(&mut self) {
        // 1. Hot path: no keep-set → just the policy (no snapshot scan).
        if self.preserve.is_empty() {
            self.policy.truncate(&mut self.storage);
            return;
        }

        // 2. Snapshot the preserved keys' pre-truncate coefficients.
        let mut saved: Vec<(S::Key, S::Coeff)> = Vec::with_capacity(self.preserve.len());
        for key in &self.preserve {
            if let Some(coeff) = self.storage.get(key) {
                saved.push((key.clone(), coeff));
            }
        }

        // 3. The configured policy, verbatim.
        self.policy.truncate(&mut self.storage);

        // 4. Restore whatever it dropped. The membership guard keeps `add_term`
        //    (an accumulate) from doubling a survivor.
        for (key, coeff) in saved {
            if self.storage.get(&key).is_none() {
                self.storage.add_term(key, coeff);
            }
        }
    }
}

// --- The Clifford fast path: needs `RekeyBijective` (move-based re-key). ------
impl<S, P> Sum<S, P>
where
    S: Accumulate + RekeyBijective<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Re-key every term by a **bijection** `f: (k, c) ↦ (φ(k), c')` in place —
    /// moving each term through `f` with no key or coefficient clones and reusing
    /// the backing allocation.
    ///
    /// This is the fast path for a Clifford conjugation (a Pauli bijection). It
    /// deliberately bypasses [`apply`](Self::apply)'s batch round-trip, whose
    /// read-side `iter()` clone and merge-side `entry(k.clone())` clone cost
    /// against the old crate's in-place aux-map swap
    /// ([`RekeyBijective`](crate::store::RekeyBijective) friction note). Because
    /// `f` is injective on keys, no re-keyed terms collide and `reduce` is
    /// unnecessary — a `±1` sign never zeroes a coefficient.
    ///
    /// **Truncation is not run here**, matching the old crate's deferred
    /// semantics — even though a Clifford *can* change a key's Pauli weight
    /// (`CNOT`: `IX ↦ XX`) and so make a term newly truncatable. Acting on that
    /// is the caller's [`truncate`](Self::truncate).
    ///
    /// Design: §"apply" and §"Pauli algebra traits" (the sum "applies the one-row
    /// action pointwise and drains each term's phase delta to its coefficient").
    #[inline]
    pub(crate) fn rekey_bijective<F>(&mut self, f: F)
    where
        F: FnMut(S::Key, S::Coeff) -> (S::Key, S::Coeff),
    {
        self.storage.rekey_bijective(f);
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
    /// Mutate every term's coefficient in place via `f(&k, &mut c)`, keyed on
    /// that term's own key — **no** map rebuild, no key movement, no
    /// reallocation. A term `f` does not touch is left exactly as it was.
    ///
    /// This is the fast path for a *diagonal* unital Pauli channel
    /// ([`PauliError`](ppvm_traits_2::PauliError)), whose Pauli words are fixed so
    /// only the coefficients pick up the channel's real transfer eigenvalue. Like
    /// [`flip_sign_by_key`](Self::flip_sign_by_key) it skips the move-based
    /// [`rekey_bijective`](Self::rekey_bijective) entirely — there is nothing to
    /// re-key — restoring the old crate's in-place `scale`. The channel is
    /// contractive (`|λ_P| ≤ 1`), so it can never grow a key's Pauli weight. A
    /// term the channel scales to *exactly* zero (a zero eigenvalue, e.g.
    /// `pauli_error(q, [0.0, 0.25, 0.25])` → `λ_X = 0`) **stays in the support**
    /// with coefficient zero — the old crate has no `reduce` and never removes a
    /// zero term, and its exact-map equality depends on that. No whole-map
    /// `reduce`/truncation pass runs here.
    ///
    /// Design: §"Behavioral traits" (`PauliError`); the eigenvalue is
    /// machine-checked in `lean/PPVM/Algebra/Noise.lean`
    /// (`pauli_channel_eigenvalue`), and so is the contractivity this fast path
    /// rests on: `|λ_P| ≤ 1` for a sub-stochastic `[p_X, p_Y, p_Z]`
    /// (`pauli_channel_eigenvalue_abs_le_one`), hence an `ℓ¹` contraction
    /// (`l1_contractive`, which is what makes `Truncation.l1_bound` compose over
    /// a noisy circuit) whose support never gains or moves a key
    /// (`scaleByKey_support_subset`) — so no weight re-check is owed. An
    /// over-normalized probability vector voids both claims
    /// (`eigenvalue_abs_le_one_needs_substochastic`).
    #[inline]
    pub(crate) fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&S::Key, &mut S::Coeff),
    {
        self.storage.scale_by_key(f);
    }
}

// --- The rotation fast path: needs `RotateInPlace` (fused branching re-key). ---
impl<S, P> Sum<S, P>
where
    S: Accumulate + RotateInPlace<S::Key, S::Coeff>,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + Indexable,
{
    /// Propagate a single-qubit rotation through the support in one fused pass:
    /// `f(&k, &mut c)` scales each diagonal coefficient **in place** (the `cosθ`
    /// factor) and returns the optional anticommuting branch term `(iGP, c·sinθ·ε)`
    /// to merge.
    ///
    /// **Truncation is not run here** — a freshly branched term below the
    /// coefficient threshold must survive to be merged by later gates (the exact
    /// case `ppvm-pauli-sum/tests/truncation_semantics.rs` pins: two `rx(θ)` with
    /// no truncate between them). The policy acts only on
    /// [`truncate`](Self::truncate).
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
    /// branches and exact-zero cancellations are **kept**, as in old's
    /// `map_insert`, which merges every produced term through `add_assign`.
    ///
    /// Lean: `accumulate_rotBatch` in `lean/PPVM/Instantiations/Rotation.lean` —
    /// the two-pass fast path (scale **all** diagonals, *then* merge **all**
    /// branches) equals folding the whole `≤2N`-term produced batch into the map
    /// through `accumulate_batch`, for every walk order. The per-term
    /// `anticommute_new_key` does *not* cover this: a branch key can collide with
    /// a **different**, not-yet-scaled key (`rx` on `{Z, Y}` at one site swaps
    /// them), so a variant that merged eagerly inside the walk would be wrong —
    /// `eagerWalk_ne_twoPass` exhibits exactly that divergence.
    ///
    /// Design: §"apply", §"Every gate is a producer feeding `accumulate`", and
    /// §"Behavioral traits" (`RotationOne`).
    #[inline]
    pub(crate) fn rotate_in_place<F>(&mut self, f: F)
    where
        F: FnMut(&S::Key, &mut S::Coeff) -> Option<(S::Key, S::Coeff)>,
    {
        self.storage.rotate_in_place(f);
    }
}
