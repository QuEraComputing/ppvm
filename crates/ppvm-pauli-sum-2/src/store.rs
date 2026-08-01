// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The storage-lifecycle helper the [`Sum`](crate::Sum) engine needs *beyond*
//! the graded algebra: allocate a fresh empty store at a capacity hint, and
//! reset a store to empty support **keeping its allocation**.
//!
//! # Friction: `apply` needs a clear/replace the graded traits do not expose
//!
//! The design's `apply` sketch (`traits-2-configuration-and-hashing.md`
//! §"apply") reads:
//!
//! ```text
//! for (k, c) in self.storage.iter() { producer.produce(&k, &c, &mut batch); }
//! self.storage.accumulate_batch(&batch);   // <-- merges onto the *existing* support
//! self.storage.reduce();
//! ```
//!
//! For a **bijective Clifford re-key** the producer emits `(φ(k), c)` for every
//! stored `(k, c)`; merging that batch onto the un-cleared support leaves both
//! the old `k` *and* the new `φ(k)`, double-counting. The old
//! `ppvm-pauli-sum::map_add` avoided this by writing into a **cleared** auxiliary
//! map and swapping. `apply` must reproduce that "replace, not merge" semantics,
//! so it resets `self.storage` between producing and accumulating.
//!
//! [`Accumulate`](ppvm_traits_2::Accumulate) exposes no `clear`, so this trait
//! supplies one. It is a **local** trait `impl`'d on the foreign `Vec`/`HashMap`
//! containers (legal: the trait is local), mirroring the graded-algebra impls in
//! `ppvm-traits-2`. `reset` clears in place so the allocation is reused across
//! millions of gates — matching the old aux-map reuse and avoiding a per-gate
//! reallocation.

use std::collections::HashMap;

use ppvm_traits_2::{
    Accumulate, Coefficient, Conjugate, IdentityBuildHasher, Indexable, KeyBatch, Pair, Retain,
    Scale, Support, TermBatch,
};

/// The provided hash-join storage backend for [`Sum`](crate::Sum): a primary
/// support map plus a persistent **auxiliary** map and **scratch** buffer, all
/// keyed through the pass-through [`IdentityBuildHasher`] (so a key's finalized
/// `key_hash()` digest reaches hashbrown untouched — Design: §"The pass-through
/// storage contract").
///
/// This relocates the old `ppvm-pauli-sum` double-buffer (`PauliSum` held
/// `map: (primary, aux)` + a reusable `scratch` Vec) **into the storage backend**,
/// so `Sum` stays a pure engine over an abstract store — and the same backend
/// composes into the generalized tableau — while the hot fast paths still reuse
/// their buffers across millions of gates instead of reallocating per gate:
///
/// * [`RekeyBijective`] (Clifford re-key) drains the primary through the re-key
///   closure into the cleared `aux`, then swaps — the old `map_add`'s
///   clear→write→swap: zero per-gate allocation, zero key clones. The two map
///   allocations ping-pong across gates and are never freed.
/// * [`RotateInPlace`] (rotation branch) buffers the ≤`N` branch terms in the
///   persistent `scratch` before merging — the old `map_insert`'s reused scratch.
///
/// `aux` and `scratch` are **transient**: both are empty between operations
/// (cleared on entry; left empty by the drain/swap on exit), so [`Clone`] copies
/// only the primary support and equality/iteration observe only it.
#[derive(Debug)]
pub struct HashMapStore<K, C> {
    /// The reduced-canonical support the graded traits read and write.
    primary: HashMap<K, C, IdentityBuildHasher>,
    /// Persistent double-buffer for [`RekeyBijective`] (clear → write → swap).
    /// Empty between operations.
    aux: HashMap<K, C, IdentityBuildHasher>,
    /// Persistent branch/term scratch for [`RotateInPlace`]. Empty between
    /// operations.
    scratch: Vec<(K, C)>,
}

impl<K: Clone, C: Clone> Clone for HashMapStore<K, C> {
    /// Clone only the primary support; `aux`/`scratch` are transient (empty
    /// between ops), so a clone starts them fresh rather than copying scratch.
    #[inline]
    fn clone(&self) -> Self {
        Self {
            primary: self.primary.clone(),
            aux: HashMap::with_hasher(IdentityBuildHasher),
            scratch: Vec::new(),
        }
    }
}

/// Allocate and reset a graded-sum storage container. Not a graded-algebra
/// operation — see the module docs on why `apply` needs the reset the graded
/// traits omit.
pub trait StoreAlloc {
    /// A fresh, empty store pre-sized for roughly `cap` terms.
    fn with_capacity(cap: usize) -> Self;

    /// Reset to empty support, keeping the backing allocation for reuse.
    fn reset(&mut self);
}

/// An in-place, **move-based** bijective re-key — the fast path a Clifford
/// conjugation drives instead of the [`apply`](crate::Sum::apply) batch
/// round-trip.
///
/// # Friction: `apply`'s batch clones every key twice; a Clifford cannot afford it
///
/// The generic `apply` reads the support through [`Support::iter`], which
/// **clones** every `(k, c)`, produces the re-keyed terms into a `TermBatch`,
/// then merges the batch with [`Accumulate::accumulate_batch`], which clones
/// each key a **second** time (`entry(k.clone())`). For a Pauli-keyed sum the new
/// word dropped `Copy` (its lazy `OnceLock` hash cache), so each clone is a real
/// cost; against the old crate's single in-place aux-map swap
/// (`ppvm-pauli-sum::PauliSum::map_add`, one clone per term) the batch path
/// benchmarked ~4.7× slower on a Clifford gate.
///
/// A Clifford conjugation is a **bijection** on Pauli words, so the transformed
/// support is just the old support with each key relabelled and its coefficient
/// scaled by `±1`: no new keys collide and no coefficient is zeroed. That lets us
/// skip the batch entirely and rebuild the support by **moving** each term
/// through the re-key closure — zero key clones, zero coefficient clones — while
/// reusing the backing allocation. This is a `ppvm-pauli-sum-2`-local trait
/// `impl`'d on the foreign `Vec`/`HashMap` containers (legal: the trait is
/// local), mirroring [`StoreAlloc`] and the graded-algebra impls in
/// `ppvm-traits-2`.
///
/// Design: `traits-2-configuration-and-hashing.md` §"apply" (the "driver-owned
/// reusable batch" / in-place-swap note) and §"Pauli algebra traits" (the sum
/// "applies the one-row action pointwise and drains each term's phase delta to
/// its coefficient"). The bijectivity this relies on is machine-checked in
/// `lean/PPVM/Pauli/Symplectic.lean` (`*_bijective`).
pub trait RekeyBijective<K, C> {
    /// Re-key every term by moving it through `f` and re-inserting the result,
    /// reusing the backing allocation. `f` maps `(k, c) ↦ (φ(k), c')`.
    ///
    /// The caller guarantees `f` is **injective on keys** (a Clifford
    /// conjugation is a Pauli bijection, machine-checked as `*_bijective` in
    /// `lean/PPVM/Pauli/Symplectic.lean`), so re-keyed terms never collide.
    ///
    /// Implementations may therefore insert without aggregating. That is a real
    /// precondition, not a hint: violating it **drops** a term rather than summing
    /// it. The `HashMapStore` impl `debug_assert!`s that no insert displaces an
    /// existing entry, so a caller-side bug trips loudly in debug builds. The
    /// defensive `and_modify`/`or_insert` this replaced cost ~1.3× on the whole
    /// re-key path (it compiled out of line), which is too much to pay for
    /// guarding a case the type-level contract already excludes.
    fn rekey_bijective<F>(&mut self, f: F)
    where
        F: FnMut(K, C) -> (K, C);
}

/// An in-place, **whole-map** per-key sign flip — the fast path a *pure-sign*
/// Clifford conjugation (`X`/`Y`/`Z`) drives instead of the move-based
/// [`RekeyBijective`] re-key.
///
/// # Friction: `X`/`Y`/`Z` leave the word fixed, so a re-key is pure waste
///
/// A single-qubit `X`/`Y`/`Z` conjugation is a **pure sign**: `XPX = (−1)^z P`,
/// `YPY = (−1)^{x⊕z} P`, `ZPZ = (−1)^x P` — the Pauli **word is unchanged**, only
/// each term's coefficient picks up a `±1` read off that term's own X/Z bits. The
/// generic [`RekeyBijective`] path still rebuilds the whole `HashMap` (fresh
/// alloc, move every entry into a new bucket) even though every key maps to
/// itself; the old crate (`ppvm-pauli-sum::PauliSum::x`/`y`/`z`) instead scaled
/// each existing entry's coefficient in place (`scale` with a conditional
/// `*= -1.0`), touching **zero** allocations and never moving a key.
///
/// This trait restores that zero-allocation fast path: walk the stored entries in
/// place and multiply each coefficient by `f(&k) ∈ {−1, +1}`, computed from that
/// key's bits. It is deliberately a **whole-map** capability rather than a
/// per-slot `&mut (K, C)` / `&mut [C]` callback so a columnar backend (bits in one
/// column, coefficients in another) stays expressible — the impl owns the walk and
/// may read the key column and scale the coefficient column independently. This is
/// a `ppvm-pauli-sum-2`-local trait `impl`'d on the foreign `Vec`/`HashMap`
/// containers (legal: the trait is local), mirroring [`RekeyBijective`] and
/// [`StoreAlloc`].
///
/// Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits" (the
/// sum "applies the one-row action pointwise and drains each term's phase delta to
/// its coefficient"). The pure-sign conjugation signs are machine-checked in
/// `lean/PPVM/Pauli/Conjugation.lean` and the phased word's
/// `PhaseTrack::x_phase`/`y_phase`/`z_phase`
/// (`ppvm-phased-pauli-word-2/src/clifford.rs`); the word being fixed is the
/// identity symplectic bit map of `lean/PPVM/Pauli/Symplectic.lean`.
pub trait SignFlipByKey<K, C> {
    /// Multiply every term's coefficient in place by `f(&k) ∈ {−1, +1}`, keyed on
    /// that term's own key — no key movement, no reallocation.
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8;
}

/// An in-place, **whole-map** per-key coefficient scale by an arbitrary ring
/// factor — the fast path a *diagonal* unital Pauli channel ([`PauliError`]) drives
/// instead of a batch rebuild.
///
/// # Friction: noise scales by a real λ_P, not the `±1` [`SignFlipByKey`] allows
///
/// A unital single-qubit Pauli channel acts **diagonally** in the Pauli basis:
/// each term `P` is multiplied by a real transfer eigenvalue
/// `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q` that depends only on `P`'s Pauli at
/// the target qubit — no key moves, no branch, no rebuild (the old
/// `ppvm-pauli-sum::PauliSum::pauli_error`'s in-place `scale`). [`SignFlipByKey`]
/// is the `{−1, +1}` special case of exactly this walk, but a noise factor is a
/// general `C`, so it cannot go through `mul_sign(i8)`; this trait generalizes the
/// per-key walk to multiply by an arbitrary `f(&k) ∈ C`. Returning `None` leaves a
/// term untouched (the identity-at-qubit and lost-qubit cases), so the walk never
/// fabricates a `C::one()` for a no-op slot.
///
/// Like its siblings it is deliberately a **whole-map** capability rather than a
/// per-slot `&mut (K, C)` callback so a columnar backend stays expressible, and it
/// is a `ppvm-pauli-sum-2`-local trait `impl`'d on the foreign `Vec`/`HashMap`
/// containers (legal: the trait is local).
///
/// Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits"
/// (`PauliError`) and §"Pauli algebra traits". The eigenvalue this scales by is
/// machine-checked in `lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`,
/// `pauli_channel_eigenvalue_omega`).
pub trait ScaleByKey<K, C> {
    /// Multiply every term's coefficient in place by `f(&k)`, keyed on that term's
    /// own key — no key movement, no reallocation. A `None` result leaves the
    /// term unchanged (an exact no-op, not a multiply by one).
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> Option<C>;
}

/// An in-place, **fused branching** re-key — the fast path a single-qubit
/// rotation drives instead of the general [`apply`](crate::Sum::apply) batch
/// round-trip.
///
/// # Friction: `apply` re-hashes the whole 2N-term fan-out; a rotation's diagonal never moves
///
/// A single-qubit rotation `exp(−i·θ/2·G)` conjugates each term `(P, c)` to
/// `(P, c·cosθ)` **plus**, when `G` anticommutes with `P`, a genuinely-new branch
/// `(iGP, c·sinθ·ε)`. The **diagonal** `P` keeps its key — only its coefficient
/// is scaled by `cosθ` — so its cached structural hash stays valid and it should
/// never leave its bucket. Routing a rotation through `apply`, though, clones
/// every key on the read side ([`Support::iter`]), `reset`s the whole map, then
/// re-accumulates **all** ~`2N` produced terms — re-hashing and re-inserting even
/// the `N` untouched diagonals. Against the old crate's fused single pass
/// (`ppvm-pauli-sum::sum::rot1` over `PauliSum::map_insert`), which mutates each
/// diagonal coefficient in place and hashes only the ≤`N` branch terms, the batch
/// path benchmarked ~2.4× slower on an `rx` over a ~1000-term sum.
///
/// This trait restores that fused single pass: walk the stored entries in place,
/// letting `f(&k, &mut c)` scale each diagonal coefficient **in place** (no key
/// clone, no re-hash) and return the optional branch term; buffer the branch
/// terms (cheap, un-hashed pushes) and merge only those back through the map's
/// hash-join, aggregating any collision (a branch may land on another term). It
/// is deliberately a **whole-map** capability rather than a per-slot callback so a
/// columnar backend stays expressible, and it is a `ppvm-pauli-sum-2`-local trait
/// `impl`'d on the foreign `Vec`/`HashMap` containers (legal: the trait is local),
/// mirroring [`SignFlipByKey`], [`ScaleByKey`], and [`RekeyBijective`].
///
/// A zero branch coefficient (an identity rotation `θ = 0` has `sinθ = 0`) is
/// never inserted, so `R₀` adds no spurious key. Beyond that no whole-map `reduce`
/// scan runs: a generic rotation's *collision* cancellations are measure-zero, so
/// — like the old crate's `map_insert`, which leaves any residue for a later
/// truncation — the fast path skips the retain [`apply`](crate::Sum::apply) would
/// run. A physical near-zero is dropped by the policy's truncation the caller runs
/// afterward.
///
/// Design: `traits-2-configuration-and-hashing.md` §"apply" and §"Every gate is a
/// producer feeding `accumulate`" (the fused in-place alternative to the batch
/// round-trip), and §"Behavioral traits" (`RotationOne`). The branch
/// `c·P → cos·c·P + sin·c·(iGP)` — a genuinely-new key, norm-preserving and
/// angle-additive — is machine-checked in `lean/PPVM/Instantiations/Rotation.lean`
/// (`anticommute_new_key`, `rot_norm_sq`, `rot_rot`).
pub trait RotateInPlace<K, C> {
    /// Walk every term `(k, c)` in place: `f(&k, &mut c)` scales the diagonal
    /// coefficient in place (the `cosθ` factor) and returns `Some((k′, c′))` for
    /// the branch term to merge, or `None` when the term commutes. Colliding
    /// branch keys are aggregated; a zero branch is skipped (identity rotation) and
    /// no whole-map `reduce` scan runs (see the trait docs — collision
    /// cancellations are measure-zero and left for the policy's truncation).
    fn rotate_in_place<F>(&mut self, f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>;
}

impl<K, C> ScaleByKey<K, C> for Vec<(K, C)>
where
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> Option<C>,
    {
        // Drop any term whose scaled coefficient becomes exactly zero (e.g. a
        // zero channel eigenvalue, `pauli_error(q, [0.0, 0.25, 0.25])` → λ_X = 0),
        // preserving the reduced-canonical-form invariant (`Accumulate::reduce`:
        // no zero coefficient stays in the support).
        self.retain_mut(|(k, c)| match f(k) {
            None => true,
            Some(factor) => {
                *c *= factor;
                !c.is_zero()
            }
        });
    }
}

impl<K, C> ScaleByKey<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> Option<C>,
    {
        // Walk the existing buckets in place: read each key, multiply its
        // coefficient by the per-key factor. No bucket is moved and the backing
        // allocation is untouched — the old crate's in-place `scale`, restored.
        // A term scaled to exactly zero (a zero channel eigenvalue) is dropped,
        // preserving the reduced-canonical-form invariant (`Accumulate::reduce`).
        self.retain(|k, c| match f(k) {
            None => true,
            Some(factor) => {
                *c *= factor;
                !c.is_zero()
            }
        });
    }
}

impl<K, C> RotateInPlace<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>,
    {
        // Pass 1: scale each diagonal coefficient in place, buffering the
        // anticommuting branch terms. The diagonal keeps its slot — no clone.
        let mut branches: Vec<(K, C)> = Vec::with_capacity(self.len());
        for (k, c) in self.iter_mut() {
            if let Some(term) = f(k, c) {
                branches.push(term);
            }
        }
        // Pass 2: merge only the branch terms back (linear-scan hash-join, the
        // small-support cost model this backend targets), aggregating collisions.
        // A zero branch (e.g. `sinθ = 0` at `θ = 0`, the identity rotation) is
        // never inserted, so an identity rotation adds no spurious key. Beyond
        // that no whole-map `reduce` scan runs — a generic rotation's collision
        // cancellations are measure-zero and left for the policy's truncation,
        // matching the old crate's `map_insert`.
        for (nk, nc) in branches {
            if nc.is_zero() {
                continue;
            }
            // `as_slice` disambiguates from the now-in-scope `Support::iter`
            // (which would yield owned pairs) — we want the borrowing slice iter.
            if let Some(pos) = self.as_slice().iter().position(|(ek, _)| *ek == nk) {
                self[pos].1 += nc;
            } else {
                self.push((nk, nc));
            }
        }
    }
}

impl<K, C> SignFlipByKey<K, C> for Vec<(K, C)>
where
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        for (k, c) in self.iter_mut() {
            *c = c.mul_sign(f(k));
        }
    }
}

impl<K, C> SignFlipByKey<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        // Walk the existing buckets in place: read each key, scale its coefficient
        // by the `±1` that key's bits demand. No bucket is moved and the backing
        // allocation is untouched — the old crate's in-place `scale`, restored.
        for (k, c) in self.iter_mut() {
            *c = c.mul_sign(f(k));
        }
    }
}

impl<K, C> RekeyBijective<K, C> for Vec<(K, C)>
where
    K: Eq + Clone,
    C: Coefficient,
{
    #[inline]
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C),
    {
        // Coordinate-list backend (small support): move every term out with
        // `mem::take` (no clones) and rebuild by re-keying each in turn. A
        // bijection maps distinct keys to distinct keys, so a positional rebuild
        // needs no aggregation.
        let taken = std::mem::take(self);
        *self = taken.into_iter().map(|(k, c)| f(k, c)).collect();
    }
}

impl<K, C> StoreAlloc for Vec<(K, C)> {
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        Vec::with_capacity(cap)
    }

    #[inline]
    fn reset(&mut self) {
        self.clear();
    }
}

// ===========================================================================
// HashMapStore — the double-buffered hash-join backend (primary + aux + scratch).
//
// The graded algebra (`Support`/`Accumulate`/`Scale`/`Pair`/`Retain`) and the
// pure-diagonal capabilities (`ScaleByKey`/`SignFlipByKey`) delegate to the
// `primary` map's existing impls; only the buffer-using fast paths
// (`StoreAlloc`/`RekeyBijective`/`RotateInPlace`) are bespoke.
// ===========================================================================

impl<K, C> Support for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    type Key = K;
    type Coeff = C;

    #[inline]
    fn len(&self) -> usize {
        Support::len(&self.primary)
    }

    #[inline]
    fn get(&self, key: &K) -> Option<C> {
        Support::get(&self.primary, key)
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = (K, C)> {
        Support::iter(&self.primary)
    }
}

impl<K, C> Accumulate for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn accumulate_batch(&mut self, terms: &TermBatch<K, C>) {
        Accumulate::accumulate_batch(&mut self.primary, terms);
    }

    #[inline]
    fn reduce(&mut self) {
        Accumulate::reduce(&mut self.primary);
    }
}

impl<K, C> Scale for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale(&mut self, s: &C) {
        Scale::scale(&mut self.primary, s);
    }
}

impl<K, C> Pair for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn probe_batch(&self, keys: &KeyBatch<K>, out: &mut [Option<C>]) {
        Pair::probe_batch(&self.primary, keys, out);
    }

    #[inline]
    fn overlap(&self, other: &Self) -> C {
        Pair::overlap(&self.primary, &other.primary)
    }

    #[inline]
    fn hermitian_overlap(&self, other: &Self) -> C
    where
        C: Conjugate,
    {
        Pair::hermitian_overlap(&self.primary, &other.primary)
    }
}

impl<K, C> Retain<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn retain(&mut self, keep: impl Fn(&K, &C) -> bool) {
        Retain::retain(&mut self.primary, keep);
    }
}

impl<K, C> StoreAlloc for HashMapStore<K, C> {
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        Self {
            primary: HashMap::with_capacity_and_hasher(cap, IdentityBuildHasher),
            // Size the aux to the same hint so the first re-key never resizes.
            aux: HashMap::with_capacity_and_hasher(cap, IdentityBuildHasher),
            scratch: Vec::with_capacity(cap),
        }
    }

    /// Reset to empty support, keeping the primary's allocation (the aux/scratch
    /// are already empty between ops).
    #[inline]
    fn reset(&mut self) {
        self.primary.clear();
    }
}

impl<K, C> ScaleByKey<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn scale_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> Option<C>,
    {
        ScaleByKey::scale_by_key(&mut self.primary, f);
    }
}

impl<K, C> SignFlipByKey<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn sign_flip_by_key<F>(&mut self, f: F)
    where
        F: Fn(&K) -> i8,
    {
        SignFlipByKey::sign_flip_by_key(&mut self.primary, f);
    }
}

impl<K, C> RekeyBijective<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C),
    {
        // The old crate's `map_add`, with the aux living in the store: clear the
        // persistent aux, **move** every term through the re-key into it
        // (`drain` yields owned keys — zero key clones), then swap. `drain` empties
        // the primary but keeps its allocation, and the swap hands it to `aux` — so
        // the two allocations ping-pong across gates and are never freed (the
        // double-buffer, recovered).
        //
        // The merge is a plain `insert`, not `entry(..).and_modify(..).or_insert(..)`.
        // `f` is required to be **injective on keys** (this trait's contract; for a
        // Clifford conjugation it is the symplectic bijection machine-checked as
        // `*_bijective` in `lean/PPVM/Pauli/Symplectic.lean`), so two terms can never
        // land on the same re-keyed slot and there is nothing to accumulate — the
        // `entry` form was guarding a collision that cannot occur.
        //
        // That guard was not free. On this crate's monomorphization the `entry`
        // chain compiled **out of line** (a standalone `hashbrown::rustc_entry`
        // frame took 41% of this loop's samples, where the old crate's identical
        // source construct inlines fully into `map_add_assign`). It cost ~1.44x
        // against old *and* made the whole re-key bimodal across processes: with
        // `entry` ~40% of processes landed in a ~1.5x-slower mode, with `insert`
        // none of 34 did. See `docs/log.md` (`ps2.cnot.rekey.perf`,
        // `ps2.rekey.bimodal`).
        self.aux.clear();
        self.aux.reserve(self.primary.len());
        for (k, c) in self.primary.drain() {
            let (nk, nc) = f(k, c);
            let displaced = self.aux.insert(nk, nc);
            debug_assert!(
                displaced.is_none(),
                "RekeyBijective requires an injective re-key: two distinct keys \
                 collided, so a term was silently dropped"
            );
        }
        std::mem::swap(&mut self.primary, &mut self.aux);
    }
}

impl<K, C> RotateInPlace<K, C> for HashMapStore<K, C>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn rotate_in_place<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut C) -> Option<(K, C)>,
    {
        // Pass 1: scale each diagonal coefficient in place (the key is unchanged,
        // so its cached hash stays valid and the entry never moves), buffering the
        // anticommuting branch terms in the persistent `scratch` — the old crate's
        // reused `map_insert` scratch, now living in the store (no per-gate Vec
        // allocation).
        self.scratch.clear();
        for (k, c) in self.primary.iter_mut() {
            if let Some(term) = f(k, c) {
                // Pre-warm the fresh branch key's structural digest HERE, in pass 1,
                // before it is buffered. `with_bits_toggled` builds the key with an
                // empty (`HASH_UNCACHED`) cache, so without this the 3-round finalize
                // fold would fire lazily inside pass-2's `entry()` — *on the bucket-
                // index critical path*, where the mul-chain latency stalls the
                // dependent bucket load with nothing to hide it. Computing the (same)
                // digest here instead lets it overlap with the walk's other work
                // (the next term's diagonal scale + branch build), so pass-2 probes a
                // cache hit. Semantic no-op (identical digest); measured ~8% off `rx`
                // on a ~1000-term sum, matching the old crate's eager-hash `map_insert`
                // (which likewise hashes in its first pass).
                let _ = term.0.key_hash();
                self.scratch.push(term);
            }
        }
        // Pass 2: merge only the branch terms through the hash-join, aggregating
        // any collision with an existing (already-scaled) diagonal or a sibling
        // branch — one hash pass over ≤N terms, not the whole 2N fan-out.
        //
        // A zero branch (e.g. `sinθ = 0` at `θ = 0`, the identity rotation) is
        // never inserted, so an identity rotation adds no spurious key and needs no
        // whole-map cleanup. Beyond that no `reduce` scan runs: a generic
        // rotation's collision cancellations are measure-zero, so — like the old
        // crate's `map_insert` (`ppvm-pauli-sum::sum::rot1`), which leaves any
        // residue for a later truncation — we skip the whole-map retain the general
        // `apply` would run. A physical near-zero is dropped by the policy's
        // truncation (or the caller's magnitude floor), not here.
        //
        // `drain(..)` empties `scratch` but keeps its allocation for the next gate.
        for (nk, nc) in self.scratch.drain(..) {
            if nc.is_zero() {
                continue;
            }
            self.primary
                .entry(nk)
                .and_modify(|e| *e += nc.clone())
                .or_insert(nc);
        }
    }
}
