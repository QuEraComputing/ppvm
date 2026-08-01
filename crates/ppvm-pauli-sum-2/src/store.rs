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

use ppvm_traits_2::{Coefficient, IdentityBuildHasher, Indexable};

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
    /// conjugation is a Pauli bijection), so re-keyed terms never collide. The
    /// merge still aggregates defensively — matching the old `map_add`'s
    /// `and_modify`/`or_insert` — so a caller-side bug degrades a lost term into
    /// a summed one rather than silently dropping coefficient mass.
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

impl<K, C> RekeyBijective<K, C> for HashMap<K, C, IdentityBuildHasher>
where
    K: Indexable,
    C: Coefficient,
{
    #[inline]
    fn rekey_bijective<F>(&mut self, mut f: F)
    where
        F: FnMut(K, C) -> (K, C),
    {
        // Swap the whole map out for a fresh one sized to the same support, then
        // **consume** the old map — moving each `(k, c)` out exactly once — re-key
        // it with zero key clones, and insert by value into the fresh map. This
        // is the old crate's two-map swap (`map_add`), minus the per-term key
        // clone it still paid, and it moves each key only twice (old bucket →
        // register → new bucket) versus the three moves a drain-into-scratch
        // staging would cost. Sizing the fresh map to the current support — which
        // a bijection preserves — means the re-insert never resizes.
        let fresh = HashMap::with_capacity_and_hasher(self.len(), IdentityBuildHasher);
        let old = std::mem::replace(self, fresh);
        for (k, c) in old {
            let (nk, nc) = f(k, c);
            self.entry(nk)
                .and_modify(|e| *e += nc.clone())
                .or_insert(nc);
        }
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

impl<K, C> StoreAlloc for HashMap<K, C, IdentityBuildHasher> {
    #[inline]
    fn with_capacity(cap: usize) -> Self {
        HashMap::with_capacity_and_hasher(cap, IdentityBuildHasher)
    }

    #[inline]
    fn reset(&mut self) {
        self.clear();
    }
}
