// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::hash::BuildHasher;

use crate::traits::{Coefficient, PauliStorage, PauliWordTrait};

/// Minimal interface for any "associative coefficient map" backing a
/// `PauliSum` — construction, length, and clear.
///
/// Implementations exist for `HashMap`, `IndexMap`, and `DashMap`; pick
/// one via a [`Config`](crate::config::Config).
pub trait ACMapBase {
    /// Construct an empty map with at least `capacity` slots reserved.
    fn with_capacity(capacity: usize) -> Self;
    /// Number of stored `(key, value)` pairs.
    fn len(&self) -> usize;
    /// `true` if the map is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Remove every entry, retaining allocated capacity.
    fn clear(&mut self);
}

/// Borrowing iteration over an `ACMap`. Lives in its own trait so that
/// implementations may pick their own item / iterator types.
pub trait ACMapIter<'a> {
    /// Yielded item type.
    type Item;
    /// Iterator type.
    type Iter: Iterator<Item = Self::Item>;
    /// Iterate over `(key, value)` pairs (or their representation).
    fn iter(&'a self) -> Self::Iter;
}

/// `+=` semantics for an `ACMap`: insert a new entry or accumulate into
/// the existing one with the same key.
pub trait ACMapAddAssign<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// Add `value` into the entry at `key`, creating the entry if absent.
    fn add_assign(&mut self, key: W, value: V);
    /// For every entry, compute `f(key, value)` and add the result into
    /// `dest` with [`add_assign`](Self::add_assign).
    fn map_add_assign<F>(&self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &V) -> (W, V) + Sync + Send;
}

/// Scalar `*=` semantics for an `ACMap`: multiply every value by a
/// constant.
pub trait ACMapMulAssign<V: Coefficient, H: BuildHasher + Clone + Default> {
    /// Scale every value in place.
    fn mul_assign(&mut self, value: V);
}

/// In-place modify + insert pattern used to express branching gates
/// (where one input entry can produce zero, one, or several output
/// entries).
pub trait ACMapInsert<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// Modify each existing entry in place; if `f` returns `Some((k', v'))`,
    /// append that new entry into `dest`, a plain `Vec` (a push per entry,
    /// no hashing). This is the primary hot-path entry point: the caller
    /// merges the buffer into the destination map afterwards, avoiding a
    /// second hashmap probe per produced entry.
    fn map_insert_vec<F>(&mut self, dest: &mut Vec<(W, V)>, f: F)
    where
        F: Fn(&W, &mut V) -> Option<(W, V)> + Sync + Send;

    /// Like [`map_insert_vec`](Self::map_insert_vec) but `f` may return a
    /// `Vec` of new entries per existing entry, inserted directly into the
    /// destination map `dest`.
    fn map_insert_multiple<F>(&mut self, dest: &mut Self, f: F)
    where
        F: Fn(&W, &mut V) -> Option<Vec<(W, V)>> + Sync + Send;
}

/// Membership queries — `(key, value)` exact match or with a custom
/// predicate on the value.
pub trait ACMapContains<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// `true` if an entry with this exact `(key, value)` is present.
    fn contains(&self, key: &W, value: &V) -> bool {
        self.contains_with(key, |v| v == value)
    }
    /// `true` if an entry for `key` exists whose value satisfies `f`.
    fn contains_with<F>(&self, key: &W, f: F) -> bool
    where
        F: Fn(&V) -> bool;
}

/// Merge two maps with accumulation: drain `dest` into `self`,
/// summing values that share a key.
pub trait ACMapConsume {
    /// Drain `dest` into `self`, accumulating values on key collision.
    fn consume(&mut self, dest: &mut Self);
}

/// In-place per-entry transformation of values.
pub trait ACMapScale<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// Apply `f(key, value)` to every entry; only the value is mutable.
    fn scale<F>(&mut self, f: F)
    where
        F: Fn(&W, &mut V) + Sync + Send;
}

/// Drop entries that don't satisfy a predicate — used by truncation
/// strategies.
pub trait ACMapRetain<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>
{
    /// Keep only entries for which `f(key, value)` returns `true`.
    fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&W, &V) -> bool;
}

/// Aggregate trait combining every operation a backing map must support
/// to be usable as the storage for a `PauliSum`.
///
/// You don't normally implement `ACMap` directly: the blanket impl below
/// covers any type that implements all the constituent traits.
///
/// # Map-backend implementor's guide
///
/// A new backing map (a GPU-resident map, a sharded map, anything beyond
/// the existing `HashMap` / `IndexMap` / `AHashMap` / `DashMap` backends)
/// must implement exactly these traits — nothing more:
///
/// * [`ACMapBase`] — construction, length, clear.
/// * [`ACMapIter`] — borrowing iteration over `(key, value)`.
/// * [`ACMapAddAssign`] — `+=` insert-or-accumulate.
/// * [`ACMapMulAssign`] — scalar `*=`.
/// * [`ACMapInsert`] — just `map_insert_vec` and `map_insert_multiple`
///   (there is no third variant to implement).
/// * [`ACMapContains`] — membership queries.
/// * [`ACMapScale`] — in-place per-entry transform.
/// * [`ACMapRetain`] — predicate-based entry removal.
/// * [`ACMapConsume`] — drain-and-accumulate merge of two maps.
///
/// Of these, `map_insert_vec`, `map_insert_multiple`, `scale`,
/// `mul_assign`, `retain`, and `consume` are whole-collection batch entry
/// points: a concurrent or GPU backend can dispatch each of them as one
/// parallel pass over every entry. `add_assign` is the exception — it
/// operates on a single key at a time, so a thread-pool backend typically
/// only parallelizes its batch sibling, `map_add_assign`. See the
/// `DashMap` impl in `crates/ppvm-traits/src/map/dashmap.rs` for a working
/// concurrent precedent: it implements every batch method with rayon
/// `par_iter`.
pub trait ACMap<
    S: PauliStorage,
    V: Coefficient,
    H: BuildHasher + Clone + Default,
    W: PauliWordTrait,
>:
    Clone
    + ACMapBase
    + ACMapAddAssign<S, V, H, W>
    + ACMapMulAssign<V, H>
    + ACMapInsert<S, V, H, W>
    + ACMapContains<S, V, H, W>
    + ACMapScale<S, V, H, W>
    + ACMapRetain<S, V, H, W>
    + ACMapConsume
{
}

impl<T, Storage, Coeff, Hasher, Word> ACMap<Storage, Coeff, Hasher, Word> for T
where
    Storage: PauliStorage,
    Coeff: Coefficient,
    Hasher: BuildHasher + Clone + Default,
    Word: PauliWordTrait,
    T: Clone
        + ACMapBase
        + ACMapAddAssign<Storage, Coeff, Hasher, Word>
        + ACMapMulAssign<Coeff, Hasher>
        + ACMapInsert<Storage, Coeff, Hasher, Word>
        + ACMapScale<Storage, Coeff, Hasher, Word>
        + ACMapContains<Storage, Coeff, Hasher, Word>
        + ACMapRetain<Storage, Coeff, Hasher, Word>
        + ACMapConsume,
{
}
