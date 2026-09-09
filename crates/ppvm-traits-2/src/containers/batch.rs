// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::containers::Indexable;
use crate::loss::LossState;
use crate::word::PauliBits;

/// A key that can be laid out as a structure-of-arrays column. Separate from
/// [`Indexable`] so the minimal hashing contract is unchanged: a batched key is
/// both `Indexable` (a valid map key) and `Columnar` (has a column layout).
///
/// Design: §"The batch contract".
pub trait Columnar: Indexable {
    /// The concrete structure-of-arrays column for this key type.
    type Column: KeyColumn<Key = Self>;
}

/// A structure-of-arrays column of keys, owned by the concrete key type (only it
/// knows its planes). Operates plane by plane, never scalar on the hot path.
///
/// Design: §"The batch contract".
pub trait KeyColumn: Default + Clone {
    /// The key type this column stores.
    type Key: Columnar;

    /// Number of keys currently in the column.
    fn len(&self) -> usize;

    /// Number of keys the existing plane allocations can hold without growing.
    ///
    /// Live stores use this when cloning persistent workspaces: cloning only the
    /// populated rows would silently discard a caller's capacity hint.
    fn capacity(&self) -> usize;

    /// Whether the column is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// A column pre-sized for `n` keys.
    fn with_capacity(n: usize) -> Self;

    /// Append one produced key; the column keeps each plane contiguous.
    fn push(&mut self, key: Self::Key);

    /// Reserve room for `additional` more keys, keeping the ones already stored.
    ///
    /// The column spelling of `Vec::reserve`, and the counterpart of
    /// [`with_capacity`](Self::with_capacity) for a column that is a *live
    /// support* rather than a throwaway batch: `ColumnStore`'s branch-merge pass
    /// knows its worst-case append count up front (the scratch length), and
    /// pre-sizing from it collapses a doubling chain of plane reallocations —
    /// plus the parallel bucket-table `reindex`es — into one.
    ///
    /// **Default: a no-op.** Pre-sizing is a pure optimization, so a column that
    /// cannot express it (a device-backed or fixed-extent one) stays legal and
    /// simply reallocates on `push`.
    #[inline]
    fn reserve(&mut self, additional: usize) {
        let _ = additional;
    }

    /// Bulk structural hash of the whole column into a parallel hash column.
    /// `out[i]` must equal the `i`-th key's [`Indexable::key_hash`].
    fn hash_into(&self, out: &mut [u64]);

    /// Join confirm: compare element `i` against a build-side key after a hash
    /// or tag match, without materializing the whole element.
    fn key_eq(&self, i: usize, other: &Self::Key) -> bool;

    /// Select or permute elements into a new column (radix partitioning,
    /// compaction, device staging). Operates plane by plane, never scalar.
    fn gather(&self, indices: &[u32]) -> Self;

    /// Scalar materialization of one element — a naive backend's fallback,
    /// never the hot path.
    fn get(&self, i: usize) -> Self::Key;

    /// Read one row's X bit without requiring scalar materialization when the
    /// concrete column can address its packed plane directly.
    #[inline]
    fn x_bit(&self, row: usize, qubit: usize) -> bool
    where
        Self::Key: PauliBits,
    {
        self.get(row).x_bit(qubit)
    }

    /// Read one row's Z bit directly when supported.
    #[inline]
    fn z_bit(&self, row: usize, qubit: usize) -> bool
    where
        Self::Key: PauliBits,
    {
        self.get(row).z_bit(qubit)
    }

    /// Read one row's loss bit directly when supported.
    #[inline]
    fn is_lost(&self, row: usize, qubit: usize) -> bool
    where
        Self::Key: LossState,
    {
        self.get(row).is_lost(qubit)
    }

    /// Materialize one row while toggling selected bits. Packed columns can
    /// build the branch key directly from their planes.
    #[inline]
    fn toggled_bits(&self, row: usize, qubit: usize, toggle_x: bool, toggle_z: bool) -> Self::Key
    where
        Self::Key: PauliBits,
    {
        self.get(row).toggled_bits(qubit, toggle_x, toggle_z)
    }

    /// Materialize one row while toggling bits at two sites.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn toggled_bits2(
        &self,
        row: usize,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self::Key
    where
        Self::Key: PauliBits,
    {
        self.get(row)
            .toggled_bits2(i, toggle_x_i, toggle_z_i, j, toggle_x_j, toggle_z_j)
    }

    /// Reset to an empty column, **keeping the backing plane allocations**.
    ///
    /// # Friction: a column that is a *store* needs in-place mutation, and the
    /// batch-only surface has none
    ///
    /// [`gather`](Self::gather) allocates a fresh column, and the design's
    /// batch contract needs nothing more: a `TermBatch`'s column is built by
    /// `push` and thrown away. The `ColumnStore` backend
    /// (implementation-plan Phase 6) makes the *same* column type the live
    /// support, and every one of its buffer-reusing fast paths — the old
    /// crate's `map_add` clear→write→swap (architecture feature 1), the retain
    /// compaction, the in-place Clifford re-key — is defined by mutating a
    /// column it already owns. Expressed through `gather` alone each of those
    /// allocates a whole new key column **per gate**, which is exactly the
    /// per-gate allocation churn the double-buffer exists to remove.
    ///
    /// So this trio ([`clear`](Self::clear), [`set`](Self::set),
    /// [`truncate`](Self::truncate)) is the minimal in-place surface: they are
    /// the column spellings of `Vec::clear`/`IndexMut`/`Vec::truncate`, they
    /// stay plane-oriented (a SIMD/GPU column implements them as plane writes),
    /// and they expose no `&mut Key` — so design rule 4 of §"Backends are
    /// containers" ("no signature exposes `&mut (W, C)` or `&mut [C]`") is
    /// untouched and the AoS layout still cannot leak.
    fn clear(&mut self);

    /// Overwrite element `i` in place, keeping every other element and the
    /// backing allocation. Panics (or debug-panics) if `i >= len()`.
    ///
    /// The write side of the in-place re-key: a Clifford conjugation is a
    /// bijection, so a columnar backend rewrites the key planes at each slot and
    /// leaves the parallel coefficient column completely untouched.
    fn set(&mut self, i: usize, key: Self::Key);

    /// Shorten to `len` elements, keeping the backing allocation. A no-op if
    /// `len >= self.len()`. The tail step of a stable retain compaction
    /// (`set` the survivors down, then cut).
    fn truncate(&mut self, len: usize);

    /// Remove row `i` by moving the final row into its slot.
    fn swap_remove(&mut self, i: usize) -> Self::Key {
        let len = self.len();
        let removed = self.get(i);
        if i + 1 != len {
            self.set(i, self.get(len - 1));
        }
        self.truncate(len - 1);
        removed
    }
}

/// Keys plus their precomputed structural hashes, in parallel columns. The
/// probe side of the join; it carries no coefficients.
///
/// See the module-level friction note: the key column is a scalar `Vec<W>`
/// fallback rather than the design's `W::Column`, so the batch is expressible
/// for any `W: Eq + Clone` (not only `Columnar` keys).
///
/// Design: §"The batch contract".
#[derive(Debug, Clone)]
pub struct KeyBatch<W> {
    keys: Vec<W>,
    hashes: Vec<u64>,
    hashes_valid: bool,
}

impl<W> Default for KeyBatch<W> {
    fn default() -> Self {
        Self {
            keys: Vec::new(),
            hashes: Vec::new(),
            hashes_valid: false,
        }
    }
}

impl<W> KeyBatch<W> {
    /// An empty key batch.
    pub fn new() -> Self {
        Self::default()
    }

    /// A key batch pre-sized for `n` keys.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            keys: Vec::with_capacity(n),
            hashes: Vec::with_capacity(n),
            hashes_valid: false,
        }
    }

    /// Number of keys in the batch.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Number of keys the existing columns can hold without growing.
    pub fn capacity(&self) -> usize {
        self.keys.capacity().min(self.hashes.capacity())
    }

    /// Whether the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// The key column.
    pub fn keys(&self) -> &[W] {
        &self.keys
    }

    /// The complete parallel hash column, or `None` if it has not been filled
    /// since the last mutation. A filled empty batch returns `Some(&[])`.
    pub fn hashes(&self) -> Option<&[u64]> {
        self.hashes_valid.then_some(self.hashes.as_slice())
    }

    /// Iterate keys in insertion order.
    pub fn iter(&self) -> impl Iterator<Item = &W> {
        self.keys.iter()
    }

    /// Clear both columns without releasing capacity (for buffer reuse).
    pub fn clear(&mut self) {
        self.hashes_valid = false;
        self.keys.clear();
        self.hashes.clear();
    }

    /// Append a key and invalidate cached hashes without releasing capacity.
    /// Call [`KeyBatch::fill_hashes`] to make the full hash column available again.
    pub fn push(&mut self, key: W) {
        self.hashes_valid = false;
        self.hashes.clear();
        self.keys.push(key);
    }
}

impl<W: Indexable> KeyBatch<W> {
    /// Fill the parallel hash column from each key's [`Indexable::key_hash`], so
    /// `hashes().unwrap()[i] == keys()[i].key_hash()`.
    /// The cache becomes available only after all keys have been hashed.
    pub fn fill_hashes(&mut self) {
        self.hashes_valid = false;
        self.hashes.clear();
        self.hashes
            .extend(self.keys.iter().map(Indexable::key_hash));
        self.hashes_valid = true;
    }
}

/// A [`KeyBatch`] with the coefficient column attached: the produced terms
/// awaiting merge. Coefficients are a separate column, touched only when a
/// probe resolves to an aggregate.
///
/// Design: §"The batch contract".
#[derive(Debug, Clone)]
pub struct TermBatch<W, C> {
    keys: KeyBatch<W>,
    coeffs: Vec<C>,
}

impl<W, C> Default for TermBatch<W, C> {
    fn default() -> Self {
        Self {
            keys: KeyBatch::new(),
            coeffs: Vec::new(),
        }
    }
}

impl<W, C> TermBatch<W, C> {
    /// An empty term batch.
    pub fn new() -> Self {
        Self::default()
    }

    /// A term batch pre-sized for `n` terms.
    pub fn with_capacity(n: usize) -> Self {
        Self {
            keys: KeyBatch::with_capacity(n),
            coeffs: Vec::with_capacity(n),
        }
    }

    /// Number of terms in the batch.
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Number of terms all three columns can hold without growing.
    pub fn capacity(&self) -> usize {
        self.keys.capacity().min(self.coeffs.capacity())
    }

    /// Whether the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The probe-side key batch.
    pub fn keys(&self) -> &KeyBatch<W> {
        &self.keys
    }

    /// The coefficient column.
    pub fn coeffs(&self) -> &[C] {
        &self.coeffs
    }

    /// Iterate `(key, coeff)` pairs — the read side an `accumulate_batch` merge
    /// loop consumes. Synthesizes the pairs from the two columns; the layout
    /// stays structure-of-arrays.
    pub fn iter(&self) -> impl Iterator<Item = (&W, &C)> {
        self.keys.keys().iter().zip(self.coeffs.iter())
    }

    /// Clear all columns without releasing capacity (for buffer reuse).
    pub fn clear(&mut self) {
        self.keys.clear();
        self.coeffs.clear();
    }
}

/// The append side of a term batch: a producer pushes `(key, coeff)` terms into
/// a sink, filling the key and coefficient columns. A naive sink collects into a
/// scalar `Vec`; a columnar sink appends into planes.
///
/// Design: §"Every gate is a producer feeding `accumulate`".
pub trait TermSink<K, C> {
    /// Append one produced term.
    fn push(&mut self, key: K, coeff: C);
}

impl<W, C> TermSink<W, C> for TermBatch<W, C> {
    #[inline]
    fn push(&mut self, key: W, coeff: C) {
        self.keys.push(key);
        self.coeffs.push(coeff);
    }
}

/// A monomorphized, inlinable term producer — never `dyn`, since this is the
/// hot loop and the abstraction must compile to nothing.
///
/// # `Send + Sync` is part of the contract (architecture feature 12)
///
/// A producer is *read-only* over its own state (`produce` takes `&self`), so a
/// storage backend is free to split the produce walk across threads — which is
/// the whole point of keeping the backend a **configuration choice**: the old
/// crate bounded every driver closure `F: Fn(..) + Sync + Send`
/// (`ppvm-traits/src/map/hashmap.rs`, `map_add_assign`/`map_insert*`/`scale`) so
/// that a concurrent map (it shipped a `DashMap`-backed config benchmarked beside
/// the `HashMap`/`IndexMap` ones) was a **backend swap, not an engine rewrite**.
/// Requiring it here — and on the `ppvm-pauli-sum-2` in-place walk closures
/// (`ScaleByKey`/`SignFlipByKey`/`RekeyBijective`/`RotateInPlace`/
/// `BranchInPlace`) — keeps that door open: widening the bound later would mean
/// touching every trait signature *and* every impl, i.e. exactly the coupling the
/// feature exists to prevent. Every real producer is a closure over gate indices
/// and ring elements, and `Coefficient` is already `Send + Sync`, so the bound
/// costs nothing today.
pub trait TermProducer<K, C>: Send + Sync {
    /// Push the produced terms for one existing `(key, coeff)` into the sink.
    fn produce<S: TermSink<K, C>>(&self, key: &K, coeff: &C, sink: &mut S);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, PartialEq, Eq)]
    struct Key(u64);

    impl std::hash::Hash for Key {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            state.write_u64(self.key_hash());
        }
    }

    impl Indexable for Key {
        fn key_hash(&self) -> u64 {
            self.0.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        }
    }

    #[test]
    fn push_invalidates_hashes_and_refill_covers_every_key() {
        let mut batch = KeyBatch::with_capacity(4);
        batch.push(Key(1));
        assert_eq!(batch.hashes(), None);
        batch.fill_hashes();
        assert_eq!(batch.hashes(), Some([Key(1).key_hash()].as_slice()));
        let capacity = batch.capacity();
        let snapshot = batch.clone();

        batch.push(Key(2));
        assert_eq!(batch.hashes(), None);
        assert_eq!(batch.capacity(), capacity);
        assert_eq!(snapshot.hashes(), Some([Key(1).key_hash()].as_slice()));

        batch.fill_hashes();
        assert_eq!(
            batch.hashes(),
            Some([Key(1).key_hash(), Key(2).key_hash()].as_slice()),
        );
    }

    #[test]
    fn empty_and_cleared_batches_have_explicit_cache_state() {
        let mut batch = KeyBatch::<Key>::new();
        assert_eq!(batch.hashes(), None);
        batch.fill_hashes();
        assert_eq!(batch.hashes(), Some([].as_slice()));

        batch.push(Key(3));
        batch.fill_hashes();
        let capacity = batch.capacity();
        batch.clear();
        assert!(batch.is_empty());
        assert_eq!(batch.hashes(), None);
        assert_eq!(batch.capacity(), capacity);
        batch.fill_hashes();
        assert_eq!(batch.hashes(), Some([].as_slice()));
    }
}
