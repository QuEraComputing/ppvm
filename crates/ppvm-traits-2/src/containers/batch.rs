// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::containers::Indexable;
use crate::loss::LossState;
use crate::word::PauliBits;

/// An [`Indexable`] key with a structure-of-arrays column representation.
/// Hashing and column layout remain separate capabilities.
pub trait Columnar: Indexable {
    /// The concrete structure-of-arrays column for this key type.
    type Column: KeyColumn<Key = Self>;
}

/// Read access and value-producing operations on a structure-of-arrays key column.
/// See [`KeyColumnMut`] for construction and mutation.
pub trait KeyColumn: Default + Clone {
    /// The key type this column stores.
    type Key: Columnar;

    /// Number of keys currently in the column.
    fn len(&self) -> usize;

    /// Number of keys the existing plane allocations can hold without growing.
    fn capacity(&self) -> usize;

    /// Whether the column is empty.
    fn is_empty(&self) -> bool {
        self.len() == 0
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

    /// Materialize one row while toggling two sites with `[toggle_x, toggle_z]` masks.
    #[inline]
    fn toggled_bits2(
        &self,
        row: usize,
        i: usize,
        toggle_i: [bool; 2],
        j: usize,
        toggle_j: [bool; 2],
    ) -> Self::Key
    where
        Self::Key: PauliBits,
    {
        self.get(row).toggled_bits2(i, toggle_i, j, toggle_j)
    }
}

/// Construction and mutation of a [`KeyColumn`].
pub trait KeyColumnMut: KeyColumn {
    /// Create an empty column with capacity for `n` keys.
    fn with_capacity(n: usize) -> Self;

    /// Append one key while keeping each plane contiguous.
    fn push(&mut self, key: Self::Key);

    /// Reserve room for `additional` keys while retaining existing entries.
    /// The default is a no-op for backends without reservation support.
    #[inline]
    fn reserve(&mut self, _additional: usize) {}

    /// Clear the column while retaining its backing allocations.
    fn clear(&mut self);

    /// Overwrite element `i`, preserving other elements and backing allocations.
    /// Panics (or debug-panics) if `i >= len()`.
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

/// Keys and cached structural hashes in parallel columns, without coefficients.
/// Uses `Vec<W>` so keys need not implement [`Columnar`].
/// Mutations invalidate hashes; [`Self::hashes`] exposes cache validity.
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

/// A [`KeyBatch`] plus a separate coefficient column, ready for accumulation.
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

/// Append produced key/coefficient pairs to a sink.
/// Backends may collect scalar pairs or append to separate columns.
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
