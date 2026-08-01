// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The packed [`PauliWord`] struct, its inherent constructors/accessors, the
//! read-only [`Word`] inspection impl, the [`PauliBits`] sub-site mutation impl,
//! and the structural [`PartialEq`]/[`Eq`]/[`Clone`]/[`Display`]/parsing that
//! agree on the logical identity `(nqubits, X bits, Z bits)`.
//!
//! Design: `word-data-structures.md` §"`PauliWord` packed representation" and
//! §"Logical Pauli model"; `traits-2-configuration-and-hashing.md`
//! §"Representation types" and §"Pauli algebra traits" (`PauliBits`).

use std::fmt;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU64, Ordering};

use bitvec::array::BitArray;
use ppvm_traits_2::{Pauli, PauliBits, Word};

use crate::storage::PauliStorage;

/// The [`hash_cache`](PauliWord::hash_cache) sentinel meaning "not yet computed".
///
/// A word whose finalized structural digest happens to equal this value is
/// simply recomputed on every `key_hash()` — the result is still the true
/// digest, only the caching is skipped, so this is a `1`-in-`2⁶⁴` perf non-event
/// with **no** correctness consequence (`key_hash()` and
/// `KeyColumn::hash_into` stay bit-for-bit identical either way).
pub(crate) const HASH_UNCACHED: u64 = 0;

/// A fixed-width Pauli word stored as parallel packed `x` and `z` bit planes.
///
/// Each qubit slot is encoded by two logical bits `(x, z)`, giving the four
/// Paulis `I = (0,0)`, `X = (1,0)`, `Z = (0,1)`, `Y = (1,1)` — the
/// stabilizer-formalism encoding of `word-data-structures.md` §"Logical Pauli
/// model". `A` is the backing-storage blob (e.g. `u64` packs up to 64 qubits,
/// `[u64; 4]` up to 256) and `H` is the private internal digest algorithm used
/// to compute `key_hash()`; **both `A` and `H` are private representation
/// parameters, neither exposed through `Word` or `Indexable`**
/// (`word-data-structures.md` §"`PauliWord` packed representation").
///
/// The structural identity is `(nqubits, logical X bits, logical Z bits)`.
/// Equality and hashing exclude unused capacity, the cache, and the
/// `PhantomData` marker; this is upheld by the **canonical-unused-bits
/// invariant**: every mutator touches only in-range qubit slots (and the twisted
/// product XORs zero against zero), so unused high bits are permanently `0` and
/// comparing the full backing blob is equivalent to comparing only the logical
/// bits.
///
/// The `hash_cache` is a lazy `AtomicU64` (Design: §"Lazy hashing and interior
/// mutability" — the design sketches this with an `OnceLock<u64>`; a sentinel
/// `AtomicU64` realizes the *same* contract — lazy, interior-mutable, `Send +
/// Sync` — in half the width and with a plain relaxed load/store instead of
/// `Once`'s CAS init path, which measurably dominated the Clifford re-key hot
/// loop where every freshly built key hit the cold init once): `key_hash()` may
/// populate it through `&self` with a relaxed store, and each structural mutator
/// resets it to [`HASH_UNCACHED`] through `&mut self`. `Copy` is intentionally
/// dropped for correct lazy caching.
///
/// # Examples
///
/// ```
/// use ppvm_pauli_word_2::PauliWord;
/// use ppvm_traits_2::{Pauli, Word, PauliBits};
///
/// let w: PauliWord = "XYZI".into();
/// assert_eq!(w.n_sites(), 4);
/// assert_eq!(w.get(0), Pauli::X);
/// assert_eq!(w.get(3), Pauli::I);
/// assert_eq!(w.weight(), 3);
///
/// let mut w2: PauliWord = PauliWord::new(2);
/// w2.set_x_bit(0, true); // X on qubit 0
/// w2.set_z_bit(1, true); // Z on qubit 1
/// assert_eq!(w2.to_string(), "XZ");
/// ```
pub struct PauliWord<A: PauliStorage = u64, H = fxhash::FxBuildHasher> {
    /// X-bit plane (one logical bit per qubit; unused high bits are `0`).
    pub(crate) xbits: BitArray<A>,
    /// Z-bit plane (one logical bit per qubit; unused high bits are `0`).
    pub(crate) zbits: BitArray<A>,
    /// Number of qubits (logical width).
    pub(crate) nqubits: usize,
    /// Lazy structural-hash cache (Design: §"Lazy hashing and interior
    /// mutability"). Holds [`HASH_UNCACHED`] until `key_hash()` first populates
    /// it with the finalized digest.
    pub(crate) hash_cache: AtomicU64,
    /// The private internal digest algorithm; never a runtime value.
    /// `fn() -> H` keeps `PauliWord` `Send + Sync` for any `H`.
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> PauliWord<A, H> {
    /// Construct the identity word `I…I` on `nqubits` qubits (all planes zero).
    ///
    /// Constructors leave the cache empty (Design: `word-data-structures.md`
    /// §"Invalidation rules").
    #[inline]
    pub fn new(nqubits: usize) -> Self {
        debug_assert!(
            nqubits <= 8 * std::mem::size_of::<A>(),
            "nqubits {nqubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        Self {
            xbits: BitArray::ZERO,
            zbits: BitArray::ZERO,
            nqubits,
            hash_cache: AtomicU64::new(HASH_UNCACHED),
            _hasher: PhantomData,
        }
    }

    /// Assemble from already-packed planes and a width. Callers must uphold the
    /// canonical-unused-bits invariant (bits at index `>= nqubits` are `0`); the
    /// cache is left empty. Crate-internal — used by the twisted product and the
    /// key column, which both preserve the invariant.
    #[inline]
    pub(crate) fn from_planes(xbits: BitArray<A>, zbits: BitArray<A>, nqubits: usize) -> Self {
        Self {
            xbits,
            zbits,
            nqubits,
            hash_cache: AtomicU64::new(HASH_UNCACHED),
            _hasher: PhantomData,
        }
    }

    /// Clear the lazy structural-hash cache after a structural mutation
    /// (Design: §"Indexable values" — mutators clear the affected private cache
    /// through `&mut self`).
    #[inline]
    pub(crate) fn invalidate_hash(&mut self) {
        // Exclusive `&mut self` access — reset without an atomic RMW; the next
        // `key_hash()` recomputes.
        *self.hash_cache.get_mut() = HASH_UNCACHED;
    }

    /// A copy with the X and/or Z bit at `i` toggled and a fresh **uncached**
    /// hash (recomputed on the next `key_hash()`).
    ///
    /// This is the rotation-branch key builder (`iGP` from a diagonal `P`). The
    /// branch always toggles a bit, so `clone()` — which copies the source's
    /// *cached* digest — followed by `set_x_bit`/`set_z_bit` — which immediately
    /// invalidates it — does a wasted atomic load + store per branch key. Building
    /// the toggled key directly (copy the plane words, flip the bit, leave the
    /// cache empty) skips both.
    #[inline]
    pub fn with_bits_toggled(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self {
        debug_assert!(i < self.nqubits, "qubit {i} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        if toggle_x {
            let b = xbits[i];
            xbits.set(i, !b);
        }
        if toggle_z {
            let b = zbits[i];
            zbits.set(i, !b);
        }
        Self::from_planes(xbits, zbits, self.nqubits)
    }
}

impl<A: PauliStorage, H> Word for PauliWord<A, H> {
    type Site = Pauli;

    #[inline]
    fn n_sites(&self) -> usize {
        self.nqubits
    }

    #[inline]
    fn get(&self, index: usize) -> Pauli {
        assert!(index < self.nqubits, "index {index} out of bounds");
        match (self.xbits[index], self.zbits[index]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        }
    }

    /// Number of non-identity factors: a fused popcount of `x | z` over the
    /// packed planes. Ported verbatim from `ppvm-pauli-word` (one `popcnt` per
    /// machine word) to keep `weight` at parity.
    #[inline]
    fn weight(&self) -> usize {
        let xs: &[u8] = bytemuck::bytes_of(&self.xbits.data);
        let zs: &[u8] = bytemuck::bytes_of(&self.zbits.data);
        debug_assert_eq!(xs.len(), zs.len());

        let mut total: u32 = 0;
        let (mut i, n) = (0usize, xs.len());

        while i + 8 <= n {
            let x = u64::from_ne_bytes(xs[i..i + 8].try_into().unwrap());
            let z = u64::from_ne_bytes(zs[i..i + 8].try_into().unwrap());
            total += (x | z).count_ones();
            i += 8;
        }
        if i + 4 <= n {
            let x = u32::from_ne_bytes(xs[i..i + 4].try_into().unwrap());
            let z = u32::from_ne_bytes(zs[i..i + 4].try_into().unwrap());
            total += (x | z).count_ones();
            i += 4;
        }
        if i + 2 <= n {
            let x = u16::from_ne_bytes(xs[i..i + 2].try_into().unwrap());
            let z = u16::from_ne_bytes(zs[i..i + 2].try_into().unwrap());
            total += (x | z).count_ones();
            i += 2;
        }
        if i < n {
            total += (xs[i] | zs[i]).count_ones();
        }

        total as usize
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        (0..self.nqubits).map(move |i| self.get(i))
    }
}

impl<A: PauliStorage, H> PauliBits for PauliWord<A, H> {
    #[inline]
    fn x_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits[i]
    }

    #[inline]
    fn z_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits[i]
    }

    /// Set the X bit at `i`, then lazily invalidate the hash (Design:
    /// §"Representation types" — the structural mutation boundary that clears the
    /// affected hash component). Only in-range slots are touched, preserving the
    /// canonical-unused-bits invariant.
    #[inline]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits.set(i, v);
        self.invalidate_hash();
    }

    #[inline]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits.set(i, v);
        self.invalidate_hash();
    }
}

impl<A: PauliStorage, H> Clone for PauliWord<A, H> {
    /// Cloning copies the (possibly cached) hash because the clone has identical
    /// structural contents (Design: `word-data-structures.md` §"Invalidation
    /// rules"). Hand-written so no spurious `H: Clone` bound is imposed.
    #[inline]
    fn clone(&self) -> Self {
        Self {
            xbits: self.xbits,
            zbits: self.zbits,
            nqubits: self.nqubits,
            // Copy the (possibly still-uncached) digest verbatim; the clone has
            // identical structural contents.
            hash_cache: AtomicU64::new(self.hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}

impl<A: PauliStorage, H> fmt::Debug for PauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PauliWord")
            .field("nqubits", &self.nqubits)
            .field("word", &self.to_string())
            .finish()
    }
}

/// Structural equality over `(nqubits, X bits, Z bits)`; the cache and the
/// `PhantomData` marker are excluded. The canonical-unused-bits invariant makes
/// the full-blob comparison equivalent to comparing only the logical bits.
impl<A: PauliStorage, H> PartialEq for PauliWord<A, H> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.nqubits == other.nqubits
            && self.xbits.data == other.xbits.data
            && self.zbits.data == other.zbits.data
    }
}

impl<A: PauliStorage, H> Eq for PauliWord<A, H> {}

impl<A: PauliStorage, H> fmt::Display for PauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nqubits {
            let c = match self.get(i) {
                Pauli::I => 'I',
                Pauli::X => 'X',
                Pauli::Y => 'Y',
                Pauli::Z => 'Z',
            };
            f.write_str(c.encode_utf8(&mut [0u8; 1]))?;
        }
        Ok(())
    }
}

impl<A: PauliStorage, H> From<&str> for PauliWord<A, H> {
    /// Parse a Pauli string of `I`/`X`/`Y`/`Z` symbols (underscores are ignored
    /// separators), mirroring `ppvm-pauli-word`'s `From<&str>`. Panics on any
    /// other character or if the width exceeds the backing storage.
    fn from(value: &str) -> Self {
        let mut xbits = BitArray::<A>::ZERO;
        let mut zbits = BitArray::<A>::ZERO;
        let mut i = 0usize;
        for ch in value.chars() {
            match ch {
                'I' => {}
                'X' => xbits.set(i, true),
                'Z' => zbits.set(i, true),
                'Y' => {
                    xbits.set(i, true);
                    zbits.set(i, true);
                }
                '_' => continue,
                other => panic!("invalid Pauli character: {other}"),
            }
            i += 1;
        }
        Self::from_planes(xbits, zbits, i)
    }
}

impl<A: PauliStorage, H> From<String> for PauliWord<A, H> {
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_get_roundtrip() {
        let w: PauliWord = "XZYI".into();
        assert_eq!(w.n_sites(), 4);
        assert_eq!(w.get(0), Pauli::X);
        assert_eq!(w.get(1), Pauli::Z);
        assert_eq!(w.get(2), Pauli::Y);
        assert_eq!(w.get(3), Pauli::I);
        assert_eq!(w.to_string(), "XZYI");
    }

    #[test]
    fn weight_counts_nonidentity() {
        let w: PauliWord = "XIYZI".into();
        assert_eq!(w.weight(), 3);
        assert_eq!(PauliWord::<u64>::new(5).weight(), 0);
    }

    #[test]
    fn iter_matches_get() {
        let w: PauliWord = "XYZI".into();
        let via_iter: Vec<Pauli> = w.iter().collect();
        let via_get: Vec<Pauli> = (0..w.n_sites()).map(|i| w.get(i)).collect();
        assert_eq!(via_iter, via_get);
    }

    #[test]
    fn set_bits_build_word() {
        let mut w: PauliWord = PauliWord::new(3);
        w.set_x_bit(0, true);
        w.set_z_bit(0, true); // Y
        w.set_z_bit(2, true); // Z
        assert_eq!(w.get(0), Pauli::Y);
        assert_eq!(w.get(1), Pauli::I);
        assert_eq!(w.get(2), Pauli::Z);
        assert!(w.x_bit(0) && w.z_bit(0));
        assert!(!w.x_bit(2) && w.z_bit(2));
    }

    #[test]
    fn underscore_is_ignored() {
        let a: PauliWord = "X_Y_Z".into();
        let b: PauliWord = "XYZ".into();
        assert_eq!(a, b);
    }

    #[test]
    fn equality_excludes_width_mismatch() {
        let a: PauliWord = "XY".into();
        let b: PauliWord = "XYI".into();
        assert_ne!(a, b, "different widths are structurally distinct");
    }

    #[test]
    fn send_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<PauliWord>();
    }
}
