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
use std::hash::BuildHasher;
use std::marker::PhantomData;

use bitvec::array::BitArray;
use bitvec::view::BitView;
use num::{One, Zero};
use ppvm_traits_2::{Pauli, PauliBits, Word};

use crate::hash::structural_hash;
use crate::storage::{DefaultStorage, HashFinalize, PauliStorage};

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
/// The `hash_cache` is an **eager** plain `u64`: every constructor and every
/// structural mutator recomputes the finalized digest immediately, so
/// [`Indexable::key_hash`](ppvm_traits_2::Indexable::key_hash) is a field read
/// that never mutates. Because nothing is interior-mutable, the word stays
/// `Copy`.
///
/// This is a deliberate, measured departure from the design's §"Lazy hashing and
/// interior mutability" (which specifies a lazy cache and lists preserving `Copy`
/// as a non-goal). The lazy forms were both shipped first and then measured out:
///
/// * `OnceLock<u64>` (the design's sketch) — `Once`'s CAS init path dominated the
///   `PauliSum` Clifford re-key loop, where every freshly built key hits the cold
///   init exactly once.
/// * a relaxed-atomic sentinel — better, but the digest still fired *lazily*
///   inside the accumulate probe's `entry()`, i.e. on the bucket-index critical
///   path, where the finalize mul-chain stalls the dependent bucket load. Old
///   hashes eagerly in its first pass, so its probe hits a cached `u64` and the
///   hashing overlaps other terms' work. Hoisting the digest earlier took
///   `rotation_rx` from 1.07× to ~0.99× (`53ebc66e`); computing it eagerly at the
///   mutation boundary is the limit of that same move.
///
/// Losing interior mutability also restores `Copy`, which matters on the re-key
/// path: an atomic cache makes `clone` measurably more expensive (about `2.1×` a
/// plain-`u64` copy at 256 qubits) and forced a borrowed-source builder — see
/// [`PauliBits::PREFER_BORROWED_REKEY`](ppvm_traits_2::PauliBits::PREFER_BORROWED_REKEY),
/// which this word leaves `false` precisely because it is `Copy`.
///
/// `LossyPauliWord` and `Tableau` keep the lazy sentinel atomic, and that split is
/// intentional: their invalidation-to-read ratio is the opposite (every Clifford
/// gate invalidates a tableau; loss-only mutations must not rehash the X/Z
/// planes), so paying the digest up front would be wasted work there.
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
pub struct PauliWord<A: PauliStorage = DefaultStorage, H = fxhash::FxBuildHasher> {
    /// X-bit plane (one logical bit per qubit; unused high bits are `0`).
    pub(crate) xbits: BitArray<A>,
    /// Z-bit plane (one logical bit per qubit; unused high bits are `0`).
    pub(crate) zbits: BitArray<A>,
    /// Number of qubits (logical width).
    pub(crate) nqubits: usize,
    /// Eager finalized structural digest.
    pub(crate) hash_cache: u64,
    /// The private internal digest algorithm; never a runtime value.
    /// `fn() -> H` keeps `PauliWord` `Send + Sync` for any `H`.
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

/// The storage word holding logical bit `i`, and `i`'s offset inside it.
#[inline(always)]
fn word_of<A: PauliStorage>(i: usize) -> (usize, usize) {
    let bits = std::mem::size_of::<<A as BitView>::Store>() * 8;
    (i / bits, i % bits)
}

/// `1 << offset` when `toggle`, else `0` — a select on a register, never a
/// branch on memory.
#[inline(always)]
fn bit_mask<A: PauliStorage>(offset: usize, toggle: bool) -> <A as BitView>::Store {
    let one = <A as BitView>::Store::one();
    let zero = <A as BitView>::Store::zero();
    if toggle { one << offset } else { zero }
}

/// XOR `toggle_i` into bit `i` and `toggle_j` into bit `j` of the **same**
/// plane, touching each storage word once.
///
/// This is the *Clifford re-key* shape, and only that: the two-qubit re-keys
/// reach it with both toggles on one plane and both derived from bits of the
/// term's own word — `CZ` is `z ⊕= x` at each of its two sites — so `if toggle {
/// plane.set(i, !plane[i]) }` is a branch the predictor cannot learn, and when
/// both qubits share a storage word the second read-modify-write serializes
/// behind the first store. Building both masks and folding them into one XOR
/// removes both; the `wi == wj` test is loop-invariant across the re-key, so it
/// predicts perfectly.
///
/// XOR, not OR, when the words coincide: with `i == j` the two masks are the
/// same bit and must cancel, exactly as toggling one bit twice does. For
/// distinct bits the masks are disjoint and `^` is `|`.
///
/// The *rotation* branch builders deliberately do **not** come here — see
/// [`PauliWord::with_bits_toggled2`].
#[inline(always)]
fn xor_bits2<A: PauliStorage>(
    plane: &mut BitArray<A>,
    i: usize,
    toggle_i: bool,
    j: usize,
    toggle_j: bool,
) {
    let (wi, oi) = word_of::<A>(i);
    let (wj, oj) = word_of::<A>(j);
    let mask_i = bit_mask::<A>(oi, toggle_i);
    let mask_j = bit_mask::<A>(oj, toggle_j);
    let raw = plane.data.as_raw_mut_slice();
    if wi == wj {
        raw[wi] = raw[wi] ^ (mask_i ^ mask_j);
    } else {
        raw[wi] = raw[wi] ^ mask_i;
        raw[wj] = raw[wj] ^ mask_j;
    }
}

impl<A, H> PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    /// Construct the identity word `I…I` on `nqubits` qubits (all planes zero).
    ///
    /// Constructors compute the digest eagerly, so the word is immediately a
    /// valid map key (Design: `word-data-structures.md` §"Invalidation rules").
    #[inline]
    pub fn new(nqubits: usize) -> Self {
        debug_assert!(
            nqubits <= 8 * std::mem::size_of::<A>(),
            "nqubits {nqubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        Self::from_planes(BitArray::ZERO, BitArray::ZERO, nqubits)
    }

    /// Assemble from already-packed planes and a width.
    #[inline]
    pub(crate) fn from_planes(xbits: BitArray<A>, zbits: BitArray<A>, nqubits: usize) -> Self {
        let hash_cache = structural_hash::<A, H>(&xbits.data, &zbits.data, nqubits);
        Self {
            xbits,
            zbits,
            nqubits,
            hash_cache,
            _hasher: PhantomData,
        }
    }

    /// Recompute the structural-hash cache after a structural mutation
    /// (Design: §"Indexable values" — mutators maintain the affected private
    /// cache through `&mut self`).
    ///
    /// Named `invalidate_*` for continuity with the lazy design, but the cache is
    /// eager: this refreshes the digest rather than clearing it. See the type-level
    /// note on why the digest is computed at the mutation boundary.
    #[inline]
    pub(crate) fn invalidate_hash(&mut self) {
        self.hash_cache = structural_hash::<A, H>(&self.xbits.data, &self.zbits.data, self.nqubits);
    }

    /// A copy with the X and/or Z bit at `i` toggled, digest refreshed once.
    ///
    /// This is the rotation-branch key builder (`iGP` from a diagonal `P`). The
    /// branch always toggles a bit, so `clone()` followed by
    /// `set_x_bit`/`set_z_bit` would copy the source's digest and then immediately
    /// recompute it. Building the toggled key directly — copy the plane words,
    /// flip the bit, hash once — pays for exactly one digest.
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

    /// A copy with the X and/or Z bits at **two** sites toggled, digest refreshed
    /// once — the *two-qubit* rotation-branch key builder
    /// (`rzz`/`rxx`/`ryy`/`rotate_2`).
    ///
    /// One plane copy, up to four bit flips, one rebuild. Chaining
    /// [`with_bits_toggled`](PauliWord::with_bits_toggled) twice instead copies
    /// **both** planes and rebuilds the word *twice* per produced branch term, and
    /// the redundant half scales with the storage tier (64 extra bytes moved per
    /// branch at `[u8; 32]`). Old built one `k.clone()` and wrote four bits into it
    /// (`ppvm-pauli-sum/src/sum/rot2.rs`); this is that shape, minus the cache
    /// load+invalidate the clone would pay.
    ///
    /// Note it does **not** share [`xor_bits2`] with the Clifford re-key builder
    /// [`PauliBits::into_toggled_bits2`], even though the bit action is the same
    /// shape. Every rotation reaches this with *constant* toggles (`rzz` is
    /// `toggled_bits2(a, false, true, b, false, true)`), so the branch the
    /// masked form exists to remove is already gone at compile time and its
    /// `wi == wj` test is pure overhead. Sharing it moved
    /// `pauli_sum/integration_trotter` — ten Trotter steps of `rx`/`rzz` over a
    /// growing support — from 0.985x against old to 1.08–1.15x, reproducibly and
    /// across three executable layouts.
    #[inline]
    pub fn with_bits_toggled2(
        &self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self {
        debug_assert!(i < self.nqubits, "qubit {i} out of bounds");
        debug_assert!(j < self.nqubits, "qubit {j} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        if toggle_x_i {
            let b = xbits[i];
            xbits.set(i, !b);
        }
        if toggle_z_i {
            let b = zbits[i];
            zbits.set(i, !b);
        }
        if toggle_x_j {
            let b = xbits[j];
            xbits.set(j, !b);
        }
        if toggle_z_j {
            let b = zbits[j];
            zbits.set(j, !b);
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

impl<A, H> PauliBits for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline(always)]
    fn x_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits[i]
    }

    #[inline(always)]
    fn z_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits[i]
    }

    #[inline(always)]
    fn pauli_code(&self, i: usize) -> u8 {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        #[cfg(target_endian = "little")]
        {
            let byte = i >> 3;
            let shift = i & 7;
            let x = (bytemuck::bytes_of(&self.xbits.data)[byte] >> shift) & 1;
            let z = (bytemuck::bytes_of(&self.zbits.data)[byte] >> shift) & 1;
            x | (z << 1)
        }
        #[cfg(target_endian = "big")]
        {
            (self.xbits[i] as u8) | ((self.zbits[i] as u8) << 1)
        }
    }

    /// Set the X bit at `i`, then lazily invalidate the hash (Design:
    /// §"Representation types" — the structural mutation boundary that clears the
    /// affected hash component). Only in-range slots are touched, preserving the
    /// canonical-unused-bits invariant.
    ///
    /// Writing a bit its **current** value is a no-op, so the cached digest —
    /// a pure function of `(nqubits, X, Z)` — is still valid and is deliberately
    /// *kept*. This is not a micro-optimization on a cold path: the Clifford
    /// kernels write both target bits unconditionally (`CNOT` does
    /// `x_tgt ⊕= x_ctrl` even when `x_ctrl = 0`), and in a real circuit most
    /// terms are `I` at the gate's qubits, so an unconditional invalidation
    /// forced a full structural re-hash of nearly the whole support on every
    /// gate — recomputing a digest guaranteed to be bit-identical.
    #[inline(always)]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != v {
            self.xbits.set(i, v);
            self.invalidate_hash();
        }
    }

    #[inline(always)]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.zbits[i] != v {
            self.zbits.set(i, v);
            self.invalidate_hash();
        }
    }

    #[inline(always)]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != x || self.zbits[i] != z {
            self.xbits.set(i, x);
            self.zbits.set(i, z);
            self.invalidate_hash();
        }
    }

    #[inline(always)]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        if self.xbits[i] != xi || self.zbits[i] != zi || self.xbits[j] != xj || self.zbits[j] != zj
        {
            self.xbits.set(i, xi);
            self.zbits.set(i, zi);
            self.xbits.set(j, xj);
            self.zbits.set(j, zj);
            self.invalidate_hash();
        }
    }

    #[inline(always)]
    fn set_x_bit_and_z_bit(&mut self, x_i: usize, x: bool, z_i: usize, z: bool) {
        debug_assert!(
            x_i < self.nqubits && z_i < self.nqubits,
            "index out of bounds"
        );
        if self.xbits[x_i] != x || self.zbits[z_i] != z {
            self.xbits.set(x_i, x);
            self.zbits.set(z_i, z);
            self.invalidate_hash();
        }
    }

    #[inline(always)]
    fn set_z_bit_pair(&mut self, i: usize, zi: bool, j: usize, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        if i == j {
            // The trait default is two scalar sets, so a repeated index is
            // last-write-wins; the fused arm below would instead treat the two
            // requests as independent toggles.
            self.set_z_bit(j, zj);
            return;
        }
        let toggle_i = self.zbits[i] != zi;
        let toggle_j = self.zbits[j] != zj;
        if toggle_i || toggle_j {
            xor_bits2(&mut self.zbits, i, toggle_i, j, toggle_j);
            self.invalidate_hash();
        }
    }

    /// The direct plane-copy branch key builder — see
    /// [`with_bits_toggled`](PauliWord::with_bits_toggled). Overrides the trait's
    /// clone-then-flip default, which would load and then immediately invalidate
    /// the source's cached digest.
    #[inline]
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self {
        PauliWord::with_bits_toggled(self, i, toggle_x, toggle_z)
    }

    /// The direct plane-copy **two-site** branch key builder — see
    /// [`with_bits_toggled2`](PauliWord::with_bits_toggled2). Overrides the trait's
    /// clone-then-flip default so the two-qubit rotations copy the planes exactly
    /// once.
    #[inline]
    fn toggled_bits2(
        &self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self {
        PauliWord::with_bits_toggled2(self, i, toggle_x_i, toggle_z_i, j, toggle_x_j, toggle_z_j)
    }

    #[inline(always)]
    fn into_toggled_bits2(
        mut self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self {
        debug_assert!(i < self.nqubits, "qubit {i} out of bounds");
        debug_assert!(j < self.nqubits, "qubit {j} out of bounds");
        xor_bits2(&mut self.xbits, i, toggle_x_i, j, toggle_x_j);
        xor_bits2(&mut self.zbits, i, toggle_z_i, j, toggle_z_j);
        self.invalidate_hash();
        self
    }
}

impl<A: PauliStorage, H> Clone for PauliWord<A, H> {
    /// Cloning copies the (possibly cached) hash because the clone has identical
    /// structural contents (Design: `word-data-structures.md` §"Invalidation
    /// rules"). Hand-written so no spurious `H: Clone` bound is imposed.
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<A: PauliStorage, H> Copy for PauliWord<A, H> {}

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

impl<A, H> From<&str> for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
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

impl<A, H> From<String> for PauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::Indexable;

    /// The fused two-site toggle must equal two sequential one-site toggles for
    /// every index pair, including a repeated index (where the two requests
    /// cancel) and one straddling a storage word.
    #[test]
    fn toggled_bits2_matches_two_single_toggles() {
        let text = "XYZI".repeat(32);
        let base: PauliWord<[u8; 16]> = PauliWord::from(text.as_str());
        for (i, j) in [(1usize, 5usize), (5, 1), (3, 3), (2, 70), (70, 2), (70, 71)] {
            for bits in 0..16u8 {
                let (xi, zi) = (bits & 1 != 0, bits & 2 != 0);
                let (xj, zj) = (bits & 4 != 0, bits & 8 != 0);
                let fused = base.with_bits_toggled2(i, xi, zi, j, xj, zj);
                let chained = base
                    .with_bits_toggled(i, xi, zi)
                    .with_bits_toggled(j, xj, zj);
                assert_eq!(fused, chained, "({i},{j}) toggles {bits:04b}");
                assert_eq!(fused.key_hash(), chained.key_hash(), "digest {bits:04b}");
                assert_eq!(
                    base.into_toggled_bits2(i, xi, zi, j, xj, zj),
                    chained,
                    "owned ({i},{j}) toggles {bits:04b}"
                );
            }
        }
    }

    /// `set_z_bit_pair` must agree with two scalar `set_z_bit` calls — including
    /// when both sites share one storage word (the fused single-XOR arm) and
    /// when they straddle two (the split arm), and when neither bit moves (no
    /// write, and the digest must be unchanged rather than merely equal).
    #[test]
    fn set_z_bit_pair_matches_scalar_setters() {
        for (i, j) in [(1usize, 5usize), (5, 1), (3, 3), (2, 70), (70, 2), (70, 71)] {
            for zi in [false, true] {
                for zj in [false, true] {
                    // 128 sites: `[u8; 16]` is exactly two 64-bit halves, so
                    // sites 2 and 70 straddle a storage-word boundary.
                    let text = "XYZI".repeat(32);
                    let base: PauliWord<[u8; 16]> = PauliWord::from(text.as_str());
                    let mut fused = base;
                    fused.set_z_bit_pair(i, zi, j, zj);
                    let mut scalar = base;
                    scalar.set_z_bit(i, zi);
                    scalar.set_z_bit(j, zj);
                    assert_eq!(fused, scalar, "({i},{j}) <- ({zi},{zj})");
                    assert_eq!(
                        fused.key_hash(),
                        scalar.key_hash(),
                        "digest ({i},{j}) <- ({zi},{zj})"
                    );
                    assert_eq!(fused.x_bit(i), base.x_bit(i), "X plane must not move");
                    assert_eq!(fused.x_bit(j), base.x_bit(j), "X plane must not move");
                }
            }
        }
    }

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
