// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli words with construction, inspection, and mutation.
//! Equality compares qubit count and X/Z planes; mutations refresh the hash.

use std::hash::BuildHasher;
use std::marker::PhantomData;

use bitvec::array::BitArray;
use bitvec::view::BitView;
use num::{One, Zero};
use ppvm_traits_2::{Pauli, PauliBits, Word};

use crate::hash::structural_hash;
use crate::{DefaultStorage, HashFinalize, PauliStorage};

/// A Pauli word in packed X/Z planes: I=(0,0), X=(1,0), Z=(0,1), Y=(1,1).
/// `A` selects backing storage; `H` selects the internal hash algorithm.
/// Unused bits stay zero; equality compares qubit count and both planes.
/// An eagerly refreshed `u64` hash makes lookup a field read and preserves `Copy`.
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

/// Toggles two bits in one plane, updating each affected storage word once.
/// Combines masks with XOR so toggling the same bit twice cancels.
/// Used for data-dependent Clifford toggles to avoid separate conditional writes.
/// Rotation builders use direct flips instead; see `PauliWord::with_bits_toggled2`.
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
    /// Constructs the identity word on `nqubits` qubits with zeroed planes.
    /// Computes the hash eagerly so the word is ready for map insertion.
    #[inline]
    pub fn new(nqubits: usize) -> Self {
        Self::from_planes(BitArray::ZERO, BitArray::ZERO, nqubits)
    }

    /// Assembles packed planes, rejecting widths beyond the backing storage capacity.
    #[inline]
    pub(crate) fn from_planes(xbits: BitArray<A>, zbits: BitArray<A>, nqubits: usize) -> Self {
        debug_assert!(
            nqubits <= 8 * std::mem::size_of::<A>(),
            "nqubits {nqubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        let hash_cache = structural_hash::<A, H>(&xbits.data, &zbits.data, nqubits);
        Self {
            xbits,
            zbits,
            nqubits,
            hash_cache,
            _hasher: PhantomData,
        }
    }

    /// Refreshes the cached hash immediately after a structural mutation.
    #[inline]
    pub(crate) fn refresh_hash(&mut self) {
        self.hash_cache = structural_hash::<A, H>(&self.xbits.data, &self.zbits.data, self.nqubits);
    }

    /// Copies the planes, toggles the requested X/Z bits at `i`, and hashes once.
    /// Builds a rotation-branch key without copying and then replacing the old hash.
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

    /// Copies both planes once, toggles up to four bits at two sites, and hashes once.
    /// Avoids rebuilding the word twice when constructing a two-qubit rotation branch.
    /// Keeps direct flips separate from `xor_bits2`: constant rotation toggles let
    /// the compiler remove branches without testing whether storage words coincide.
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
        debug_assert!(index < self.nqubits, "index {index} out of bounds");
        match (self.xbits[index], self.zbits[index]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        }
    }

    /// Counts non-identity factors by counting set bits in `x | z` in packed chunks.
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

    /// Sets the X bit at `i` and eagerly refreshes the hash only if the bit changes.
    /// Unchanged bits retain the cached hash, avoiding redundant work during gates.
    #[inline(always)]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != v {
            self.xbits.set(i, v);
            self.refresh_hash();
        }
    }

    #[inline(always)]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.zbits[i] != v {
            self.zbits.set(i, v);
            self.refresh_hash();
        }
    }

    #[inline(always)]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != x || self.zbits[i] != z {
            self.xbits.set(i, x);
            self.zbits.set(i, z);
            self.refresh_hash();
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
            self.refresh_hash();
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
            self.refresh_hash();
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
            self.refresh_hash();
        }
    }

    /// Uses `with_bits_toggled` to copy planes and hash once instead of clone-then-flip.
    #[inline]
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self {
        PauliWord::with_bits_toggled(self, i, toggle_x, toggle_z)
    }

    /// Uses `with_bits_toggled2` to copy planes once for a two-site rotation branch.
    #[inline]
    fn toggled_bits2(
        &self,
        i: usize,
        [toggle_x_i, toggle_z_i]: [bool; 2],
        j: usize,
        [toggle_x_j, toggle_z_j]: [bool; 2],
    ) -> Self {
        PauliWord::with_bits_toggled2(self, i, toggle_x_i, toggle_z_i, j, toggle_x_j, toggle_z_j)
    }

    #[inline(always)]
    fn into_toggled_bits2(
        mut self,
        i: usize,
        [toggle_x_i, toggle_z_i]: [bool; 2],
        j: usize,
        [toggle_x_j, toggle_z_j]: [bool; 2],
    ) -> Self {
        debug_assert!(i < self.nqubits, "qubit {i} out of bounds");
        debug_assert!(j < self.nqubits, "qubit {j} out of bounds");
        xor_bits2(&mut self.xbits, i, toggle_x_i, j, toggle_x_j);
        xor_bits2(&mut self.zbits, i, toggle_z_i, j, toggle_z_j);
        self.refresh_hash();
        self
    }
}

impl<A: PauliStorage, H> Clone for PauliWord<A, H> {
    /// Copies the word and its valid cached hash without requiring `H: Clone`.
    #[inline]
    fn clone(&self) -> Self {
        *self
    }
}

impl<A: PauliStorage, H> Copy for PauliWord<A, H> {}

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

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::Indexable;

    #[test]
    fn constructors_enforce_storage_capacity() {
        assert!(std::panic::catch_unwind(|| PauliWord::<u8>::new(9)).is_err());
        assert!(
            std::panic::catch_unwind(|| {
                PauliWord::<u8>::from_planes(BitArray::ZERO, BitArray::ZERO, 9)
            })
            .is_err()
        );
        for width in [0, 8] {
            let word = PauliWord::<u8>::new(width);
            assert_eq!(word.n_sites(), width);
            assert_eq!(word.to_string(), "I".repeat(width));
        }
    }

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
                let borrowed = base.toggled_bits2(i, [xi, zi], j, [xj, zj]);
                assert_eq!(borrowed, chained, "borrowed ({i},{j}) toggles {bits:04b}");
                assert_eq!(borrowed.key_hash(), chained.key_hash());
                assert_eq!(
                    base.into_toggled_bits2(i, [xi, zi], j, [xj, zj]),
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
                    // With byte-backed storage, sites 2 and 70 lie in different bytes.
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
