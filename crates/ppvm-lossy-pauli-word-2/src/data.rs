// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{
    fmt,
    marker::PhantomData,
    sync::atomic::{AtomicU64, Ordering},
};

use bitvec::array::BitArray;
use ppvm_pauli_word_2::{DefaultStorage, PauliStorage};
use ppvm_traits_2::{LossState, LossySite, Pauli, PauliBits, Word};

pub struct LossyPauliWord<A: PauliStorage = DefaultStorage, H = fxhash::FxBuildHasher> {
    pub(crate) xbits: BitArray<A>,
    pub(crate) zbits: BitArray<A>,
    /// Loss-bit plane: `1` at index `q` marks qubit `q` lost.
    pub(crate) lbits: BitArray<A>,
    /// Number of qubits (logical width).
    pub(crate) nqubits: usize,
    /// Lazy cache for the lossy structural digest.
    pub(crate) loss_hash_cache: AtomicU64,
    /// Lazy cache for the x/z structure digest.
    pub(crate) xz_hash_cache: AtomicU64,
    /// The private internal digest algorithm; `fn() -> H` keeps `Send + Sync`.
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> LossyPauliWord<A, H> {
    #[inline]
    pub fn new(nqubits: usize) -> Self {
        debug_assert!(
            nqubits <= 8 * std::mem::size_of::<A>(),
            "nqubits {nqubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        Self::from_planes(BitArray::ZERO, BitArray::ZERO, BitArray::ZERO, nqubits)
    }

    #[inline]
    pub(crate) fn from_planes(
        xbits: BitArray<A>,
        zbits: BitArray<A>,
        lbits: BitArray<A>,
        nqubits: usize,
    ) -> Self {
        Self {
            xbits,
            zbits,
            lbits,
            nqubits,
            loss_hash_cache: AtomicU64::new(0),
            xz_hash_cache: AtomicU64::new(0),
            _hasher: PhantomData,
        }
    }

    // /// clear the structural digest after word mutation
    #[inline]
    pub(crate) fn invalidate_xz_hash(&mut self) {
        *self.xz_hash_cache.get_mut() = 0;
    }

    #[inline]
    pub(crate) fn invalidate_loss_hash(&mut self) {
        *self.loss_hash_cache.get_mut() = 0;
    }

    /// check for whether a qubit was lost
    #[inline]
    pub fn is_lost(&self, qb: usize) -> bool {
        debug_assert!(qb < self.nqubits, "qubit {qb} out of bounds");
        self.lbits[qb]
    }

    #[inline]
    pub fn set_lost(&mut self, qb: usize) {
        debug_assert!(qb < self.nqubits, "qubit {qb} out of bounds");
        if self.xbits[qb] || self.zbits[qb] {
            self.xbits.set(qb, false);
            self.zbits.set(qb, false);
            self.invalidate_xz_hash();
        }
        if !self.lbits[qb] {
            self.lbits.set(qb, true);
            self.invalidate_loss_hash();
        }
    }

    // clear loss at qubit qb, returning it to identity I. The X/Z component
    // is preserved.
    #[inline]
    pub fn clear_loss(&mut self, qb: usize) {
        debug_assert!(qb < self.nqubits, "qubit {qb} out of bounds");
        if self.lbits[qb] {
            self.lbits.set(qb, false);
            self.invalidate_loss_hash();
        }
    }

    /// Build the reset-channel branch with `q` marked lost.
    pub fn with_lost(&self, qb: usize) -> Self {
        debug_assert!(qb < self.nqubits, "qubit {qb} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        let mut lbits = self.lbits;
        xbits.set(qb, false);
        zbits.set(qb, false);
        lbits.set(qb, true);
        Self::from_planes(xbits, zbits, lbits, self.nqubits)
    }

    pub fn set(&mut self, qb: usize, site: LossySite<Pauli>) {
        debug_assert!(qb < self.nqubits, "qubit {qb} out of bounds");
        match site {
            LossySite::Present(pauli) => {
                self.clear_loss(qb);
                let (x, z) = match pauli {
                    Pauli::I => (false, false),
                    Pauli::X => (true, false),
                    Pauli::Z => (false, true),
                    Pauli::Y => (true, true),
                };
                self.set_xz_bits(qb, x, z);
            }
            LossySite::Lost => self.set_lost(qb),
        }
    }

    /// Number of lost bits
    #[inline]
    pub fn loss_weight(&self) -> usize {
        self.lbits.count_ones()
    }
}

impl<A: PauliStorage, H> Word for LossyPauliWord<A, H> {
    type Site = LossySite<Pauli>;

    #[inline]
    fn n_sites(&self) -> usize {
        self.nqubits
    }

    /// Retrieve Pauli type for some index
    fn get(&self, idx: usize) -> Self::Site {
        debug_assert!(idx < self.nqubits, "index {idx} out of bounds");
        if self.lbits[idx] {
            return LossySite::Lost;
        }
        let pauli = match (self.xbits[idx], self.zbits[idx]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        };
        LossySite::Present(pauli)
    }

    /// Number of non-identity factors.
    #[inline]
    fn weight(&self) -> usize {
        let xs = bytemuck::bytes_of(&self.xbits.data);
        let zs = bytemuck::bytes_of(&self.zbits.data);
        let ls = bytemuck::bytes_of(&self.lbits.data);

        let mut total: u32 = 0;
        let (mut i, n) = (0usize, xs.len());

        // Count 64 qubits at a time, then handle the remaining bytes.
        while i + 8 <= n {
            let x = u64::from_ne_bytes(xs[i..i + 8].try_into().unwrap());
            let z = u64::from_ne_bytes(zs[i..i + 8].try_into().unwrap());
            let l = u64::from_ne_bytes(ls[i..i + 8].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 8;
        }
        if i + 4 <= n {
            let x = u32::from_ne_bytes(xs[i..i + 4].try_into().unwrap());
            let z = u32::from_ne_bytes(zs[i..i + 4].try_into().unwrap());
            let l = u32::from_ne_bytes(ls[i..i + 4].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 4;
        }
        if i + 2 <= n {
            let x = u16::from_ne_bytes(xs[i..i + 2].try_into().unwrap());
            let z = u16::from_ne_bytes(zs[i..i + 2].try_into().unwrap());
            let l = u16::from_ne_bytes(ls[i..i + 2].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 2;
        }
        if i < n {
            total += (xs[i] | zs[i] | ls[i]).count_ones();
        }

        total as usize
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = Self::Site> {
        (0..self.nqubits).map(move |i| self.get(i))
    }
}

impl<A: PauliStorage, H> Clone for LossyPauliWord<A, H> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            xbits: self.xbits,
            zbits: self.zbits,
            lbits: self.lbits,
            nqubits: self.nqubits,
            xz_hash_cache: AtomicU64::new(self.xz_hash_cache.load(Ordering::Relaxed)),
            loss_hash_cache: AtomicU64::new(self.loss_hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}

impl<A: PauliStorage, H> fmt::Display for LossyPauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nqubits {
            let c = match self.get(i) {
                LossySite::Present(Pauli::I) => 'I',
                LossySite::Present(Pauli::X) => 'X',
                LossySite::Present(Pauli::Y) => 'Y',
                LossySite::Present(Pauli::Z) => 'Z',
                LossySite::Lost => 'L',
            };
            f.write_str(c.encode_utf8(&mut [0u8; 1]))?;
        }
        Ok(())
    }
}

impl<A: PauliStorage, H> fmt::Debug for LossyPauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LossyPauliWord")
            .field("nqubits", &self.nqubits)
            .field("word", &self.to_string())
            .finish()
    }
}

impl<A: PauliStorage, H> PartialEq for LossyPauliWord<A, H> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.nqubits == other.nqubits
            && self.xbits.data == other.xbits.data
            && self.zbits.data == other.zbits.data
            && self.lbits.data == other.lbits.data
    }
}

impl<A: PauliStorage, H> Eq for LossyPauliWord<A, H> {}

impl<A: PauliStorage, H> From<&str> for LossyPauliWord<A, H> {
    /// Parse I/X/Y/Z/L symbols, ignoring underscores.
    /// Panics on invalid characters or if the backing capacity is exceeded.
    fn from(value: &str) -> Self {
        let mut xbits = BitArray::<A>::ZERO;
        let mut zbits = BitArray::<A>::ZERO;
        let mut lbits = BitArray::<A>::ZERO;

        let symbols = value.chars().filter(|&ch| ch != '_').collect::<Vec<char>>();
        for (i, ch) in symbols.iter().enumerate() {
            match ch {
                'I' => {}
                'X' => xbits.set(i, true),
                'Z' => zbits.set(i, true),
                'Y' => {
                    xbits.set(i, true);
                    zbits.set(i, true);
                }
                'L' => lbits.set(i, true),
                other => panic!("invalid lossy Pauli character: {other}"),
            }
        }
        let nqubits = symbols.len();

        Self::from_planes(xbits, zbits, lbits, nqubits)
    }
}

impl<A: PauliStorage, H> From<String> for LossyPauliWord<A, H> {
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

impl<A: PauliStorage, H> PauliBits for LossyPauliWord<A, H> {
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

    #[inline]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != v {
            self.xbits.set(i, v);
            self.invalidate_xz_hash();
        }
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
    }

    #[inline]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.zbits[i] != v {
            self.zbits.set(i, v);
            self.invalidate_xz_hash();
        }
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
    }

    #[inline]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        if self.xbits[i] != x || self.zbits[i] != z {
            self.xbits.set(i, x);
            self.zbits.set(i, z);
            self.invalidate_xz_hash();
        }
        if (x || z) && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
    }

    #[inline]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        let xz_changed = if i == j {
            self.xbits[j] != xj || self.zbits[j] != zj
        } else {
            self.xbits[i] != xi || self.zbits[i] != zi || self.xbits[j] != xj || self.zbits[j] != zj
        };
        if xz_changed {
            self.xbits.set(i, xi);
            self.zbits.set(i, zi);
            self.xbits.set(j, xj);
            self.zbits.set(j, zj);
            self.invalidate_xz_hash();
        }
        if (xi || zi) && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
        if (xj || zj) && self.lbits[j] {
            self.lbits.set(j, false);
            self.invalidate_loss_hash();
        }
    }

    #[inline]
    fn set_z_bit_pair(&mut self, i: usize, zi: bool, j: usize, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        let z_changed = if i == j {
            self.zbits[j] != zj
        } else {
            self.zbits[i] != zi || self.zbits[j] != zj
        };
        if z_changed {
            self.zbits.set(i, zi);
            self.zbits.set(j, zj);
            self.invalidate_xz_hash();
        }
        if zi && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
        if zj && self.lbits[j] {
            self.lbits.set(j, false);
            self.invalidate_loss_hash();
        }
    }

    // Build a rotation branch directly from the packed planes.
    #[inline]
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self
    where
        Self: Sized + Clone,
    {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        if toggle_x {
            let bit = xbits[i];
            xbits.set(i, !bit);
        }
        if toggle_z {
            let bit = zbits[i];
            zbits.set(i, !bit);
        }
        let mut lbits = self.lbits;
        if xbits[i] || zbits[i] {
            lbits.set(i, false);
        }
        Self::from_planes(xbits, zbits, lbits, self.nqubits)
    }

    /// Build a two-site rotation branch with one packed-plane copy.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn toggled_bits2(&self, i: usize, toggle_i: [bool; 2], j: usize, toggle_j: [bool; 2]) -> Self {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        let toggle_x_i = toggle_i[0];
        let toggle_x_j = toggle_j[0];
        let toggle_z_i = toggle_i[1];
        let toggle_z_j = toggle_j[1];
        if toggle_x_i {
            let bit = xbits[i];
            xbits.set(i, !bit);
        }
        if toggle_z_i {
            let bit = zbits[i];
            zbits.set(i, !bit);
        }
        if toggle_x_j {
            let bit = xbits[j];
            xbits.set(j, !bit);
        }
        if toggle_z_j {
            let bit = zbits[j];
            zbits.set(j, !bit);
        }
        let mut lbits = self.lbits;
        if xbits[i] || zbits[i] {
            lbits.set(i, false);
        }
        if xbits[j] || zbits[j] {
            lbits.set(j, false);
        }
        Self::from_planes(xbits, zbits, lbits, self.nqubits)
    }
}

impl<A: PauliStorage, H> LossState for LossyPauliWord<A, H> {
    #[inline]
    fn is_lost(&self, i: usize) -> bool {
        LossyPauliWord::is_lost(self, i)
    }

    #[inline]
    fn set_lost(&mut self, i: usize) {
        LossyPauliWord::set_lost(self, i);
    }

    #[inline]
    fn clear_lost(&mut self, i: usize) {
        LossyPauliWord::clear_loss(self, i);
    }

    #[inline]
    fn loss_cleared(&self, i: usize) -> Self
    where
        Self: Sized + Clone,
    {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        let mut lbits = self.lbits;
        lbits.set(i, false);
        Self::from_planes(self.xbits, self.zbits, lbits, self.nqubits)
    }

    #[inline]
    fn loss_weight(&self) -> usize {
        LossyPauliWord::loss_weight(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unchanged_setters_preserve_both_caches() {
        let mut word = LossyPauliWord::<u8>::from("YZL");
        // Seed nonzero cache values to observe invalidation without hashing.
        *word.xz_hash_cache.get_mut() = 11;
        *word.loss_hash_cache.get_mut() = 22;

        word.set(0, LossySite::Present(Pauli::Y));
        word.set_x_bit(0, true);
        word.set_z_bit(0, true);
        word.set_xz_bits(0, true, true);
        word.set_xz_bits2(0, true, true, 1, false, true);
        word.set_z_bit_pair(0, true, 1, true);
        word.set_x_bit_and_z_bit(0, true, 1, true);
        word.set_xz_bits(2, false, false);
        word.set_lost(2);
        word.clear_loss(0);

        assert_eq!(word.to_string(), "YZL");
        assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
        assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
    }

    #[test]
    fn changed_setters_invalidate_only_affected_components() {
        let setters: [fn(&mut LossyPauliWord<u8>); 6] = [
            |w| w.set(0, LossySite::Present(Pauli::Y)),
            |w| w.set_x_bit(0, true),
            |w| w.set_z_bit(0, true),
            |w| w.set_xz_bits(0, true, true),
            |w| w.set_xz_bits2(0, true, true, 1, false, true),
            |w| w.set_z_bit_pair(0, true, 1, true),
        ];
        for setter in setters {
            for input in ["II", "LL"] {
                let mut word = LossyPauliWord::<u8>::from(input);
                *word.xz_hash_cache.get_mut() = 11;
                *word.loss_hash_cache.get_mut() = 22;

                setter(&mut word);

                assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 0);
                let expected_loss_cache = if input == "LL" { 0 } else { 22 };
                assert_eq!(
                    word.loss_hash_cache.load(Ordering::Relaxed),
                    expected_loss_cache
                );
            }
        }
        let mut word = LossyPauliWord::<u8>::from("L");
        *word.xz_hash_cache.get_mut() = 11;
        *word.loss_hash_cache.get_mut() = 22;
        word.set(0, LossySite::Present(Pauli::I));
        assert_eq!(word.to_string(), "I");
        assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
        assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn repeated_indices_preserve_cache_when_final_bits_are_unchanged() {
        let mut word = LossyPauliWord::<u8>::from("I");
        *word.xz_hash_cache.get_mut() = 11;

        word.set_xz_bits2(0, true, true, 0, false, false);
        word.set_z_bit_pair(0, true, 0, false);

        assert_eq!(word.to_string(), "I");
        assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
    }

    #[test]
    fn weight_counts_correctly_across_storage_and_chunk_boundaries() {
        fn check<const N: usize>() {
            for width in 0..=N * 8 {
                let input: String = "IXYZL".chars().cycle().take(width).collect();
                let word = LossyPauliWord::<[u8; N]>::from(input.as_str());
                let expected = input.chars().filter(|&ch| ch != 'I').count();

                assert_eq!(word.weight(), expected, "storage={N}, width={width}");
            }
        }

        check::<1>();
        check::<2>();
        check::<4>();
        check::<8>();
        check::<15>();
        check::<16>();
        check::<17>();
    }

    #[test]
    fn toggling_a_lost_site_restores_it_only_when_bits_change() {
        let word = LossyPauliWord::<u8>::from("LL");

        assert_eq!(word.toggled_bits(0, true, false).to_string(), "XL");
        assert_eq!(word.toggled_bits(0, false, true).to_string(), "ZL");
        assert_eq!(word.toggled_bits(0, true, true).to_string(), "YL");
        assert_eq!(word.toggled_bits(0, false, false), word);
        assert_eq!(word.to_string(), "LL");
    }

    #[test]
    fn toggling_two_sites_clears_only_the_affected_loss_flags() {
        let word = LossyPauliWord::<u8>::from("LLL");

        let branch = word.toggled_bits2(0, [true, false], 2, [false, true]);
        assert_eq!(branch.to_string(), "XLZ");
        assert_eq!(branch.loss_weight(), 1);

        let branch = word.toggled_bits2(0, [false, false], 2, [true, true]);
        assert_eq!(branch.to_string(), "LLY");
        assert_eq!(word.to_string(), "LLL");
    }

    #[test]
    fn new_word_is_all_identity() {
        let word = LossyPauliWord::<u8>::new(3);

        assert_eq!(word.n_sites(), 3);
        assert_eq!(word.to_string(), "III");
        assert_eq!(word.weight(), 0);
        assert_eq!(word.loss_weight(), 0);
    }

    #[test]
    fn parsing_ignores_separators_and_preserves_site_order() {
        let word = LossyPauliWord::<u8>::from("_IX_YZ_L_");

        assert_eq!(word.n_sites(), 5);
        assert_eq!(word.to_string(), "IXYZL");
        assert_eq!(
            word.iter().collect::<Vec<_>>(),
            vec![
                LossySite::Present(Pauli::I),
                LossySite::Present(Pauli::X),
                LossySite::Present(Pauli::Y),
                LossySite::Present(Pauli::Z),
                LossySite::Lost,
            ]
        );
    }

    #[test]
    fn weights_count_nonidentity_and_lost_sites() {
        let word = LossyPauliWord::<u8>::from("IXYZLLI");

        assert_eq!(word.weight(), 5);
        assert_eq!(word.loss_weight(), 2);
    }

    #[test]
    fn marking_a_site_lost_removes_its_pauli() {
        let mut word = LossyPauliWord::<u8>::from("XYZ");

        word.set_lost(1);

        assert_eq!(word.to_string(), "XLZ");
        assert!(word.is_lost(1));
        assert!(!word.x_bit(1));
        assert!(!word.z_bit(1));
    }

    #[test]
    fn clearing_loss_returns_identity_and_preserves_other_sites() {
        let mut word = LossyPauliWord::<u8>::from("XLZ");

        word.clear_loss(1);
        word.clear_loss(0);

        assert_eq!(word.to_string(), "XIZ");
        assert!(!word.is_lost(1));
        assert_eq!(word.loss_weight(), 0);
    }

    #[test]
    fn setting_a_site_replaces_loss_or_pauli() {
        let mut word = LossyPauliWord::<u8>::from("LX");

        word.set(0, LossySite::Present(Pauli::Y));
        word.set(1, LossySite::Lost);

        assert_eq!(word.to_string(), "YL");
        assert!(!word.is_lost(0));
        assert!(word.is_lost(1));
    }

    #[test]
    fn with_lost_leaves_the_original_unchanged() {
        let word = LossyPauliWord::<u8>::from("XYZ");

        let branch = word.with_lost(1);

        assert_eq!(branch.to_string(), "XLZ");
        assert_eq!(word.to_string(), "XYZ");
    }

    #[test]
    fn changing_a_clone_leaves_the_original_unchanged() {
        let word = LossyPauliWord::<u8>::from("XLZ");
        let mut copy = word.clone();
        assert_eq!(copy, word);

        copy.clear_loss(1);
        copy.set(0, LossySite::Present(Pauli::Y));

        assert_eq!(copy.to_string(), "YIZ");
        assert_eq!(word.to_string(), "XLZ");
    }

    #[test]
    fn setting_pauli_bits_restores_lost_sites() {
        let mut word = LossyPauliWord::<u8>::from("LLL");

        word.set_x_bit(0, true);
        word.set_z_bit(1, true);
        word.set_xz_bits(2, true, true);

        assert_eq!(word.to_string(), "XZY");
        assert_eq!(word.loss_weight(), 0);
    }

    #[test]
    fn equality_distinguishes_paulis_loss_and_width() {
        let word = LossyPauliWord::<u8>::from("XL");

        assert_eq!(word, LossyPauliWord::from(String::from("X_L")));
        assert_ne!(word, LossyPauliWord::from("ZL"));
        assert_ne!(word, LossyPauliWord::from("XI"));
        assert_ne!(word, LossyPauliWord::from("XLI"));
    }
}
