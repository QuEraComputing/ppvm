// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{
    fmt,
    marker::PhantomData,
    sync::atomic::{AtomicBool, AtomicU64, Ordering},
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
    /// Whether the cached digest is valid, including a digest of zero.
    pub(crate) loss_hash_valid: AtomicBool,
    /// Lazy cache for the x/z structure digest.
    pub(crate) xz_hash_cache: AtomicU64,
    /// Whether the cached digest is valid, including a digest of zero.
    pub(crate) xz_hash_valid: AtomicBool,
    /// The private internal digest algorithm; `fn() -> H` keeps `Send + Sync`.
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> LossyPauliWord<A, H> {
    #[inline]
    fn check_qubit(&self, qubit: usize) {
        if qubit >= self.nqubits {
            panic!("qubit {qubit} out of bounds for {} sites", self.nqubits);
        }
    }

    #[inline]
    pub fn new(nqubits: usize) -> Self {
        Self::from_planes(BitArray::ZERO, BitArray::ZERO, BitArray::ZERO, nqubits)
    }

    #[inline]
    pub(crate) fn from_planes(
        xbits: BitArray<A>,
        zbits: BitArray<A>,
        lbits: BitArray<A>,
        nqubits: usize,
    ) -> Self {
        if nqubits > xbits.len() {
            panic!(
                "nqubits {nqubits} exceeds the {}-bit backing storage",
                xbits.len(),
            );
        }
        debug_assert!(
            lbits.iter_ones().all(|i| !xbits[i] && !zbits[i]),
            "lost sites must have both X and Z bits cleared",
        );
        Self {
            xbits,
            zbits,
            lbits,
            nqubits,
            loss_hash_cache: AtomicU64::new(0),
            loss_hash_valid: AtomicBool::new(false),
            xz_hash_cache: AtomicU64::new(0),
            xz_hash_valid: AtomicBool::new(false),
            _hasher: PhantomData,
        }
    }

    /// Mark the X/Z digest invalid after its bits change.
    #[inline]
    pub(crate) fn invalidate_xz_hash(&mut self) {
        *self.xz_hash_valid.get_mut() = false;
    }

    #[inline]
    pub(crate) fn invalidate_loss_hash(&mut self) {
        *self.loss_hash_valid.get_mut() = false;
    }

    /// check for whether a qubit was lost
    #[inline]
    pub fn is_lost(&self, qb: usize) -> bool {
        self.check_qubit(qb);
        self.lbits[qb]
    }

    #[inline]
    pub fn set_lost(&mut self, qb: usize) {
        self.check_qubit(qb);
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
        self.check_qubit(qb);
        if self.lbits[qb] {
            self.lbits.set(qb, false);
            self.invalidate_loss_hash();
        }
    }

    /// Build the reset-channel branch with `q` marked lost.
    pub fn with_lost(&self, qb: usize) -> Self {
        let mut branch = self.clone();
        branch.set_lost(qb);
        branch
    }

    /// Replace a site explicitly; `Present` restores a lost site.
    pub fn set(&mut self, qb: usize, site: LossySite<Pauli>) {
        self.check_qubit(qb);
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
        self.check_qubit(idx);
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
        let xz_hash_valid = self.xz_hash_valid.load(Ordering::Acquire);
        let loss_hash_valid = self.loss_hash_valid.load(Ordering::Acquire);
        Self {
            xbits: self.xbits,
            zbits: self.zbits,
            lbits: self.lbits,
            nqubits: self.nqubits,
            xz_hash_cache: AtomicU64::new(self.xz_hash_cache.load(Ordering::Relaxed)),
            xz_hash_valid: AtomicBool::new(xz_hash_valid),
            loss_hash_cache: AtomicU64::new(self.loss_hash_cache.load(Ordering::Relaxed)),
            loss_hash_valid: AtomicBool::new(loss_hash_valid),
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
        let nqubits = symbols.len();
        if nqubits > xbits.len() {
            panic!(
                "nqubits {nqubits} exceeds the {}-bit backing storage",
                xbits.len(),
            );
        }
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
        Self::from_planes(xbits, zbits, lbits, nqubits)
    }
}

impl<A: PauliStorage, H> From<String> for LossyPauliWord<A, H> {
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

/// Bit writes and toggles skip lost sites without changing their loss state.
impl<A: PauliStorage, H> PauliBits for LossyPauliWord<A, H> {
    #[inline]
    fn x_bit(&self, i: usize) -> bool {
        self.check_qubit(i);
        self.xbits[i]
    }

    #[inline]
    fn z_bit(&self, i: usize) -> bool {
        self.check_qubit(i);
        self.zbits[i]
    }

    #[inline]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        self.check_qubit(i);
        if self.lbits[i] {
            return;
        }
        if self.xbits[i] != v {
            self.xbits.set(i, v);
            self.invalidate_xz_hash();
        }
    }

    #[inline]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        self.check_qubit(i);
        if self.lbits[i] {
            return;
        }
        if self.zbits[i] != v {
            self.zbits.set(i, v);
            self.invalidate_xz_hash();
        }
    }

    #[inline]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        self.check_qubit(i);
        if self.lbits[i] {
            return;
        }
        if self.xbits[i] != x || self.zbits[i] != z {
            self.xbits.set(i, x);
            self.zbits.set(i, z);
            self.invalidate_xz_hash();
        }
    }

    #[inline]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        self.check_qubit(i);
        self.check_qubit(j);
        if i == j {
            self.set_xz_bits(j, xj, zj);
            return;
        }
        self.set_xz_bits(i, xi, zi);
        self.set_xz_bits(j, xj, zj);
    }

    #[inline]
    fn set_z_bit_pair(&mut self, i: usize, zi: bool, j: usize, zj: bool) {
        self.check_qubit(i);
        self.check_qubit(j);
        if i == j {
            self.set_z_bit(j, zj);
            return;
        }
        self.set_z_bit(i, zi);
        self.set_z_bit(j, zj);
    }

    #[inline]
    fn set_x_bit_and_z_bit(&mut self, x_i: usize, x: bool, z_i: usize, z: bool) {
        self.check_qubit(x_i);
        self.check_qubit(z_i);
        self.set_x_bit(x_i, x);
        self.set_z_bit(z_i, z);
    }

    #[inline]
    fn into_toggled_bits(mut self, i: usize, toggle: [bool; 2]) -> Self {
        self.check_qubit(i);
        if toggle[0] {
            self.set_x_bit(i, !self.x_bit(i));
        }
        if toggle[1] {
            self.set_z_bit(i, !self.z_bit(i));
        }
        self
    }

    #[inline]
    fn into_toggled_bits2(
        self,
        i: usize,
        toggle_i: [bool; 2],
        j: usize,
        toggle_j: [bool; 2],
    ) -> Self {
        self.check_qubit(i);
        self.check_qubit(j);
        if i == j {
            return self
                .into_toggled_bits(i, [toggle_i[0] ^ toggle_j[0], toggle_i[1] ^ toggle_j[1]]);
        }
        self.into_toggled_bits(i, toggle_i)
            .into_toggled_bits(j, toggle_j)
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
        let mut branch = self.clone();
        branch.clear_loss(i);
        branch
    }

    #[inline]
    fn loss_weight(&self) -> usize {
        LossyPauliWord::loss_weight(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_traits_2::Indexable;

    #[test]
    fn branch_helpers_invalidate_only_changed_components() {
        let masks = [[false, false], [true, false], [false, true], [true, true]];
        for input in ["IL", "LI", "LL", "XY"] {
            for warm_cache in [false, true] {
                let word = LossyPauliWord::<u8>::from(input);
                if warm_cache {
                    word.key_hash();
                }
                let mut branches = Vec::new();
                for i in 0..2 {
                    branches.push(word.with_lost(i));
                    branches.push(word.loss_cleared(i));
                    for mask_i in masks {
                        branches.push(word.toggled_bits(i, mask_i[0], mask_i[1]));
                        for j in 0..2 {
                            for mask_j in masks {
                                branches.push(word.toggled_bits2(i, mask_i, j, mask_j));
                                branches
                                    .push(word.clone().into_toggled_bits2(i, mask_i, j, mask_j));
                            }
                        }
                    }
                }
                for branch in branches {
                    let xz_unchanged = branch.xbits == word.xbits && branch.zbits == word.zbits;
                    let loss_unchanged = branch.lbits == word.lbits;
                    assert_eq!(
                        branch.xz_hash_valid.load(Ordering::Relaxed),
                        warm_cache && xz_unchanged,
                    );
                    assert_eq!(
                        branch.loss_hash_valid.load(Ordering::Relaxed),
                        warm_cache && loss_unchanged,
                    );
                    if warm_cache && xz_unchanged {
                        assert_eq!(
                            branch.xz_hash_cache.load(Ordering::Relaxed),
                            word.xz_hash_cache.load(Ordering::Relaxed),
                        );
                    }
                    if warm_cache && loss_unchanged {
                        assert_eq!(
                            branch.loss_hash_cache.load(Ordering::Relaxed),
                            word.loss_hash_cache.load(Ordering::Relaxed),
                        );
                    }
                    let fresh = LossyPauliWord::<u8>::from(branch.to_string());
                    assert_eq!(branch.key_hash(), fresh.key_hash());
                }
                assert_eq!(word.to_string(), input);
                assert_eq!(word.xz_hash_valid.load(Ordering::Relaxed), warm_cache);
                assert_eq!(word.loss_hash_valid.load(Ordering::Relaxed), warm_cache);
            }
        }
    }

    #[test]
    fn zero_hash_validity_survives_cloning_and_invalidates_independently() {
        let mut word = LossyPauliWord::<u8>::from("II");
        assert!(!word.xz_hash_valid.load(Ordering::Relaxed));
        assert!(!word.loss_hash_valid.load(Ordering::Relaxed));
        *word.xz_hash_valid.get_mut() = true;
        *word.loss_hash_valid.get_mut() = true;

        let mut cloned = word.clone();
        assert!(cloned.xz_hash_valid.load(Ordering::Relaxed));
        assert!(cloned.loss_hash_valid.load(Ordering::Relaxed));
        assert_eq!(cloned.xz_hash_cache.load(Ordering::Relaxed), 0);
        assert_eq!(cloned.loss_hash_cache.load(Ordering::Relaxed), 0);

        cloned.set_x_bit(0, true);
        assert!(!cloned.xz_hash_valid.load(Ordering::Relaxed));
        assert!(cloned.loss_hash_valid.load(Ordering::Relaxed));
        cloned.set_lost(1);
        assert!(!cloned.loss_hash_valid.load(Ordering::Relaxed));
        assert!(word.xz_hash_valid.load(Ordering::Relaxed));
        assert!(word.loss_hash_valid.load(Ordering::Relaxed));
    }

    #[test]
    fn bit_setters_skip_lost_sites_and_preserve_the_loss_cache() {
        let setters: [fn(&mut LossyPauliWord<u8>); 6] = [
            |w| w.set_x_bit(0, true),
            |w| w.set_z_bit(0, true),
            |w| w.set_xz_bits(0, true, true),
            |w| w.set_xz_bits2(0, true, true, 1, true, true),
            |w| w.set_z_bit_pair(0, true, 1, true),
            |w| w.set_x_bit_and_z_bit(0, true, 1, true),
        ];
        for setter in setters {
            let mut word = LossyPauliWord::<u8>::from("LL");
            *word.xz_hash_cache.get_mut() = 11;
            *word.xz_hash_valid.get_mut() = true;
            *word.loss_hash_cache.get_mut() = 22;
            *word.loss_hash_valid.get_mut() = true;

            setter(&mut word);

            assert_eq!(word, LossyPauliWord::from("LL"));
            assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
            assert!(word.xz_hash_valid.load(Ordering::Relaxed));
            assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
            assert!(word.loss_hash_valid.load(Ordering::Relaxed));
        }
        for input in ["LI", "IL"] {
            for setter in &setters[3..] {
                let mut word = LossyPauliWord::<u8>::from(input);
                *word.loss_hash_cache.get_mut() = 22;
                *word.loss_hash_valid.get_mut() = true;
                setter(&mut word);
                let lost = if input == "LI" { 0 } else { 1 };
                assert_eq!(word.get(lost), LossySite::Lost);
                assert!(!word.x_bit(lost) && !word.z_bit(lost));
                assert_ne!(word.get(1 - lost), LossySite::Present(Pauli::I));
                assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
                assert!(word.loss_hash_valid.load(Ordering::Relaxed));
            }
        }
    }

    #[test]
    fn owned_and_borrowed_toggles_preserve_loss_and_agree() {
        let masks = [[false, false], [true, false], [false, true], [true, true]];
        for input in ["LI", "IL", "LL", "XY"] {
            let word = LossyPauliWord::<u8>::from(input);
            for mask_i in masks {
                for mask_j in masks {
                    for j in [0, 1] {
                        let branch = word.toggled_bits2(0, mask_i, j, mask_j);
                        assert_eq!(
                            branch,
                            word.clone().into_toggled_bits2(0, mask_i, j, mask_j)
                        );
                        for q in 0..2 {
                            assert_eq!(branch.is_lost(q), word.is_lost(q));
                            if word.is_lost(q) {
                                assert!(!branch.x_bit(q) && !branch.z_bit(q));
                            }
                        }
                    }
                }
                assert_eq!(
                    word.toggled_bits(0, mask_i[0], mask_i[1]),
                    word.clone().into_toggled_bits(0, mask_i)
                );
            }
            assert_eq!(word.to_string(), input);
        }
        let word = LossyPauliWord::<u8>::from("LI");
        assert_eq!(
            word.toggled_bits2(0, [true, true], 1, [true, false])
                .to_string(),
            "LX"
        );
    }

    #[test]
    fn invalid_single_site_operations_panic_without_mutation() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let operations: [fn(&mut LossyPauliWord<u8>, usize); 14] = [
            |w, i| {
                w.get(i);
            },
            |w, i| {
                w.is_lost(i);
            },
            |w, i| {
                w.x_bit(i);
            },
            |w, i| {
                w.z_bit(i);
            },
            |w, i| w.set(i, LossySite::Present(Pauli::X)),
            |w, i| w.set_lost(i),
            |w, i| w.clear_loss(i),
            |w, i| w.set_x_bit(i, true),
            |w, i| w.set_z_bit(i, true),
            |w, i| w.set_xz_bits(i, true, true),
            |w, i| {
                w.with_lost(i);
            },
            |w, i| {
                w.loss_cleared(i);
            },
            |w, i| {
                w.toggled_bits(i, false, false);
            },
            |w, i| {
                w.clone().into_toggled_bits(i, [false, false]);
            },
        ];
        for input in ["", "I", "L"] {
            for index in [input.len(), 7, 8] {
                for operation in operations {
                    let mut word = LossyPauliWord::<u8>::from(input);
                    *word.xz_hash_cache.get_mut() = 11;
                    *word.xz_hash_valid.get_mut() = true;
                    *word.loss_hash_cache.get_mut() = 22;
                    *word.loss_hash_valid.get_mut() = true;
                    let original = word.clone();

                    assert!(
                        catch_unwind(AssertUnwindSafe(|| operation(&mut word, index))).is_err()
                    );
                    assert_eq!(word, original);
                    assert_eq!(word.weight(), original.weight());
                    assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
                    assert!(word.xz_hash_valid.load(Ordering::Relaxed));
                    assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
                    assert!(word.loss_hash_valid.load(Ordering::Relaxed));
                }
            }
        }
    }

    #[test]
    fn two_site_operations_check_both_indices_before_mutating() {
        use std::panic::{AssertUnwindSafe, catch_unwind};

        let operations: [fn(&mut LossyPauliWord<u8>, usize, usize); 5] = [
            |w, i, j| w.set_xz_bits2(i, true, true, j, true, true),
            |w, i, j| w.set_z_bit_pair(i, true, j, true),
            |w, i, j| w.set_x_bit_and_z_bit(i, true, j, true),
            |w, i, j| {
                w.toggled_bits2(i, [false, false], j, [false, false]);
            },
            |w, i, j| {
                w.clone()
                    .into_toggled_bits2(i, [false, false], j, [false, false]);
            },
        ];
        for (i, j) in [(0, 1), (1, 0), (0, 7), (7, 0)] {
            for operation in operations {
                let mut word = LossyPauliWord::<u8>::from("I");
                *word.xz_hash_cache.get_mut() = 11;
                *word.xz_hash_valid.get_mut() = true;
                *word.loss_hash_cache.get_mut() = 22;
                *word.loss_hash_valid.get_mut() = true;
                let original = word.clone();

                assert!(catch_unwind(AssertUnwindSafe(|| operation(&mut word, i, j))).is_err());
                assert_eq!(word, original);
                assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
                assert!(word.xz_hash_valid.load(Ordering::Relaxed));
                assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
                assert!(word.loss_hash_valid.load(Ordering::Relaxed));
            }
        }
    }

    #[test]
    fn parsing_rejects_words_exceeding_storage_capacity() {
        for input in ["IIIIIIIII", "XXXXXXXXI", "XXXXXXXXX", "IIII_IIIII"] {
            assert!(
                std::panic::catch_unwind(|| LossyPauliWord::<u8>::from(input)).is_err(),
                "accepted oversized input: {input}"
            );
            assert!(
                std::panic::catch_unwind(|| LossyPauliWord::<u8>::from(input.to_owned())).is_err()
            );
        }
    }

    #[test]
    #[should_panic(expected = "exceeds the 8-bit backing storage")]
    fn new_rejects_width_exceeding_storage_capacity() {
        LossyPauliWord::<u8>::new(9);
    }

    #[test]
    #[should_panic(expected = "exceeds the 8-bit backing storage")]
    fn from_planes_rejects_width_exceeding_storage_capacity() {
        LossyPauliWord::<u8>::from_planes(BitArray::ZERO, BitArray::ZERO, BitArray::ZERO, 9);
    }

    #[test]
    #[cfg(debug_assertions)]
    fn from_planes_rejects_pauli_bits_at_lost_sites() {
        for (x, z) in [(1u8, 0u8), (0, 1), (1, 1)] {
            assert!(
                std::panic::catch_unwind(|| {
                    LossyPauliWord::<u8>::from_planes(
                        BitArray::new(x),
                        BitArray::new(z),
                        BitArray::new(1),
                        1,
                    )
                })
                .is_err()
            );
        }
    }

    #[test]
    fn parsing_accepts_exact_capacity_and_ignores_separators() {
        for input in ["IIIIIIII", "XXXXXXXI", "_IIII_IIII_", "", "___"] {
            let word = LossyPauliWord::<u8>::from(input);
            let expected = input.replace('_', "");

            assert_eq!(word.n_sites(), expected.len());
            assert_eq!(word.to_string(), expected);
            assert_eq!(word.iter().count(), expected.len());
        }
    }

    #[test]
    fn unchanged_setters_preserve_both_caches() {
        let mut word = LossyPauliWord::<u8>::from("YZL");
        // Seed nonzero cache values to observe invalidation without hashing.
        *word.xz_hash_cache.get_mut() = 11;
        *word.xz_hash_valid.get_mut() = true;
        *word.loss_hash_cache.get_mut() = 22;
        *word.loss_hash_valid.get_mut() = true;

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
        assert!(word.xz_hash_valid.load(Ordering::Relaxed));
        assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
        assert!(word.loss_hash_valid.load(Ordering::Relaxed));
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
            let mut word = LossyPauliWord::<u8>::from("II");
            *word.xz_hash_cache.get_mut() = 11;
            *word.xz_hash_valid.get_mut() = true;
            *word.loss_hash_cache.get_mut() = 22;
            *word.loss_hash_valid.get_mut() = true;

            setter(&mut word);

            assert!(!word.xz_hash_valid.load(Ordering::Relaxed));
            assert_eq!(word.loss_hash_cache.load(Ordering::Relaxed), 22);
            assert!(word.loss_hash_valid.load(Ordering::Relaxed));
        }
        let mut word = LossyPauliWord::<u8>::from("L");
        *word.xz_hash_cache.get_mut() = 11;
        *word.xz_hash_valid.get_mut() = true;
        *word.loss_hash_cache.get_mut() = 22;
        *word.loss_hash_valid.get_mut() = true;
        word.set(0, LossySite::Present(Pauli::I));
        assert_eq!(word.to_string(), "I");
        assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
        assert!(word.xz_hash_valid.load(Ordering::Relaxed));
        assert!(!word.loss_hash_valid.load(Ordering::Relaxed));
    }

    #[test]
    fn repeated_indices_preserve_cache_when_final_bits_are_unchanged() {
        let mut word = LossyPauliWord::<u8>::from("I");
        *word.xz_hash_cache.get_mut() = 11;
        *word.xz_hash_valid.get_mut() = true;

        word.set_xz_bits2(0, true, true, 0, false, false);
        word.set_z_bit_pair(0, true, 0, false);

        assert_eq!(word.to_string(), "I");
        assert_eq!(word.xz_hash_cache.load(Ordering::Relaxed), 11);
        assert!(word.xz_hash_valid.load(Ordering::Relaxed));
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
    fn toggling_a_lost_site_preserves_loss() {
        let word = LossyPauliWord::<u8>::from("LL");

        assert_eq!(word.toggled_bits(0, true, false).to_string(), "LL");
        assert_eq!(word.toggled_bits(0, false, true).to_string(), "LL");
        assert_eq!(word.toggled_bits(0, true, true).to_string(), "LL");
        assert_eq!(word.toggled_bits(0, false, false), word);
        assert_eq!(word.to_string(), "LL");
    }

    #[test]
    fn toggling_two_sites_preserves_loss() {
        let word = LossyPauliWord::<u8>::from("LLL");

        let branch = word.toggled_bits2(0, [true, false], 2, [false, true]);
        assert_eq!(branch.to_string(), "LLL");
        assert_eq!(branch.loss_weight(), 3);

        let branch = word.toggled_bits2(0, [false, false], 2, [true, true]);
        assert_eq!(branch.to_string(), "LLL");
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
    fn setting_pauli_bits_preserves_lost_sites() {
        let mut word = LossyPauliWord::<u8>::from("LLL");

        word.set_x_bit(0, true);
        word.set_z_bit(1, true);
        word.set_xz_bits(2, true, true);

        assert_eq!(word.to_string(), "LLL");
        assert_eq!(word.loss_weight(), 3);
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
