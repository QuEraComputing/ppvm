use std::{marker::PhantomData, sync::atomic::{AtomicU64, Ordering}};

use bitvec::array::BitArray;
use ppvm_pauli_word_2::{DefaultStorage, PauliStorage};
use ppvm_traits_2::{LossySite, Pauli, PauliBits, Word};

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
                self.xbits.set(qb, x);
                self.zbits.set(qb, z);
                self.invalidate_xz_hash();
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
    fn weight(&self) -> usize {
        let xs = bytemuck::bytes_of(&self.xbits.data);
        let zs = bytemuck::bytes_of(&self.zbits.data);
        let ls = bytemuck::bytes_of(&self.lbits.data);

        xs.iter()
            .zip(zs)
            .zip(ls)
            .map(|((&x, &z), &l)| (x | z | l).count_ones() as usize)
            .sum()
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
            xz_hash_cache: AtomicU64::new(self.loss_hash_cache.load(Ordering::Relaxed)),
            loss_hash_cache: AtomicU64::new(self.loss_hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}