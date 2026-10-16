use std::{fmt, marker::PhantomData, sync::atomic::{AtomicU64, Ordering}};

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
        self.xbits.set(i, v);
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
        self.invalidate_xz_hash();
    }

    #[inline]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits.set(i, v);
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss_hash();
        }
        self.invalidate_xz_hash();
    }

    #[inline]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits.set(i, x);
        self.zbits.set(i, z);
        if (x || z) && self.lbits[i] {
            self.lbits.set(i, false);
        }
        self.invalidate_xz_hash();
    }

    #[inline]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        self.xbits.set(i, xi);
        self.zbits.set(i, zi);
        self.xbits.set(j, xj);
        self.zbits.set(j, zj);
        if (xi || zi) && self.lbits[i] {
            self.lbits.set(i, false);
        }
        if (xj || zj) && self.lbits[j] {
            self.lbits.set(j, false);
        }
        self.invalidate_xz_hash();
    }

    #[inline]
    fn set_z_bit_pair(&mut self, i: usize, zi: bool, j: usize, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        self.zbits.set(i, zi);
        self.zbits.set(j, zj);
        if zi && self.lbits[i] {
            self.lbits.set(i, false);
        }
        if zj && self.lbits[j] {
            self.lbits.set(j, false);
        }
        self.invalidate_xz_hash();
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
        Self::from_planes(xbits, zbits, self.lbits, self.nqubits)
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
        Self::from_planes(xbits, zbits, self.lbits, self.nqubits)
    }
}

impl<A:PauliStorage, H> LossState for LossyPauliWord<A, H> {
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