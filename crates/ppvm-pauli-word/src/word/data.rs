// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use bitvec::array::BitArray;
use core::panic;
use ppvm_traits::char::Pauli;
use ppvm_traits::traits::{HashFinalize, PauliIter, PauliStorage, PauliWordTrait};
use std::hash::{BuildHasher, Hash};
use std::ops::Index;

/// A fixed-width Pauli string stored as parallel `x` and `z` bit arrays.
///
/// Each qubit slot is encoded by two bits, `(x, z)`, giving the four
/// Pauli operators `I = (0,0)`, `X = (1,0)`, `Z = (0,1)`, `Y = (1,1)` —
/// the same encoding stabilizer-formalism tools use. The bit arrays live
/// in a single `PauliStorage` blob (typically `[u8; N]` or `[u64; N]`),
/// so a `PauliWord<[u8; 1]>` packs up to 8 qubits, `[u8; 2]` up to 16,
/// etc.
///
/// A cached hash speeds up use as a map key; pass `REHASH = false` if
/// you build words in bulk and want to control rehashing manually.
///
/// # Examples
///
/// ```
/// use ppvm_traits::char::Pauli;
/// use ppvm_traits::traits::PauliWordTrait;
/// use ppvm_pauli_word::word::PauliWord;
///
/// // Build a word from a string …
/// let w: PauliWord<[u8; 1]> = "XYZI".into();
/// assert_eq!(w.n_qubits(), 4);
/// assert_eq!(w.get(0), Pauli::X);
/// assert_eq!(w.get(3), Pauli::I);
///
/// // … or build it slot-by-slot.
/// let mut w2: PauliWord<[u8; 1]> = PauliWord::new(2);
/// w2.set(0, Pauli::X).set(1, Pauli::Y);
/// assert_eq!(w2.to_string(), "XY");
/// assert_eq!(w2.weight(), 2);
/// ```
#[derive(Debug, Clone)]
pub struct PauliWord<A: PauliStorage, S = fxhash::FxBuildHasher, const REHASH: bool = true> {
    /// X-bit array (one bit per qubit).
    pub xbits: BitArray<A>,
    /// Z-bit array (one bit per qubit).
    pub zbits: BitArray<A>,
    /// Number of qubits
    nqubits: usize,
    hash_cache: u64,
    _phantom: std::marker::PhantomData<S>,
}

impl<A: PauliStorage, S, const REHASH: bool> Hash for PauliWord<A, S, REHASH> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // self.xbits.data.hash(state);
        // self.zbits.data.hash(state);
        state.write_u64(self.hash_cache);
    }
}

impl<A: PauliStorage, S, const REHASH: bool> Eq for PauliWord<A, S, REHASH> {}

impl<A: PauliStorage, S, const REHASH: bool> PartialEq for PauliWord<A, S, REHASH> {
    fn eq(&self, other: &Self) -> bool {
        self.xbits.data == other.xbits.data && self.zbits.data == other.zbits.data
    }
}

impl<A: PauliStorage, S: Clone, const REHASH: bool> Copy for PauliWord<A, S, REHASH> {}

impl<A, S, const REHASH: bool> PauliIter for PauliWord<A, S, REHASH>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        PauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

impl<A, S, const REHASH: bool> PauliIter for &PauliWord<A, S, REHASH>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default + HashFinalize,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        PauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

// implement PauliString where A can be converted to chunks of u8, e.g u64
impl<A: PauliStorage, S: BuildHasher + Clone + Default + HashFinalize, const REHASH: bool>
    PauliWordTrait for PauliWord<A, S, REHASH>
{
    fn new(nqubits: usize) -> Self {
        Self {
            xbits: BitArray::ZERO,
            zbits: BitArray::ZERO,
            nqubits,
            hash_cache: 0,
            _phantom: std::marker::PhantomData,
        }
    }

    fn n_qubits(&self) -> usize {
        self.nqubits
    }

    #[inline(always)]
    fn loss_weight(&self) -> usize {
        0
    }

    #[inline(always)]
    fn get_xbit(&self, index: usize) -> bool {
        debug_assert!(index < self.nqubits, "Index out of bounds");
        self.xbits[index]
    }

    #[inline(always)]
    fn get_zbit(&self, index: usize) -> bool {
        debug_assert!(index < self.nqubits, "Index out of bounds");
        self.zbits[index]
    }

    #[inline(always)]
    fn get_lbit(&self, _index: usize) -> bool {
        false
    }

    #[inline(always)]
    fn set_xbit(&mut self, index: usize, value: bool) {
        debug_assert!(index < self.nqubits, "Index out of bounds");
        self.xbits.set(index, value);
    }
    #[inline(always)]
    fn set_zbit(&mut self, index: usize, value: bool) {
        debug_assert!(index < self.nqubits, "Index out of bounds");
        self.zbits.set(index, value);
    }

    #[inline]
    fn weight(&self) -> usize {
        let xs: &[u8] = bytemuck::bytes_of(&self.xbits.data);
        let zs: &[u8] = bytemuck::bytes_of(&self.zbits.data);
        debug_assert_eq!(xs.len(), zs.len());

        let mut total: u32 = 0;
        let (mut i, n) = (0usize, xs.len());

        // u64 chunks — one popcnt per chunk (x86 +popcnt; AArch64 CNT+ADDV).
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

    fn rehash(&mut self) {
        if REHASH {
            use std::hash::Hasher;
            let mut hasher = S::default().build_hasher();
            // Feed the raw storage bytes through the Hasher's byte-slice path
            // so FxHasher consumes 8 bytes per round instead of one byte at a
            // time. `A: bytemuck::Pod` (via the `PauliStorage` bound) makes
            // the `&[u8]` view safe — no padding, all bytes initialized.
            hasher.write(bytemuck::bytes_of(&self.xbits.data));
            hasher.write(bytemuck::bytes_of(&self.zbits.data));
            // Let the hasher finalize its own digest, told how wide the key is.
            // fxhash folds for narrow storage (its low bits correlate and
            // cluster hashbrown's buckets at high fill); strong hashers like
            // gxhash use the identity default. See `HashFinalize`.
            self.hash_cache = S::finalize_hash(hasher.finish(), std::mem::size_of::<A>());
        }
    }

    #[inline(always)]
    fn get(&self, index: usize) -> Pauli {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }
        match (self.xbits[index], self.zbits[index]) {
            (false, false) => Pauli::I,
            (false, true) => Pauli::Z,
            (true, false) => Pauli::X,
            (true, true) => Pauli::Y,
        }
    }

    #[inline(always)]
    fn is(&self, index: usize, pauli: Pauli) -> bool {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }
        match pauli {
            Pauli::I => !self.xbits[index] && !self.zbits[index],
            Pauli::X => self.xbits[index] && !self.zbits[index],
            Pauli::Z => !self.xbits[index] && self.zbits[index],
            Pauli::Y => self.xbits[index] && self.zbits[index],
            _ => false,
        }
    }

    #[inline(always)]
    fn set(&mut self, index: usize, pauli: Pauli) -> &mut Self {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }
        match pauli {
            Pauli::I => {
                self.xbits.set(index, false);
                self.zbits.set(index, false);
            }
            Pauli::X => {
                self.xbits.set(index, true);
                self.zbits.set(index, false);
            }
            Pauli::Z => {
                self.xbits.set(index, false);
                self.zbits.set(index, true);
            }
            Pauli::Y => {
                self.xbits.set(index, true);
                self.zbits.set(index, true);
            }
            _ => {
                panic!("Loss not supported in PauliWord! Use LossyPauliWord instead.");
            }
        }
        self.rehash();
        self
    }
}

impl<A: PauliStorage, S, const REHASH: bool> Ord for PauliWord<A, S, REHASH> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.nqubits
            .cmp(&other.nqubits)
            .then_with(|| self.xbits.cmp(&other.xbits))
            .then_with(|| self.zbits.cmp(&other.zbits))
    }
}

impl<A: PauliStorage, S, const REHASH: bool> PartialOrd for PauliWord<A, S, REHASH> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default + HashFinalize, const REHASH: bool>
    From<&str> for PauliWord<A, S, REHASH>
{
    fn from(value: &str) -> Self {
        PauliWord::from(value.to_string())
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default + HashFinalize, const REHASH: bool>
    From<String> for PauliWord<A, S, REHASH>
{
    fn from(value: String) -> Self {
        let n_qubits = value.chars().count();
        let chars = value.chars();
        let mut x = BitArray::ZERO;
        let mut z = BitArray::ZERO;

        let mut i = 0;
        for ch in chars {
            match ch {
                'I' => {
                    x.set(i, false);
                    z.set(i, false);
                }
                'X' => {
                    x.set(i, true);
                    z.set(i, false);
                }
                'Z' => {
                    x.set(i, false);
                    z.set(i, true);
                }
                'Y' => {
                    x.set(i, true);
                    z.set(i, true);
                }
                '_' => {
                    continue;
                }
                _ => panic!("Invalid Pauli character: {}", ch),
            };
            i += 1;
        }

        let mut ret = Self {
            xbits: x,
            zbits: z,
            nqubits: n_qubits,
            hash_cache: 0,
            _phantom: std::marker::PhantomData,
        };
        ret.rehash();
        ret
    }
}

impl<A: PauliStorage, S> From<PauliWord<A, S>> for usize {
    fn from(value: PauliWord<A, S>) -> Self {
        if value.nqubits > 64 {
            panic!("Cannot convert PauliString with more than 64 qubits to usize");
        }
        let mut result: BitArray<usize> = BitArray::ZERO;
        for i in 0..value.nqubits {
            result.set(2 * i, value.zbits[i]);
            result.set(2 * i + 1, value.xbits[i]);
        }
        result.into_inner()
    }
}

impl<A: PauliStorage, S, const REHASH: bool> Index<usize> for PauliWord<A, S, REHASH> {
    type Output = Pauli;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }

        match (self.xbits[index], self.zbits[index]) {
            (false, false) => &Pauli::I,
            (false, true) => &Pauli::Z,
            (true, false) => &Pauli::X,
            (true, true) => &Pauli::Y,
        }
    }
}

/// Iterator over the individual [`Pauli`] symbols of a [`PauliWord`].
pub struct PauliWordIter<'a, A: PauliStorage, S, const REHASH: bool = true> {
    word: &'a PauliWord<A, S, REHASH>,
    curr: usize,
}

impl<'a, A: PauliStorage, S: BuildHasher + Clone + Default + HashFinalize, const REHASH: bool>
    Iterator for PauliWordIter<'a, A, S, REHASH>
{
    type Item = Pauli;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr < self.word.nqubits {
            let pauli = self.word.get(self.curr);
            self.curr += 1;
            Some(pauli)
        } else {
            None
        }
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default + HashFinalize, const REHASH: bool>
    std::fmt::Display for PauliWord<A, S, REHASH>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.nqubits {
            let pauli = self.get(i);
            write!(f, "{}", pauli)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pauli_string_creation() {
        let ps = PauliWord::<[u64; 4]>::new(4);
        assert_eq!(ps.nqubits, 4);
    }
    #[test]
    fn test_pauli_string_set_get() {
        let ps = PauliWord::<[u64; 4]>::new(4);
        let ps = ps.set_new(0, Pauli::X);
        assert_eq!(ps.get(0), Pauli::X);
        let ps = ps.set_new(1, Pauli::Y);
        assert_eq!(ps.get(1), Pauli::Y);
        let ps = ps.set_new(2, Pauli::Z);
        assert_eq!(ps.get(2), Pauli::Z);
        let ps = ps.set_new(3, Pauli::I);
        assert_eq!(ps.get(3), Pauli::I);
    }
    #[test]
    fn test_pauli_string_display() {
        let ps = PauliWord::<[u64; 4]>::new(4);
        let ps = ps
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I);
        assert_eq!(ps.to_string(), "XYZI");
    }
    #[test]
    fn test_pauli_string_from_string() {
        let ps: PauliWord<[u64; 4]> = "XZYI".to_string().into();
        assert_eq!(ps.nqubits, 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Z);
        assert_eq!(ps.get(2), Pauli::Y);
        assert_eq!(ps.get(3), Pauli::I);
    }
    #[test]
    fn test_pauli_string_phase() {
        let ps: PauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: PauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: PauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: PauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
    }

    #[test]
    fn test_pauli_string_as_usize() {
        let ps: PauliWord<[u64; 4]> = "XZYY".to_string().into();
        let value: usize = ps.into();
        assert_eq!(value, 0b11110110); // X=01, Z=10, Y=11
    }
}
