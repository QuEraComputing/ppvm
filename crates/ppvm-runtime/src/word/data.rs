use crate::char::Pauli;
use crate::traits::{PauliIter, PauliStorage, PauliWordTrait};
use bitvec::array::BitArray;
use core::panic;
use std::hash::{BuildHasher, Hash};
use std::ops::Index;

#[derive(Debug, Clone)]
pub struct PauliWord<A: PauliStorage, S = fxhash::FxBuildHasher> {
    pub xbits: BitArray<A>,
    pub zbits: BitArray<A>,
    /// Number of qubits
    nqubits: usize,
    hash_cache: u64,
    _phantom: std::marker::PhantomData<S>,
}

impl<A: PauliStorage, S> Hash for PauliWord<A, S> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // self.xbits.data.hash(state);
        // self.zbits.data.hash(state);
        state.write_u64(self.hash_cache);
    }
}

impl<A: PauliStorage, S> Eq for PauliWord<A, S> {}

impl<A: PauliStorage, S> PartialEq for PauliWord<A, S> {
    fn eq(&self, other: &Self) -> bool {
        self.xbits.data == other.xbits.data && self.zbits.data == other.zbits.data
    }
}

impl<A, S> PauliIter for PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        PauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

impl<A, S> PauliIter for &PauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        PauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

// implement PauliString where A can be converted to chunks of u8, e.g u64
impl<A: PauliStorage, S: BuildHasher + Clone + Default> PauliWordTrait for PauliWord<A, S> {
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

    fn weight(&self) -> usize {
        (0..self.nqubits)
            .filter(|&i| self.xbits[i] || self.zbits[i])
            .count()
    }

    fn rehash(&mut self) {
        use std::hash::Hasher;
        let mut hasher = S::default().build_hasher();
        self.xbits.data.hash(&mut hasher);
        self.zbits.data.hash(&mut hasher);
        self.hash_cache = hasher.finish();
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
    fn get_multiple<const Q: usize>(&self, indices: [usize; Q]) -> Self {
        let mut result = Self::new(Q);
        for (i, &idx) in indices.iter().enumerate() {
            result.set(i, self.get(idx));
        }
        result
    }

    #[inline(always)]
    fn set_multiple<const Q: usize, B: PauliStorage>(
        &mut self,
        indices: [usize; Q],
        values: &Self,
    ) {
        if values.nqubits != Q {
            panic!("Values must have the same number of qubits as indices");
        }
        for (i, &idx) in indices.iter().enumerate() {
            self.xbits.set(idx, values.xbits[i]);
            self.zbits.set(idx, values.zbits[i]);
        }
        self.rehash();
    }

    #[inline(always)]
    fn get_slice(&self, slice: std::ops::Range<usize>) -> Self {
        if slice.end > self.nqubits {
            panic!("Slice out of bounds");
        }
        let n_qubits = slice.len();
        let mut xbits = BitArray::ZERO;
        let mut zbits = BitArray::ZERO;
        for (i, idx) in slice.into_iter().enumerate() {
            xbits.set(i, self.xbits[idx]);
            zbits.set(i, self.zbits[idx]);
        }
        let mut ret = Self {
            xbits,
            zbits,
            nqubits: n_qubits,
            hash_cache: 0,
            _phantom: std::marker::PhantomData,
        };
        ret.rehash();
        ret
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

impl<A: PauliStorage, S> Ord for PauliWord<A, S> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.nqubits != other.nqubits {
            panic!("Cannot compare PauliStrings with different number of qubits");
        }
        self.xbits
            .cmp(&other.xbits)
            .then(self.zbits.cmp(&other.zbits))
    }
}

impl<A: PauliStorage, S> PartialOrd for PauliWord<A, S> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.nqubits != other.nqubits {
            return None;
        }
        Some(
            self.xbits
                .cmp(&other.xbits)
                .then(self.zbits.cmp(&other.zbits)),
        )
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default> From<&str> for PauliWord<A, S> {
    fn from(value: &str) -> Self {
        PauliWord::from(value.to_string())
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default> From<String> for PauliWord<A, S> {
    fn from(value: String) -> Self {
        let n_qubits = value.chars().count();
        let mut chars = value.chars();
        let mut x = BitArray::ZERO;
        let mut z = BitArray::ZERO;

        let mut i = 0;
        while let Some(ch) = chars.next() {
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

impl<A: PauliStorage, S> Index<usize> for PauliWord<A, S> {
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

pub struct PauliWordIter<'a, A: PauliStorage, S> {
    word: &'a PauliWord<A, S>,
    curr: usize,
}

impl<'a, A: PauliStorage, S: BuildHasher + Clone + Default> Iterator for PauliWordIter<'a, A, S> {
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

impl<A: PauliStorage, S: BuildHasher + Clone + Default> std::fmt::Display for PauliWord<A, S> {
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
