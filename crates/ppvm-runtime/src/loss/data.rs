use crate::char::Pauli;
use crate::traits::{PauliIter, PauliStorage, PauliWordTrait};
use bitvec::prelude::BitArray;
use std::hash::{BuildHasher, Hash};
use std::ops::Index;

/// A [`PauliWord`](crate::word::PauliWord)-like type that additionally
/// tracks per-qubit loss.
///
/// In neutral-atom hardware (and in any loss-aware error model), a
/// measurement can return "lost" instead of `0` or `1`. `LossyPauliWord`
/// adds a parallel `lbits` array marking which qubits have been lost; a
/// set loss bit excludes that qubit from the standard Pauli operations
/// and shows up as the [`Pauli::L`] symbol.
///
/// # Examples
///
/// ```
/// use ppvm_runtime::char::Pauli;
/// use ppvm_runtime::loss::LossyPauliWord;
/// use ppvm_runtime::traits::PauliWordTrait;
///
/// let w: LossyPauliWord<[u8; 1]> = "XLZL".into();
/// assert_eq!(w.weight(), 4);          // X, L, Z, L are all non-identity
/// assert_eq!(w.loss_weight(), 2);     // two qubits are lost
/// assert!(w.is(1, Pauli::L));
/// assert!(w.is(0, Pauli::X));
/// ```
#[derive(Debug, Clone)]
pub struct LossyPauliWord<A: PauliStorage, S = fxhash::FxBuildHasher> {
    /// X-bit array.
    pub xbits: BitArray<A>,
    /// Z-bit array.
    pub zbits: BitArray<A>,
    /// Loss bit array — `1` at index `i` means qubit `i` has been lost.
    pub lbits: BitArray<A>,
    /// Number of qubits
    nqubits: usize,
    hash_cache: u64,
    _phantom: std::marker::PhantomData<S>,
}

impl<A: PauliStorage, S> Hash for LossyPauliWord<A, S> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash_cache);
    }
}

impl<A: PauliStorage, S> Eq for LossyPauliWord<A, S> {}

impl<A: PauliStorage, S> PartialEq for LossyPauliWord<A, S> {
    fn eq(&self, other: &Self) -> bool {
        self.xbits.data == other.xbits.data
            && self.zbits.data == other.zbits.data
            && self.lbits.data == other.lbits.data
    }
}

impl<A, S> PauliIter for LossyPauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        LossyPauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

impl<A, S> PauliIter for &LossyPauliWord<A, S>
where
    A: PauliStorage,
    S: BuildHasher + Clone + Default,
{
    fn iter(&self) -> impl Iterator<Item = Pauli> {
        LossyPauliWordIter {
            word: self,
            curr: 0,
        }
    }
}

// implement PauliString where A can be converted to chunks of u8, e.g u64
impl<A: PauliStorage, S: BuildHasher + Clone + Default> PauliWordTrait for LossyPauliWord<A, S> {
    fn new(nqubits: usize) -> Self {
        Self {
            xbits: BitArray::ZERO,
            zbits: BitArray::ZERO,
            lbits: BitArray::ZERO,
            nqubits,
            hash_cache: 0,
            _phantom: std::marker::PhantomData,
        }
    }

    fn n_qubits(&self) -> usize {
        self.nqubits
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
    fn get_lbit(&self, index: usize) -> bool {
        debug_assert!(index < self.nqubits, "Index out of bounds");
        self.lbits[index]
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
        (self.xbits | self.zbits | self.lbits).count_ones()
    }

    fn loss_weight(&self) -> usize {
        self.lbits.count_ones()
    }

    fn rehash(&mut self) {
        use std::hash::Hasher;
        let mut hasher = S::default().build_hasher();
        self.xbits.data.hash(&mut hasher);
        self.zbits.data.hash(&mut hasher);
        self.lbits.data.hash(&mut hasher);
        self.hash_cache = hasher.finish();
    }

    #[inline(always)]
    fn get(&self, index: usize) -> Pauli {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }
        match (self.xbits[index], self.zbits[index], self.lbits[index]) {
            (false, false, false) => Pauli::I,
            (false, true, false) => Pauli::Z,
            (true, false, false) => Pauli::X,
            (true, true, false) => Pauli::Y,
            (false, false, true) => Pauli::L,
            _ => panic!(
                "Invalid Pauli character: LossyPauliWord cannot represent combinations of X, Y, Z with L"
            ),
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
            self.lbits.set(idx, values.lbits[i]);
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
        let mut lbits = BitArray::ZERO;
        for (i, idx) in slice.into_iter().enumerate() {
            xbits.set(i, self.xbits[idx]);
            zbits.set(i, self.zbits[idx]);
            lbits.set(i, self.lbits[idx]);
        }
        let mut ret = Self {
            xbits,
            zbits,
            lbits,
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
            Pauli::I => !self.xbits[index] && !self.zbits[index] && !self.lbits[index],
            Pauli::X => self.xbits[index] && !self.zbits[index] && !self.lbits[index],
            Pauli::Z => !self.xbits[index] && self.zbits[index] && !self.lbits[index],
            Pauli::Y => self.xbits[index] && self.zbits[index] && !self.lbits[index],
            Pauli::L => !self.xbits[index] && !self.zbits[index] && self.lbits[index],
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
                self.lbits.set(index, false);
            }
            Pauli::X => {
                self.xbits.set(index, true);
                self.zbits.set(index, false);
                self.lbits.set(index, false);
            }
            Pauli::Z => {
                self.xbits.set(index, false);
                self.zbits.set(index, true);
                self.lbits.set(index, false);
            }
            Pauli::Y => {
                self.xbits.set(index, true);
                self.zbits.set(index, true);
                self.lbits.set(index, false);
            }
            Pauli::L => {
                self.xbits.set(index, false);
                self.zbits.set(index, false);
                self.lbits.set(index, true);
            }
        }
        self.rehash();
        self
    }
}

impl<A: PauliStorage, S> Ord for LossyPauliWord<A, S> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.nqubits
            .cmp(&other.nqubits)
            .then_with(|| self.xbits.cmp(&other.xbits))
            .then_with(|| self.zbits.cmp(&other.zbits))
            .then_with(|| self.lbits.cmp(&other.lbits))
    }
}

impl<A: PauliStorage, S> PartialOrd for LossyPauliWord<A, S> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default> From<&str> for LossyPauliWord<A, S> {
    fn from(value: &str) -> Self {
        LossyPauliWord::from(value.to_string())
    }
}

impl<A: PauliStorage, S: BuildHasher + Clone + Default> From<String> for LossyPauliWord<A, S> {
    fn from(value: String) -> Self {
        let n_qubits = value.chars().count();
        let chars = value.chars();
        let mut x = BitArray::ZERO;
        let mut z = BitArray::ZERO;
        let mut l = BitArray::ZERO;

        let mut i = 0;
        for ch in chars {
            match ch {
                'I' => {
                    x.set(i, false);
                    z.set(i, false);
                    l.set(i, false);
                }
                'X' => {
                    x.set(i, true);
                    z.set(i, false);
                    l.set(i, false);
                }
                'Z' => {
                    x.set(i, false);
                    z.set(i, true);
                    l.set(i, false);
                }
                'Y' => {
                    x.set(i, true);
                    z.set(i, true);
                    l.set(i, false);
                }
                'L' => {
                    x.set(i, false);
                    z.set(i, false);
                    l.set(i, true);
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
            lbits: l,
            nqubits: n_qubits,
            hash_cache: 0,
            _phantom: std::marker::PhantomData,
        };
        ret.rehash();
        ret
    }
}

impl<A: PauliStorage, S> Index<usize> for LossyPauliWord<A, S> {
    type Output = Pauli;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.nqubits {
            panic!("Index out of bounds");
        }

        match (self.xbits[index], self.zbits[index], self.lbits[index]) {
            (false, false, false) => &Pauli::I,
            (false, true, false) => &Pauli::Z,
            (true, false, false) => &Pauli::X,
            (true, true, false) => &Pauli::Y,
            (false, false, true) => &Pauli::L,
            _ => panic!("Invalid Pauli configuration!"),
        }
    }
}

/// Iterator over the [`Pauli`] symbols of a [`LossyPauliWord`].
pub struct LossyPauliWordIter<'a, A: PauliStorage, S> {
    word: &'a LossyPauliWord<A, S>,
    curr: usize,
}

impl<'a, A: PauliStorage, S: BuildHasher + Clone + Default> Iterator
    for LossyPauliWordIter<'a, A, S>
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

impl<A: PauliStorage, S: BuildHasher + Clone + Default> std::fmt::Display for LossyPauliWord<A, S> {
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
        let ps = LossyPauliWord::<[u64; 4]>::new(4);
        assert_eq!(ps.nqubits, 4);
    }
    #[test]
    fn test_pauli_string_set_get() {
        let ps = LossyPauliWord::<[u64; 4]>::new(5);
        let ps = ps.set_new(0, Pauli::X);
        assert_eq!(ps.get(0), Pauli::X);
        let ps = ps.set_new(1, Pauli::Y);
        assert_eq!(ps.get(1), Pauli::Y);
        let ps = ps.set_new(2, Pauli::Z);
        assert_eq!(ps.get(2), Pauli::Z);
        let ps = ps.set_new(3, Pauli::I);
        assert_eq!(ps.get(3), Pauli::I);
        let ps = ps.set_new(4, Pauli::L);
        assert_eq!(ps.get(4), Pauli::L);
    }
    #[test]
    fn test_pauli_string_display() {
        let ps = LossyPauliWord::<[u64; 4]>::new(5);
        let ps = ps
            .set_new(0, Pauli::X)
            .set_new(1, Pauli::Y)
            .set_new(2, Pauli::Z)
            .set_new(3, Pauli::I)
            .set_new(4, Pauli::L);
        assert_eq!(ps.to_string(), "XYZIL");
    }
    #[test]
    fn test_pauli_string_from_string() {
        let ps: LossyPauliWord<[u64; 4]> = "XZYI".to_string().into();
        assert_eq!(ps.nqubits, 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Z);
        assert_eq!(ps.get(2), Pauli::Y);
        assert_eq!(ps.get(3), Pauli::I);

        let ps: LossyPauliWord<[u64; 4]> = "XLZL".to_string().into();
        assert_eq!(ps.nqubits, 4);
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::L);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::L);
    }
    #[test]
    fn test_pauli_string_phase() {
        let ps: LossyPauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: LossyPauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: LossyPauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
        let ps: LossyPauliWord<[u64; 4]> = "XYZI".to_string().into();
        assert_eq!(ps.get(0), Pauli::X);
        assert_eq!(ps.get(1), Pauli::Y);
        assert_eq!(ps.get(2), Pauli::Z);
        assert_eq!(ps.get(3), Pauli::I);
    }

    #[test]
    fn test_pauli_string_is_with_loss() {
        let ps: LossyPauliWord<[u64; 4]> = "XLIL".into();
        assert!(ps.is(0, Pauli::X));
        assert!(ps.is(1, Pauli::L));
        assert!(ps.is(2, Pauli::I));
        assert!(ps.is(3, Pauli::L));
        assert!(!ps.is(0, Pauli::L));
        assert!(!ps.is(1, Pauli::X));
    }

    #[test]
    fn test_pauli_string_weight_with_loss() {
        let ps: LossyPauliWord<[u64; 4]> = "IIII".into();
        assert_eq!(ps.weight(), 0);

        let ps: LossyPauliWord<[u64; 4]> = "XLIL".into();
        assert_eq!(ps.weight(), 3); // X, L, L are non-identity

        let ps: LossyPauliWord<[u64; 4]> = "LLLL".into();
        assert_eq!(ps.weight(), 4);
    }

    #[test]
    fn test_pauli_string_iter_with_loss() {
        let ps: LossyPauliWord<[u64; 4]> = "XLYZ".into();
        let paulis: Vec<Pauli> = ps.iter().collect();
        assert_eq!(paulis, vec![Pauli::X, Pauli::L, Pauli::Y, Pauli::Z]);
    }
}
