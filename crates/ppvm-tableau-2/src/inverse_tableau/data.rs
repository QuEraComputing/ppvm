// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The tableau's storage: rows, signs, and where rows sit in the arena.

use super::bitword::{
    BITWORD_BITS, BitWord, TOP_BIT, ZERO, bit_position, flip_bit, get_bit, mask_if, set_bits,
    set_bits_of,
};
use super::view::QubitMask;

/// The inverse tableau of an `n`-qubit Clifford.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InverseTableau {
    pub(super) n: usize,
    /// Words per row and plane: `n.div_ceil(BITWORD_BITS)`.
    pub(super) stride: usize,
    /// X and Z words of all `2n` rows, arranged by `offsets`.
    pub(super) words: Vec<BitWord>,
    /// Bit `r` is the sign of row `r` (set = −).
    pub(super) signs: Vec<BitWord>,
}

// NOTE: we could impl Word for RowRef to make it a true PauliWord if we ever need it
/// A borrowed row: its X words, Z words and sign.
#[derive(Clone, Copy, Debug)]
pub struct RowRef<'a> {
    pub x: &'a [BitWord],
    pub z: &'a [BitWord],
    pub negative: bool,
}

impl<'a> RowRef<'a> {
    /// X bit of the row at qubit `q`.
    pub fn x_bit(&self, q: usize) -> bool {
        get_bit(self.x, q)
    }

    /// Z bit of the row at qubit `q`.
    pub fn z_bit(&self, q: usize) -> bool {
        get_bit(self.z, q)
    }

    /// Qubits where the row has an X or Y, ascending.
    pub fn x_support(&self) -> impl Iterator<Item = usize> + 'a {
        set_bits(self.x)
    }

    /// Qubits outside `mask` where the row has an X or Y, ascending.
    pub fn x_outside(&self, mask: &'a QubitMask) -> impl Iterator<Item = usize> + 'a {
        set_bits_of(self.x.iter().zip(&mask.words).map(|(&x, &m)| x & !m))
    }

    /// Qubits in `mask` where the row has a Z or Y, ascending.
    pub fn z_inside(&self, mask: &'a QubitMask) -> impl Iterator<Item = usize> + 'a {
        set_bits_of(self.z.iter().zip(&mask.words).map(|(&z, &m)| z & m))
    }
}

impl InverseTableau {
    /// The identity Clifford, i.e. the state `|0…0⟩`: `ix_q = X_q`, `iz_q = Z_q`.
    pub fn new(n: usize) -> Self {
        let stride = n.div_ceil(BITWORD_BITS);
        let mut t = Self {
            n,
            stride,
            words: vec![ZERO; 4 * n * stride],
            signs: vec![ZERO; (2 * n).div_ceil(BITWORD_BITS)],
        };
        for q in 0..n {
            let (ix, iz) = (t.x_at(q), t.z_at(n + q));
            flip_bit(&mut t.words[ix..], q);
            flip_bit(&mut t.words[iz..], q);
        }
        t
    }

    pub fn num_qubits(&self) -> usize {
        self.n
    }

    /// Row `ix_q = C⁻¹ X_q C`.
    pub fn ix(&self, q: usize) -> RowRef<'_> {
        assert!(q < self.n, "qubit {q} out of range");
        self.row(q)
    }

    /// Row `iz_q = C⁻¹ Z_q C`.
    pub fn iz(&self, q: usize) -> RowRef<'_> {
        assert!(q < self.n, "qubit {q} out of range");
        self.row(self.n + q)
    }

    pub(super) fn row(&self, r: usize) -> RowRef<'_> {
        let (x, z) = (self.x_at(r), self.z_at(r));
        RowRef {
            x: &self.words[x..x + self.stride],
            z: &self.words[z..z + self.stride],
            negative: self.negative(r),
        }
    }

    /// Whether row `r` has sign `−`.
    #[inline(always)]
    pub(super) fn negative(&self, r: usize) -> bool {
        get_bit(&self.signs, r)
    }

    pub(super) fn flip_sign(&mut self, r: usize) {
        flip_bit(&mut self.signs, r);
    }

    // Arrangement: blocks by row kind, ix-X, ix-Z, iz-X, iz-Z (#204's order).
    // Row r is row i = r % n of kind r / n (0 for ix, 1 for iz).
    #[inline(always)]
    pub(super) fn offsets(n: usize, stride: usize, kind: usize, i: usize) -> (usize, usize) {
        (
            ((2 * kind) * n + i) * stride,
            ((2 * kind + 1) * n + i) * stride,
        )
    }

    /// The rows of `qubits` in generator order: `ix_q`, `iz_q` per qubit. `M`
    /// is twice the number of qubits.
    #[inline(always)]
    pub(super) fn qubit_rows<const M: usize>(&self, qubits: &[usize]) -> [usize; M] {
        let mut rows = [0; M];
        for (k, &q) in qubits.iter().enumerate() {
            rows[2 * k] = q;
            rows[2 * k + 1] = self.n + q;
        }
        rows
    }

    /// The X words and the Z words of the rows of distinct qubits, in the
    /// order of [`qubit_rows`](Self::qubit_rows), mutably.
    ///
    /// Panics unless every qubit is in range and no qubit repeats: distinct
    /// qubits own distinct rows, which makes the ranges disjoint, so this is
    /// the one place that checks it.
    #[inline(always)]
    pub(super) fn qubit_rows_mut<const M: usize>(
        &mut self,
        qubits: &[usize],
    ) -> ([&mut [BitWord]; M], [&mut [BitWord]; M]) {
        assert_eq!(M, 2 * qubits.len());
        for (i, &q) in qubits.iter().enumerate() {
            assert!(q < self.n, "qubit {q} out of range");
            assert!(!qubits[..i].contains(&q), "qubit {q} repeated");
        }
        let rows: [usize; M] = self.qubit_rows(qubits);
        let (xs, zs) = (rows.map(|r| self.x_at(r)), rows.map(|r| self.z_at(r)));
        let (stride, words) = (self.stride, self.words.as_mut_ptr());
        // SAFETY: distinct qubits in range own distinct rows, whose
        // `stride`-word ranges of `words` are disjoint; the slices borrow
        // `self` mutably for their lifetime.
        unsafe {
            (
                xs.map(|x| std::slice::from_raw_parts_mut(words.add(x), stride)),
                zs.map(|z| std::slice::from_raw_parts_mut(words.add(z), stride)),
            )
        }
    }

    /// Kind and index of row `r < 2n`, without dividing.
    #[inline(always)]
    fn kind_and_index(&self, r: usize) -> (usize, usize) {
        let kind = usize::from(r >= self.n);
        (kind, r - kind * self.n)
    }

    #[inline(always)]
    pub(super) fn x_at(&self, r: usize) -> usize {
        let (kind, i) = self.kind_and_index(r);
        Self::offsets(self.n, self.stride, kind, i).0
    }

    #[inline(always)]
    pub(super) fn z_at(&self, r: usize) -> usize {
        let (kind, i) = self.kind_and_index(r);
        Self::offsets(self.n, self.stride, kind, i).1
    }

    /// Visit every row in order with its X words and Z words; `f` returns
    /// whether the row's sign flips. Signs are updated a word (64 rows) at a
    /// time. Walks the blocks of `offsets` directly: rows of one kind are
    /// `stride` words apart.
    #[inline(always)]
    pub(super) fn for_each_row(
        &mut self,
        mut f: impl FnMut(&mut [BitWord], &mut [BitWord]) -> bool,
    ) {
        let (n, stride) = (self.n, self.stride);
        if n == 0 {
            return;
        }
        debug_assert_eq!(
            Self::offsets(n, stride, 1, 0),
            (2 * n * stride, 3 * n * stride)
        );
        let (ix, iz) = self.words.split_at_mut(2 * n * stride);
        let (ix_x, ix_z) = ix.split_at_mut(n * stride);
        let (iz_x, iz_z) = iz.split_at_mut(n * stride);
        let (mut flips, mut r) = (ZERO, 0);
        for (xs, zs) in [(ix_x, ix_z), (iz_x, iz_z)] {
            for (x, z) in xs.chunks_exact_mut(stride).zip(zs.chunks_exact_mut(stride)) {
                let (sw, sm) = bit_position(r);
                flips |= sm & mask_if(f(x, z));
                if sm == TOP_BIT {
                    self.signs[sw] ^= flips;
                    flips = ZERO;
                }
                r += 1;
            }
        }
        if !r.is_multiple_of(BITWORD_BITS) {
            self.signs[r / BITWORD_BITS] ^= flips;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::*;
    use super::*;

    #[test]
    fn identity_rows_are_single_paulis() {
        for n in SIZES {
            let t = InverseTableau::new(n);
            for q in 0..n {
                for j in 0..n {
                    assert_eq!(t.ix(q).x_bit(j), j == q);
                    assert!(!t.ix(q).z_bit(j));
                    assert!(!t.iz(q).x_bit(j));
                    assert_eq!(t.iz(q).z_bit(j), j == q);
                }
                assert!(!t.ix(q).negative && !t.iz(q).negative);
            }
        }
    }

    #[test]
    fn identity_has_canonical_commutation() {
        for n in SIZES {
            assert_canonical_commutation(&InverseTableau::new(n));
        }
    }
}
