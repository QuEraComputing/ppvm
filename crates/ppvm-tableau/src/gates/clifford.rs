// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::prelude::*;
use bitvec::view::BitView;
use bitvec::view::BitViewSized;
use num::complex::Complex;
use num::{One, PrimInt, Zero};
use smallvec::{SmallVec, smallvec};

/// Per-word bitmask buffer used by the batched Clifford gates.
/// Stack-allocates for up to 8 storage words; spills to heap beyond.
type MaskBuf<T> = SmallVec<[<<T as Config>::Storage as BitView>::Store; 8]>;

macro_rules! impl_tableau_clifford {
    ($name:ident, $($index:ident),*) => {
        #[inline]
        fn $name(&mut self, $($index: usize),*) {
            self.data.iter_mut().for_each(|pw| {
                pw.$name($($index),*);
            });
        }
    };
}

macro_rules! impl_generalized_tableau_clifford {
    ($name:ident, $index:ident) => {
        fn $name(&mut self, $index: usize) {
            if self.is_lost[$index] {
                return;
            }
            self.tableau.$name($index);
        }
    };
    ($name:ident, $index0:ident, $index1:ident) => {
        fn $name(&mut self, $index0: usize, $index1: usize) {
            if self.is_lost[$index0] || self.is_lost[$index1] {
                return;
            }
            self.tableau.$name($index0, $index1);
        }
    };
}

impl<T: Config> Clifford for Tableau<T> {
    impl_tableau_clifford!(x, index);
    impl_tableau_clifford!(y, index);
    impl_tableau_clifford!(z, index);
    impl_tableau_clifford!(h, index);
    impl_tableau_clifford!(cnot, control, target);
    impl_tableau_clifford!(cz, control, target);

    fn s(&mut self, index: usize) {
        // NOTE: S is the only clifford where forward and backward propagation differ
        // since it's non-hermitian
        // only difference is the phase though
        // TODO: just use the conjugate sdagger impl
        self.data.iter_mut().for_each(|pw| {
            let phase = (pw.word.xbits[index] & pw.word.zbits[index]) as u8;
            pw.word.s(index);
            pw.phase ^= phase << 1;
        });
    }
}

impl<T: Config> CliffordExtensions for Tableau<T> {
    // |    Gate    |  X  |  Y  |  Z  |
    // |:----------:|:---:|:---:|:---:|
    // |     s      |  Y  | -X  |  Z  |
    // |   s_adj    | -Y  |  X  |  Z  |
    // |   sqrt_x   |  X  |  Z  | -Y  |
    // | sqrt_x_adj |  X  | -Z  |  Y  |
    // |   sqrt_y   | -Z  |  Y  |  X  |
    // | sqrt_y_adj |  Z  |  Y  | -X  |

    fn s_adj(&mut self, addr0: usize) {
        // NOTE: the backwards prop version of S is just S_adj
        self.data.iter_mut().for_each(|pw| {
            pw.s(addr0);
        });
    }

    fn sqrt_x(&mut self, addr0: usize) {
        self.data.iter_mut().for_each(|pw| {
            let x = pw.word.xbits[addr0];
            let z = pw.word.zbits[addr0];
            pw.word.xbits.set(addr0, x ^ z);
            pw.phase ^= ((z & !x) as u8) << 1;
        });
    }

    fn sqrt_x_adj(&mut self, addr0: usize) {
        self.data.iter_mut().for_each(|pw| {
            let x = pw.word.xbits[addr0];
            let z = pw.word.zbits[addr0];
            pw.word.xbits.set(addr0, x ^ z);
            pw.phase ^= ((x & z) as u8) << 1;
        });
    }

    fn sqrt_y(&mut self, addr0: usize) {
        self.data.iter_mut().for_each(|pw| {
            let x = pw.word.xbits[addr0];
            let z = pw.word.zbits[addr0];
            pw.word.xbits.set(addr0, z);
            pw.word.zbits.set(addr0, x);
            pw.phase ^= ((x & !z) as u8) << 1;
        });
    }

    fn sqrt_y_adj(&mut self, addr0: usize) {
        self.data.iter_mut().for_each(|pw| {
            let x = pw.word.xbits[addr0];
            let z = pw.word.zbits[addr0];
            pw.word.xbits.set(addr0, z);
            pw.word.zbits.set(addr0, x);
            pw.phase ^= ((z & !x) as u8) << 1;
        });
    }

    // control: row, target: col
    // | CY  |  I  |  X  |  Y  |  Z  |
    // |:---:|:---:|:---:|:---:|:---:|
    // |  I  | II  | ZX  | IY  | ZZ  |
    // |  X  | XY  | -YZ | XI  | YX  |
    // |  Y  | YY  | XZ  | YI  | -XX |
    // |  Z  | ZI  | IX  | ZY  | IZ  |
    //
    // Bit transforms: xc'=xc, zc'=zc^xt^zt, xt'=xt^xc, zt'=zt^xc
    // Phase +2 when: xc & (xt ^ zt) & !(zc ^ zt)
    fn cy(&mut self, addr0: usize, addr1: usize) {
        self.data.iter_mut().for_each(|pw| {
            let xc = pw.word.xbits[addr0];
            let zc = pw.word.zbits[addr0];
            let xt = pw.word.xbits[addr1];
            let zt = pw.word.zbits[addr1];
            pw.word.zbits.set(addr0, zc ^ xt ^ zt);
            pw.word.xbits.set(addr1, xt ^ xc);
            pw.word.zbits.set(addr1, zt ^ xc);
            pw.phase ^= ((xc & (xt ^ zt) & !(zc ^ zt)) as u8) << 1;
        });
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Clifford for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_clifford!(x, index);
    impl_generalized_tableau_clifford!(y, index);
    impl_generalized_tableau_clifford!(z, index);
    impl_generalized_tableau_clifford!(h, index);
    impl_generalized_tableau_clifford!(s, index);
    impl_generalized_tableau_clifford!(cnot, control, target);
    impl_generalized_tableau_clifford!(cz, control, target);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordExtensions
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
{
    impl_generalized_tableau_clifford!(s_adj, addr0);
    impl_generalized_tableau_clifford!(sqrt_x, addr0);
    impl_generalized_tableau_clifford!(sqrt_x_adj, addr0);
    impl_generalized_tableau_clifford!(sqrt_y, addr0);
    impl_generalized_tableau_clifford!(sqrt_y_adj, addr0);
    impl_generalized_tableau_clifford!(cy, addr0, addr1);
}

impl<T: Config> Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    /// Build per-word bitmasks from a list of qubit indices.
    /// Returns `(masks, n_words)`. Stack-allocates for up to 8 storage words;
    /// spills to the heap beyond that, so there is no hard qubit cap.
    #[inline]
    fn build_masks(&self, indices: &[usize]) -> Option<(MaskBuf<T>, usize)> {
        if self.data.is_empty() || indices.is_empty() {
            return None;
        }
        let n_words = self.data[0].word.xbits.data.as_raw_slice().len();
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mut masks: MaskBuf<T> = smallvec![zero; n_words];
        for &addr0 in indices {
            masks[addr0 / bits_per_word] =
                masks[addr0 / bits_per_word] | (one << (addr0 % bits_per_word));
        }
        Some((masks, n_words))
    }
}

impl<T: Config> CliffordBatch for Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    /// `X` is bit-preserving: phase flips for each masked qubit where z=1.
    #[inline]
    fn x_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let zp = pw.word.zbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                popcount += (zp[wi] & mask).count_ones();
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// `Y` is bit-preserving: phase flips for each masked qubit where x⊕z=1.
    #[inline]
    fn y_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                popcount += ((xp[wi] ^ zp[wi]) & mask).count_ones();
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// `Z` is bit-preserving: phase flips for each masked qubit where x=1.
    #[inline]
    fn z_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                popcount += (xp[wi] & mask).count_ones();
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Forward `S`: phase flips where x&z=1, then z ^= x for masked qubits.
    #[inline]
    fn s_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                popcount += ((xw & zw) & mask).count_ones();
                zp[wi] = zw ^ (xw & mask);
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Apply CNOT to many pairs, operating on raw storage words to avoid
    /// per-bit `bitvec` addressing. Pairs are applied sequentially per row,
    /// so semantics match the per-pair `cnot` loop exactly.
    fn cnot_batch(&mut self, pairs: &[(usize, usize)]) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xa = (xp[wc] >> sc) & one;
                let za = (zp[wc] >> sc) & one;
                let xb = (xp[wt] >> st) & one;
                let zb = (zp[wt] >> st) & one;
                // +2 phase when x_a & z_b & !(x_b ^ z_a); 2+2 == 0 mod 4 so XOR-accumulate
                phase_flips = phase_flips ^ ((xa & zb) & (xb ^ za ^ one));
                // z_a ^= z_b ; x_b ^= x_a
                zp[wc] = zp[wc] ^ (zb << sc);
                xp[wt] = xp[wt] ^ (xa << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }

    /// Apply `H` to multiple qubits using combined bitmask operations.
    /// H swaps x<->z bits (same as sqrt_y) but with different phase:
    /// phase += 2 when x=1 & z=1 (Y goes to -Y).
    #[inline]
    fn h_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = x_bits & z_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply CZ to many pairs on raw storage words. CZ is symmetric and touches
    /// only z-bits; pairs are applied sequentially per row, so semantics match
    /// the per-pair `cz` loop exactly.
    fn cz_batch(&mut self, pairs: &[(usize, usize)]) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xa = (xp[wc] >> sc) & one;
                let za = (zp[wc] >> sc) & one;
                let xb = (xp[wt] >> st) & one;
                let zb = (zp[wt] >> st) & one;
                // +2 phase when x_a & x_b & (z_a ^ z_b)
                phase_flips = phase_flips ^ ((xa & xb) & (za ^ zb));
                // z_a ^= x_b ; z_b ^= x_a
                zp[wc] = zp[wc] ^ (xb << sc);
                zp[wt] = zp[wt] ^ (xa << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }
}

impl<T: Config> CliffordExtensionsBatch for Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    /// Backward `S` (i.e. `S†`): same bit mapping as `S`, phase rule differs.
    /// Phase flips where x&!z=1, then z ^= x for masked qubits.
    #[inline]
    fn s_adj_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let mut popcount = 0u32;
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                popcount += ((xw & !zw) & mask).count_ones();
                zp[wi] = zw ^ (xw & mask);
            }
            pw.phase ^= ((popcount & 1) as u8) << 1;
        });
    }

    /// Apply CY to many pairs on raw storage words. Pairs are applied
    /// sequentially per row, so semantics match the per-pair `cy` loop exactly.
    fn cy_batch(&mut self, pairs: &[(usize, usize)]) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let mut phase_flips = zero;
            for &(control, target) in pairs {
                let (wc, sc) = (control / bits, control % bits);
                let (wt, st) = (target / bits, target % bits);
                let xc = (xp[wc] >> sc) & one;
                let zc = (zp[wc] >> sc) & one;
                let xt = (xp[wt] >> st) & one;
                let zt = (zp[wt] >> st) & one;
                // +2 phase when x_c & (x_t ^ z_t) & !(z_c ^ z_t)
                phase_flips = phase_flips ^ ((xc & (xt ^ zt)) & (zc ^ zt ^ one));
                // z_c ^= x_t ^ z_t ; x_t ^= x_c ; z_t ^= x_c
                zp[wc] = zp[wc] ^ ((xt ^ zt) << sc);
                xp[wt] = xp[wt] ^ (xc << st);
                zp[wt] = zp[wt] ^ (xc << st);
            }
            pw.phase ^= (((phase_flips & one) != zero) as u8) << 1;
        });
    }

    /// Apply `√Y` to multiple qubits using combined bitmask operations.
    /// All qubits targeting the same word are merged into a single mask,
    /// reducing N individual operations to O(n_words) per row.
    #[inline]
    fn sqrt_y_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = x_bits & !z_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply `(√Y)†` to multiple qubits using combined bitmask operations.
    #[inline]
    fn sqrt_y_adj_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let not_mask = !mask;
                let xw = xp[wi];
                let zw = zp[wi];
                let x_bits = xw & mask;
                let z_bits = zw & mask;
                xp[wi] = (xw & not_mask) | z_bits;
                zp[wi] = (zw & not_mask) | x_bits;
                let phase_bits = z_bits & !x_bits;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
            }
        });
    }

    /// Apply `√X` to multiple qubits using combined bitmask operations.
    #[inline]
    fn sqrt_x_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                let phase_bits = (zw & !xw) & mask;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
                xp[wi] = xw ^ (zw & mask);
            }
        });
    }

    /// Apply `(√X)†` to multiple qubits using combined bitmask operations.
    #[inline]
    fn sqrt_x_adj_batch(&mut self, indices: &[usize]) {
        let (masks, n_words) = match self.build_masks(indices) {
            Some(m) => m,
            None => return,
        };
        let zero = <T::Storage as BitView>::Store::zero();

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            for wi in 0..n_words {
                let mask = masks[wi];
                if mask == zero {
                    continue;
                }
                let xw = xp[wi];
                let zw = zp[wi];
                let phase_bits = (xw & zw) & mask;
                pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;
                xp[wi] = xw ^ (zw & mask);
            }
        });
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    /// Fast path: check if any qubit in the slice is lost
    #[inline]
    fn any_lost_single(&self, indices: &[usize]) -> bool {
        indices.iter().any(|&i| self.is_lost[i])
    }

    /// Fast path: check if any qubit pair has a lost qubit
    #[inline]
    fn any_lost_pair(&self, pairs: &[(usize, usize)]) -> bool {
        pairs
            .iter()
            .any(|&(c, t)| self.is_lost[c] || self.is_lost[t])
    }
}

macro_rules! impl_gen_tableau_batch_single {
    ($name:ident) => {
        fn $name(&mut self, indices: &[usize]) {
            if !self.any_lost_single(indices) {
                self.tableau.$name(indices);
                return;
            }
            let filtered: Vec<usize> = indices
                .iter()
                .copied()
                .filter(|&i| !self.is_lost[i])
                .collect();
            self.tableau.$name(&filtered);
        }
    };
}

macro_rules! impl_gen_tableau_batch_pair {
    ($name:ident) => {
        fn $name(&mut self, pairs: &[(usize, usize)]) {
            if !self.any_lost_pair(pairs) {
                self.tableau.$name(pairs);
                return;
            }
            let filtered: Vec<(usize, usize)> = pairs
                .iter()
                .copied()
                .filter(|&(c, t)| !self.is_lost[c] && !self.is_lost[t])
                .collect();
            self.tableau.$name(&filtered);
        }
    };
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordBatch
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    impl_gen_tableau_batch_single!(x_batch);
    impl_gen_tableau_batch_single!(y_batch);
    impl_gen_tableau_batch_single!(z_batch);
    impl_gen_tableau_batch_single!(h_batch);
    impl_gen_tableau_batch_single!(s_batch);
    impl_gen_tableau_batch_pair!(cnot_batch);
    impl_gen_tableau_batch_pair!(cz_batch);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordExtensionsBatch
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    impl_gen_tableau_batch_single!(s_adj_batch);
    impl_gen_tableau_batch_single!(sqrt_x_batch);
    impl_gen_tableau_batch_single!(sqrt_x_adj_batch);
    impl_gen_tableau_batch_single!(sqrt_y_batch);
    impl_gen_tableau_batch_single!(sqrt_y_adj_batch);
    impl_gen_tableau_batch_pair!(cy_batch);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_runtime::config::fxhash::ByteF64;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;

    /// Returns (xbit, zbit, phase) for each tableau row: (destabilizer, stabilizer).
    fn rows(tab: &TestTableau) -> [(bool, bool, u8); 2] {
        [0, 1].map(|i| {
            let pw = &tab.tableau.data[i];
            (pw.word.xbits[0], pw.word.zbits[0], pw.phase)
        })
    }

    // Initial |0⟩: destabilizer = X (1,0,0), stabilizer = Z (0,1,0)

    #[test]
    fn test_sqrt_x_stabilizer() {
        // Z → -Y: forward prop √X P √X†
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x(0);
        let r = rows(&tab);
        assert_eq!(r[0], (true, false, 0), "destabilizer X should stay X");
        assert_eq!(r[1], (true, true, 2), "stabilizer Z should become -Y");
    }

    #[test]
    fn test_sqrt_x_adj_stabilizer() {
        // Z → +Y
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x_adj(0);
        let r = rows(&tab);
        assert_eq!(r[0], (true, false, 0), "destabilizer X should stay X");
        assert_eq!(r[1], (true, true, 0), "stabilizer Z should become +Y");
    }

    #[test]
    fn test_sqrt_y_stabilizer() {
        // Z → +X, X → -Z
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y(0);
        let r = rows(&tab);
        assert_eq!(r[0], (false, true, 2), "destabilizer X should become -Z");
        assert_eq!(r[1], (true, false, 0), "stabilizer Z should become +X");
    }

    #[test]
    fn test_sqrt_y_adj_stabilizer() {
        // Z → -X, X → +Z
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y_adj(0);
        let r = rows(&tab);
        assert_eq!(r[0], (false, true, 0), "destabilizer X should become +Z");
        assert_eq!(r[1], (true, false, 2), "stabilizer Z should become -X");
    }

    #[test]
    fn test_sqrt_x_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x(0);
        tab.sqrt_x_adj(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y(0);
        tab.sqrt_y_adj(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_x_fourth_power_is_identity() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        for _ in 0..4 {
            tab.sqrt_x(0);
        }
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_fourth_power_is_identity() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        for _ in 0..4 {
            tab.sqrt_y(0);
        }
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_x_on_lost_qubit_is_noop() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.is_lost[0] = true;
        tab.sqrt_x(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_on_lost_qubit_is_noop() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.is_lost[0] = true;
        tab.sqrt_y(0);
        assert_eq!(rows(&tab), initial);
    }

    /// Returns (x0, z0, x1, z1, phase) for each of the 4 tableau rows of a 2-qubit tableau.
    fn rows2(tab: &GeneralizedTableau<TestConfig>) -> [(bool, bool, bool, bool, u8); 4] {
        [0, 1, 2, 3].map(|i| {
            let pw = &tab.tableau.data[i];
            (
                pw.word.xbits[0],
                pw.word.zbits[0],
                pw.word.xbits[1],
                pw.word.zbits[1],
                pw.phase,
            )
        })
    }

    #[test]
    fn test_cy_stabilizers() {
        // CY (control=0, target=1) forward-propagates Paulis as CY P CY†.
        // From the truth table: xc'=xc, zc'=zc^xt^zt, xt'=xt^xc, zt'=zt^xc.
        // Phase +2 when xc & (xt^zt) & !(zc^zt), i.e. only for XX→-YZ and YZ→-XX.
        //   XI → +XY  (xt'=0^1=1, zt'=0^1=1; no phase since xt^zt=0)
        //   IX →  ZX  (zc'=0^1^0=1, xt'=1^0=1; no phase since xc=0)
        //   ZI →  ZI  (zc'=1^0^0=1; no phase)
        //   IZ →  ZZ  (zc'=0^0^1=1; no phase since xc=0)
        let mut tab: GeneralizedTableau<TestConfig> = GeneralizedTableau::new(2, 1e-12);
        tab.cy(0, 1);
        let r = rows2(&tab);
        assert_eq!(r[0], (true, false, true, true, 0), "XI should become +XY");
        assert_eq!(r[1], (false, true, true, false, 0), "IX should become ZX");
        assert_eq!(r[2], (false, true, false, false, 0), "ZI should stay ZI");
        assert_eq!(r[3], (false, true, false, true, 0), "IZ should become ZZ");
    }

    #[test]
    fn test_cy_round_trip() {
        // CY is self-inverse: CY² = I
        let initial = rows2(&GeneralizedTableau::new(2, 1e-12));
        let mut tab: GeneralizedTableau<TestConfig> = GeneralizedTableau::new(2, 1e-12);
        tab.cy(0, 1);
        tab.cy(0, 1);
        assert_eq!(rows2(&tab), initial);
    }

    // ---- Batch method tests ----

    mod batch_tests {
        use super::*;
        use ppvm_runtime::config::fxhash::ByteF64;

        type TC = ByteF64<2>; // 2 u8 words = up to 16 qubits
        type TTab = Tableau<TC>;

        /// Helper: extract all (xbits_raw, zbits_raw, phase) from a Tableau.
        fn snapshot(tab: &TTab) -> Vec<(Vec<u8>, Vec<u8>, u8)> {
            tab.data
                .iter()
                .map(|pw| {
                    (
                        pw.word.xbits.data.as_raw_slice().to_vec(),
                        pw.word.zbits.data.as_raw_slice().to_vec(),
                        pw.phase,
                    )
                })
                .collect()
        }

        /// Apply individual gate calls and return the resulting snapshot.
        fn apply_individual_sqrt_y(n: usize, indices: &[usize]) -> Vec<(Vec<u8>, Vec<u8>, u8)> {
            let mut tab = TTab::new(n);
            // Put tableau in a non-trivial state
            tab.h(0);
            tab.h(3);
            tab.s(1);
            for &i in indices {
                tab.sqrt_y(i);
            }
            snapshot(&tab)
        }

        fn apply_batch_sqrt_y(n: usize, indices: &[usize]) -> Vec<(Vec<u8>, Vec<u8>, u8)> {
            let mut tab = TTab::new(n);
            tab.h(0);
            tab.h(3);
            tab.s(1);
            tab.sqrt_y_batch(indices);
            snapshot(&tab)
        }

        #[test]
        fn test_sqrt_y_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 2, 5, 7];
            assert_eq!(
                apply_individual_sqrt_y(n, &indices),
                apply_batch_sqrt_y(n, &indices)
            );
        }

        #[test]
        fn test_sqrt_y_adj_batch_matches_individual() {
            let n = 8;
            let indices = vec![1, 3, 4, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.s(2);
            for &i in &indices {
                tab_ind.sqrt_y_adj(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.s(2);
            tab_batch.sqrt_y_adj_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_sqrt_x_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 1, 4, 7];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(2);
            tab_ind.s(5);
            for &i in &indices {
                tab_ind.sqrt_x(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(2);
            tab_batch.s(5);
            tab_batch.sqrt_x_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_sqrt_x_adj_batch_matches_individual() {
            let n = 8;
            let indices = vec![2, 3, 5, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(1);
            tab_ind.s(4);
            for &i in &indices {
                tab_ind.sqrt_x_adj(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(1);
            tab_batch.s(4);
            tab_batch.sqrt_x_adj_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_h_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 3, 5, 7];
            let mut tab_ind = TTab::new(n);
            tab_ind.s(1);
            tab_ind.sqrt_y(2);
            for &i in &indices {
                tab_ind.h(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.s(1);
            tab_batch.sqrt_y(2);
            tab_batch.h_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cz_batch_matches_individual() {
            let n = 8;
            let pairs = vec![(0, 1), (2, 3), (4, 5)];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.h(2);
            tab_ind.h(4);
            for &(c, t) in &pairs {
                tab_ind.cz(c, t);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.h(2);
            tab_batch.h(4);
            tab_batch.cz_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cz_batch_cross_word() {
            // Pairs straddle the two storage words (qubits 0..8 and 8..16).
            let n = 16;
            let pairs = vec![(1, 9), (8, 2), (7, 15), (10, 0)];
            let setup = [0usize, 7, 8, 9, 15];
            let mut tab_ind = TTab::new(n);
            for &q in &setup {
                tab_ind.h(q);
            }
            tab_ind.s(2);
            for &(c, t) in &pairs {
                tab_ind.cz(c, t);
            }
            let mut tab_batch = TTab::new(n);
            for &q in &setup {
                tab_batch.h(q);
            }
            tab_batch.s(2);
            tab_batch.cz_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_batch_empty_indices() {
            let n = 4;
            let initial = {
                let tab = TTab::new(n);
                snapshot(&tab)
            };
            let mut tab = TTab::new(n);
            tab.sqrt_y_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.sqrt_x_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.h_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.x_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.y_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.z_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.s_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.s_adj_batch(&[]);
            assert_eq!(snapshot(&tab), initial);
        }

        #[test]
        fn test_batch_all_qubits() {
            let n = 8;
            let all: Vec<usize> = (0..n).collect();
            let mut tab_ind = TTab::new(n);
            for &i in &all {
                tab_ind.sqrt_y(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.sqrt_y_batch(&all);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_batch_round_trip() {
            let n = 8;
            let indices = vec![1, 3, 5, 7];
            let initial = {
                let tab = TTab::new(n);
                snapshot(&tab)
            };
            let mut tab = TTab::new(n);
            tab.sqrt_y_batch(&indices);
            tab.sqrt_y_adj_batch(&indices);
            assert_eq!(snapshot(&tab), initial);
        }

        #[test]
        fn test_x_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 2, 5, 7];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.s(3);
            for &i in &indices {
                tab_ind.x(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.s(3);
            tab_batch.x_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_y_batch_matches_individual() {
            let n = 8;
            let indices = vec![1, 3, 4, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.s(2);
            for &i in &indices {
                tab_ind.y(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.s(2);
            tab_batch.y_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_z_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 1, 4, 7];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(2);
            tab_ind.s(5);
            for &i in &indices {
                tab_ind.z(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(2);
            tab_batch.s(5);
            tab_batch.z_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_s_batch_matches_individual() {
            let n = 8;
            let indices = vec![0, 2, 5, 7];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(1);
            tab_ind.h(4);
            for &i in &indices {
                tab_ind.s(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(1);
            tab_batch.h(4);
            tab_batch.s_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_s_adj_batch_matches_individual() {
            let n = 8;
            let indices = vec![1, 3, 4, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(2);
            tab_ind.h(5);
            for &i in &indices {
                tab_ind.s_adj(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(2);
            tab_batch.h(5);
            tab_batch.s_adj_batch(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cnot_batch_matches_individual() {
            let n = 8;
            let pairs = vec![(0, 1), (2, 3), (4, 5)];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.h(2);
            tab_ind.h(4);
            for &(c, t) in &pairs {
                tab_ind.cnot(c, t);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.h(2);
            tab_batch.h(4);
            tab_batch.cnot_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cnot_batch_cross_word() {
            // Pairs straddle the two storage words (qubits 0..8 and 8..16).
            let n = 16;
            let pairs = vec![(1, 9), (8, 2), (7, 15), (10, 0)];
            let setup = [0usize, 7, 8, 9, 15];
            let mut tab_ind = TTab::new(n);
            for &q in &setup {
                tab_ind.h(q);
            }
            tab_ind.s(2);
            for &(c, t) in &pairs {
                tab_ind.cnot(c, t);
            }
            let mut tab_batch = TTab::new(n);
            for &q in &setup {
                tab_batch.h(q);
            }
            tab_batch.s(2);
            tab_batch.cnot_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cy_batch_matches_individual() {
            let n = 8;
            let pairs = vec![(0, 1), (2, 3), (4, 5)];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.h(2);
            tab_ind.s(4);
            for &(c, t) in &pairs {
                tab_ind.cy(c, t);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.h(2);
            tab_batch.s(4);
            tab_batch.cy_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cy_batch_cross_word() {
            // Pairs straddle the two storage words (qubits 0..8 and 8..16).
            let n = 16;
            let pairs = vec![(1, 9), (8, 2), (7, 15), (10, 0)];
            let setup = [0usize, 7, 8, 9, 15];
            let mut tab_ind = TTab::new(n);
            for &q in &setup {
                tab_ind.h(q);
            }
            tab_ind.s(2);
            for &(c, t) in &pairs {
                tab_ind.cy(c, t);
            }
            let mut tab_batch = TTab::new(n);
            for &q in &setup {
                tab_batch.h(q);
            }
            tab_batch.s(2);
            tab_batch.cy_batch(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_batch_more_than_8_storage_words() {
            // 16 storage words exceeds the previous fixed-array limit of 8;
            // verifies the SmallVec heap-spill path in `build_masks`.
            type BigConfig = ByteF64<16>;
            type BigTab = Tableau<BigConfig>;

            fn big_snapshot(tab: &BigTab) -> Vec<(Vec<u8>, Vec<u8>, u8)> {
                tab.data
                    .iter()
                    .map(|pw| {
                        (
                            pw.word.xbits.data.as_raw_slice().to_vec(),
                            pw.word.zbits.data.as_raw_slice().to_vec(),
                            pw.phase,
                        )
                    })
                    .collect()
            }

            let n = 128;
            // Indices straddle word 0 (qubits 0..8), middle words, and the
            // last word (qubits 120..128).
            let indices: Vec<usize> = vec![0, 7, 8, 63, 64, 71, 72, 120, 127];

            let mut tab_ind = BigTab::new(n);
            tab_ind.h(3);
            tab_ind.s(60);
            for &i in &indices {
                tab_ind.h(i);
            }

            let mut tab_batch = BigTab::new(n);
            tab_batch.h(3);
            tab_batch.s(60);
            tab_batch.h_batch(&indices);

            assert_eq!(big_snapshot(&tab_ind), big_snapshot(&tab_batch));
        }

        #[test]
        fn test_batch_fourth_power_identity() {
            let n = 8;
            let indices = vec![0, 2, 4, 6];
            let initial = {
                let tab = TTab::new(n);
                snapshot(&tab)
            };
            let mut tab = TTab::new(n);
            for _ in 0..4 {
                tab.sqrt_x_batch(&indices);
            }
            assert_eq!(snapshot(&tab), initial);
        }
    }
}
