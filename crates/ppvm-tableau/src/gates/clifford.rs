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

// Single-qubit gate on a `GeneralizedTableau`: skip lost qubits, delegate to
// the inner tableau's canonical (word-level) method.
macro_rules! impl_generalized_tableau_clifford {
    ($name:ident) => {
        fn $name(&mut self, index: usize) {
            if self.is_lost[index] {
                return;
            }
            self.tableau.$name(index);
        }
    };
}

// Two-qubit gate on a `GeneralizedTableau`: skip pairs with a lost qubit.
macro_rules! impl_generalized_tableau_clifford_pair {
    ($name:ident) => {
        fn $name(&mut self, control: usize, target: usize) {
            if self.is_lost[control] || self.is_lost[target] {
                return;
            }
            self.tableau.$name(control, target);
        }
    };
}

// Single source of truth for the per-gate Clifford phase/bit logic: every
// method below loops over the rows operating directly on the packed Pauli
// words (raw integer slices via `as_raw_slice`/`as_raw_mut_slice` plus a
// hoisted `index/bits`, `index%bits`, mask) rather than going through
// `bitvec`'s bounds-checked single-bit indexing inside the per-row loop.
//
// Every caller — a bare `Tableau`, a `GeneralizedTableau` (which delegates
// here via the `impl_generalized_tableau_clifford*` macros), and the fused
// batch path — runs through this one implementation, so there is no parallel
// copy that can silently diverge.
impl<T: Config> Clifford for Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    #[inline]
    fn x(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let zp = pw.word.zbits.data.as_raw_slice();
            let zw = zp[wi];
            pw.phase ^= (((zw & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn y(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            pw.phase ^= ((((xw ^ zw) & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn z(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let xw = xp[wi];
            pw.phase ^= (((xw & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn h(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= (((xb & zb) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn s(&mut self, index: usize) {
        // NOTE: S is the only clifford where forward and backward propagation
        // differ since it's non-hermitian; only the phase rule differs.
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            pw.phase ^= ((((xw & zw) & mask) != zero) as u8) << 1;
            zp[wi] = zw ^ (xw & mask);
        });
    }

    #[inline]
    fn cnot(&mut self, control: usize, target: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let (wc, sc) = (control / bits, control % bits);
        let (wt, st) = (target / bits, target % bits);
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xa = (xp[wc] >> sc) & one;
            let za = (zp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & zb) & (xb ^ za ^ one);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ (zb << sc);
            xp[wt] = xp[wt] ^ (xa << st);
        });
    }

    #[inline]
    fn cz(&mut self, control: usize, target: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let (wc, sc) = (control / bits, control % bits);
        let (wt, st) = (target / bits, target % bits);
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xa = (xp[wc] >> sc) & one;
            let za = (zp[wc] >> sc) & one;
            let xb = (xp[wt] >> st) & one;
            let zb = (zp[wt] >> st) & one;
            let phase_flip = (xa & xb) & (za ^ zb);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ (xb << sc);
            zp[wt] = zp[wt] ^ (xa << st);
        });
    }
}

impl<T: Config> CliffordExtensions for Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    // |    Gate    |  X  |  Y  |  Z  |
    // |:----------:|:---:|:---:|:---:|
    // |     s      |  Y  | -X  |  Z  |
    // |   s_dag    | -Y  |  X  |  Z  |
    // |   sqrt_x   |  X  |  Z  | -Y  |
    // | sqrt_x_dag |  X  | -Z  |  Y  |
    // |   sqrt_y   | -Z  |  Y  |  X  |
    // | sqrt_y_dag |  Z  |  Y  | -X  |

    #[inline]
    fn s_dag(&mut self, index: usize) {
        // NOTE: the backwards-prop version of S is just S†: same bit mapping,
        // phase rule differs (flip where x & !z).
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            pw.phase ^= ((((xw & !zw) & mask) != zero) as u8) << 1;
            zp[wi] = zw ^ (xw & mask);
        });
    }

    #[inline]
    fn sqrt_x(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            pw.phase ^= ((((zw & !xw) & mask) != zero) as u8) << 1;
            xp[wi] = xw ^ (zw & mask);
        });
    }

    #[inline]
    fn sqrt_x_dag(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            pw.phase ^= ((((xw & zw) & mask) != zero) as u8) << 1;
            xp[wi] = xw ^ (zw & mask);
        });
    }

    #[inline]
    fn sqrt_y(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= ((((xw & !zw) & mask) != zero) as u8) << 1;
        });
    }

    #[inline]
    fn sqrt_y_dag(&mut self, index: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let wi = index / bits;
        let off = index % bits;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let mask = one << off;
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xw = xp[wi];
            let zw = zp[wi];
            let xb = xw & mask;
            let zb = zw & mask;
            xp[wi] = (xw & !mask) | zb;
            zp[wi] = (zw & !mask) | xb;
            pw.phase ^= ((((zw & !xw) & mask) != zero) as u8) << 1;
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
    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        let bits = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let (wc, sc) = (control / bits, control % bits);
        let (wt, st) = (target / bits, target % bits);
        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let xc = (xp[wc] >> sc) & one;
            let zc = (zp[wc] >> sc) & one;
            let xt = (xp[wt] >> st) & one;
            let zt = (zp[wt] >> st) & one;
            let phase_flip = (xc & (xt ^ zt)) & (zc ^ zt ^ one);
            pw.phase ^= (((phase_flip & one) != zero) as u8) << 1;
            zp[wc] = zp[wc] ^ ((xt ^ zt) << sc);
            xp[wt] = xp[wt] ^ (xc << st);
            zp[wt] = zp[wt] ^ (xc << st);
        });
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Clifford for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    impl_generalized_tableau_clifford!(x);
    impl_generalized_tableau_clifford!(y);
    impl_generalized_tableau_clifford!(z);
    impl_generalized_tableau_clifford!(h);
    impl_generalized_tableau_clifford!(s);
    impl_generalized_tableau_clifford_pair!(cnot);
    impl_generalized_tableau_clifford_pair!(cz);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordExtensions
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    impl_generalized_tableau_clifford!(s_dag);
    impl_generalized_tableau_clifford!(sqrt_x);
    impl_generalized_tableau_clifford!(sqrt_x_dag);
    impl_generalized_tableau_clifford!(sqrt_y);
    impl_generalized_tableau_clifford!(sqrt_y_dag);
    impl_generalized_tableau_clifford_pair!(cy);
}

/// How many sites a per-word mask covers — one bit per distinct index, so this
/// equals the index count exactly when the indices are distinct.
#[inline]
fn mask_bits<S: PrimInt>(masks: &[S]) -> u64 {
    masks.iter().map(|m| u64::from(m.count_ones())).sum()
}

/// Take the fused per-word masks, or bail out of the caller: `return` when there
/// is nothing to do, and on a **repeated** index hand the whole batch to the
/// gate's outlined `$fallback`.
///
/// A mask carries one bit per site, so a repeated index would make the fused
/// sweep apply the gate once where the per-index loop — the `CliffordBatch`
/// contract, and what a legal `X 0 0` in a `.stim` file means — conjugates by
/// `Gᵏ`. `Gᵏ` is the identity only for the involutory gates (`S² = Z`,
/// `(√X)² = X`, `(√Y)² = Y`), so neither one bit nor an XOR-cancelled bit is
/// right for every family; only the per-index loop is (G-061).
///
/// The fast path costs exactly one branch on a flag `build_masks` has already
/// computed, and the slow path is a *call* — no scalar loop is inlined into the
/// fused body. Spelling the fallback out here instead (so each gate inlined a
/// whole `2n`-row sweep per index into its own preamble) cost 22–31% on the five
/// gates that map the x/z bit planes (`h`, `√X`, `√X†`, `√Y`, `√Y†`) when it was
/// measured on the `-2` crates: the bodies with real work to do lose registers
/// and I-cache to code that never runs. The phase-only gates (`x`/`y`/`z`/`s`)
/// never move, which is why A/B'ing only those rows hides it.
macro_rules! masks_or_scalar_fallback {
    ($self:ident, $indices:ident, $fallback:ident) => {
        match $self.build_masks($indices) {
            // Distinct indices: the fused sweep is licensed.
            Some((masks, n_words, false)) => (masks, n_words),
            // No rows or no indices: the batch is a no-op.
            None => return,
            // A repeated index: owe `Gᵏ`, i.e. the per-index loop.
            Some((.., true)) => return $self.$fallback($indices),
        }
    };
}

/// Define the outlined per-index loops the batched gates fall back to on a
/// repeated index (`$fallback` runs the canonical `$one` once per occurrence).
///
/// `#[cold]` + `#[inline(never)]` is the point: each loop is a full `2n`-row
/// sweep per index, and keeping it out of line keeps it off the fused bodies'
/// register budget and out of their cache lines — the caller keeps a call.
macro_rules! repeated_index_loops {
    ($($fallback:ident => $one:ident),+ $(,)?) => {
        impl<T: Config> Tableau<T>
        where
            <T::Storage as BitView>::Store: PrimInt,
        {
            $(
                #[cold]
                #[inline(never)]
                fn $fallback(&mut self, indices: &[usize]) {
                    for &q in indices {
                        self.$one(q);
                    }
                }
            )+
        }
    };
}

repeated_index_loops! {
    x_each => x,
    y_each => y,
    z_each => z,
    h_each => h,
    s_each => s,
    s_dag_each => s_dag,
    sqrt_x_each => sqrt_x,
    sqrt_x_dag_each => sqrt_x_dag,
    sqrt_y_each => sqrt_y,
    sqrt_y_dag_each => sqrt_y_dag,
}

impl<T: Config> Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    /// Build per-word bitmasks from a list of qubit indices.
    /// Returns `(masks, n_words, repeated)`. Stack-allocates for up to 8 storage
    /// words; spills to the heap beyond that, so there is no hard qubit cap.
    ///
    /// Distinctness is **checked**, not assumed, and checked here, in the third
    /// tuple slot: the mask holds one bit per *distinct* index, so it is short of
    /// `indices.len()` bits exactly when an index repeats. That is `n_words`
    /// popcounts (typically one or two) *after* the index walk, which is left
    /// byte-for-byte as it was — folding the test *into* the walk, as a
    /// `repeated |= word & bit` accumulator, measured worse on the `-2` crates
    /// (it adds a loop-carried dependency to the walk), and outlining this whole
    /// function was worse again. `repeated` sends the caller, via
    /// [`masks_or_scalar_fallback`], to its outlined per-index loop.
    #[inline]
    fn build_masks(&self, indices: &[usize]) -> Option<(MaskBuf<T>, usize, bool)> {
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
        let repeated = mask_bits(&masks) != indices.len() as u64;
        Some((masks, n_words, repeated))
    }
}

impl<T: Config> CliffordBatch for Tableau<T>
where
    <T::Storage as BitView>::Store: PrimInt,
{
    /// `X` is bit-preserving: phase flips for each masked qubit where z=1.
    #[inline]
    fn x_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, x_each);
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
    fn y_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, y_each);
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
    fn z_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, z_each);
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
    fn s_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, s_each);
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
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
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
    fn h_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, h_each);
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
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
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
    fn s_dag_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, s_dag_each);
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
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
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
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, sqrt_y_each);
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
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, sqrt_y_dag_each);
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
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, sqrt_x_each);
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
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        let (masks, n_words) = masks_or_scalar_fallback!(self, indices, sqrt_x_dag_each);
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
    impl_gen_tableau_batch_single!(x_many);
    impl_gen_tableau_batch_single!(y_many);
    impl_gen_tableau_batch_single!(z_many);
    impl_gen_tableau_batch_single!(h_many);
    impl_gen_tableau_batch_single!(s_many);
    impl_gen_tableau_batch_pair!(cnot_many);
    impl_gen_tableau_batch_pair!(cz_many);
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> CliffordExtensionsBatch
    for GeneralizedTableau<T, I, C>
where
    Complex<<T as Config>::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    impl_gen_tableau_batch_single!(s_dag_many);
    impl_gen_tableau_batch_single!(sqrt_x_many);
    impl_gen_tableau_batch_single!(sqrt_x_dag_many);
    impl_gen_tableau_batch_single!(sqrt_y_many);
    impl_gen_tableau_batch_single!(sqrt_y_dag_many);
    impl_gen_tableau_batch_pair!(cy_many);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_pauli_sum::config::fxhash::ByteF64;

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
    fn test_sqrt_x_dag_stabilizer() {
        // Z → +Y
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x_dag(0);
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
    fn test_sqrt_y_dag_stabilizer() {
        // Z → -X, X → +Z
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y_dag(0);
        let r = rows(&tab);
        assert_eq!(r[0], (false, true, 0), "destabilizer X should become +Z");
        assert_eq!(r[1], (true, false, 2), "stabilizer Z should become -X");
    }

    #[test]
    fn test_sqrt_x_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_x(0);
        tab.sqrt_x_dag(0);
        assert_eq!(rows(&tab), initial);
    }

    #[test]
    fn test_sqrt_y_round_trip() {
        let initial = rows(&GeneralizedTableau::new(1, 1e-12));
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.sqrt_y(0);
        tab.sqrt_y_dag(0);
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
        use ppvm_pauli_sum::config::fxhash::ByteF64;

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
            tab.sqrt_y_many(indices);
            snapshot(&tab)
        }

        #[test]
        fn test_sqrt_y_many_matches_individual() {
            let n = 8;
            let indices = vec![0, 2, 5, 7];
            assert_eq!(
                apply_individual_sqrt_y(n, &indices),
                apply_batch_sqrt_y(n, &indices)
            );
        }

        #[test]
        fn test_sqrt_y_dag_many_matches_individual() {
            let n = 8;
            let indices = vec![1, 3, 4, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(0);
            tab_ind.s(2);
            for &i in &indices {
                tab_ind.sqrt_y_dag(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(0);
            tab_batch.s(2);
            tab_batch.sqrt_y_dag_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_sqrt_x_many_matches_individual() {
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
            tab_batch.sqrt_x_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_sqrt_x_dag_many_matches_individual() {
            let n = 8;
            let indices = vec![2, 3, 5, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(1);
            tab_ind.s(4);
            for &i in &indices {
                tab_ind.sqrt_x_dag(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(1);
            tab_batch.s(4);
            tab_batch.sqrt_x_dag_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_h_many_matches_individual() {
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
            tab_batch.h_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cz_many_matches_individual() {
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
            tab_batch.cz_many(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cz_many_cross_word() {
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
            tab_batch.cz_many(&pairs);
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
            tab.sqrt_y_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.sqrt_x_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.h_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.x_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.y_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.z_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.s_many(&[]);
            assert_eq!(snapshot(&tab), initial);
            tab.s_dag_many(&[]);
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
            tab_batch.sqrt_y_many(&all);
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
            tab.sqrt_y_many(&indices);
            tab.sqrt_y_dag_many(&indices);
            assert_eq!(snapshot(&tab), initial);
        }

        #[test]
        fn test_x_many_matches_individual() {
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
            tab_batch.x_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_y_many_matches_individual() {
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
            tab_batch.y_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_z_many_matches_individual() {
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
            tab_batch.z_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_s_many_matches_individual() {
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
            tab_batch.s_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_s_dag_many_matches_individual() {
            let n = 8;
            let indices = vec![1, 3, 4, 6];
            let mut tab_ind = TTab::new(n);
            tab_ind.h(2);
            tab_ind.h(5);
            for &i in &indices {
                tab_ind.s_dag(i);
            }
            let mut tab_batch = TTab::new(n);
            tab_batch.h(2);
            tab_batch.h(5);
            tab_batch.s_dag_many(&indices);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cnot_many_matches_individual() {
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
            tab_batch.cnot_many(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cnot_many_cross_word() {
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
            tab_batch.cnot_many(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cy_many_matches_individual() {
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
            tab_batch.cy_many(&pairs);
            assert_eq!(snapshot(&tab_ind), snapshot(&tab_batch));
        }

        #[test]
        fn test_cy_many_cross_word() {
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
            tab_batch.cy_many(&pairs);
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
            tab_batch.h_many(&indices);

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
                tab.sqrt_x_many(&indices);
            }
            assert_eq!(snapshot(&tab), initial);
        }
    }

    // ---- Canonical (single-gate) path on a MULTI-WORD tableau ----
    //
    // The canonical `Clifford`/`CliffordExtensions` methods on `Tableau<T>`
    // now operate directly on the packed storage words (the former `_word`
    // implementations, folded in as the single source of truth). These tests
    // lock that path against the INDEPENDENT batch (`_many`) implementation —
    // which uses a wholly separate `build_masks` + popcount code path — on a
    // tableau whose `data` spans more than one storage word (qubit indices
    // >= the per-word bit width). Single-element batch calls are
    // behaviorally identical to one canonical call, so any divergence between
    // the two implementations (especially around the `index / bits`,
    // `index % bits`, mask, and shift arithmetic for high qubit indices)
    // surfaces here.
    mod multiword_canonical_tests {
        use super::*;
        use ppvm_pauli_sum::config::fxhash::ByteF64;

        // 2 u8 words → qubits 0..8 live in word 0, qubits 8..16 in word 1.
        type TC = ByteF64<2>;
        type TTab = Tableau<TC>;

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

        // Put the tableau in a non-trivial state straddling both words before
        // exercising the gate under test.
        fn seed(tab: &mut TTab) {
            for &q in &[0usize, 3, 7, 8, 11, 15] {
                tab.h(q);
            }
            tab.s(2);
            tab.s(10);
            tab.sqrt_y(5);
            tab.sqrt_x(13);
            tab.cnot(1, 9);
            tab.cz(6, 14);
        }

        // Indices chosen to hit word 0 (0..8), word 1 (8..16), and both word
        // boundaries (7 = last bit of word 0, 8 = first bit of word 1).
        const SINGLE_INDICES: [usize; 6] = [0, 7, 8, 11, 15, 4];

        /// For each index, apply the canonical single-qubit gate and the
        /// single-element batch call to freshly-seeded multi-word tableaus and
        /// assert the resulting storage words + phases match exactly.
        macro_rules! single_gate_matches_batch {
            ($test:ident, $gate:ident, $gate_many:ident) => {
                #[test]
                fn $test() {
                    let n = 16;
                    for &i in &SINGLE_INDICES {
                        let mut tab_canon = TTab::new(n);
                        seed(&mut tab_canon);
                        tab_canon.$gate(i);

                        let mut tab_batch = TTab::new(n);
                        seed(&mut tab_batch);
                        tab_batch.$gate_many(&[i]);

                        assert_eq!(
                            snapshot(&tab_canon),
                            snapshot(&tab_batch),
                            "canonical {} disagrees with batch {} at qubit {i}",
                            stringify!($gate),
                            stringify!($gate_many),
                        );
                    }
                }
            };
        }

        single_gate_matches_batch!(canonical_x_matches_batch_multiword, x, x_many);
        single_gate_matches_batch!(canonical_y_matches_batch_multiword, y, y_many);
        single_gate_matches_batch!(canonical_z_matches_batch_multiword, z, z_many);
        single_gate_matches_batch!(canonical_h_matches_batch_multiword, h, h_many);
        single_gate_matches_batch!(canonical_s_matches_batch_multiword, s, s_many);
        single_gate_matches_batch!(canonical_s_dag_matches_batch_multiword, s_dag, s_dag_many);
        single_gate_matches_batch!(
            canonical_sqrt_x_matches_batch_multiword,
            sqrt_x,
            sqrt_x_many
        );
        single_gate_matches_batch!(
            canonical_sqrt_x_dag_matches_batch_multiword,
            sqrt_x_dag,
            sqrt_x_dag_many
        );
        single_gate_matches_batch!(
            canonical_sqrt_y_matches_batch_multiword,
            sqrt_y,
            sqrt_y_many
        );
        single_gate_matches_batch!(
            canonical_sqrt_y_dag_matches_batch_multiword,
            sqrt_y_dag,
            sqrt_y_dag_many
        );

        // Two-qubit pairs covering: both in word 0, both in word 1, control
        // in word 0 / target in word 1, and the reverse.
        const PAIRS: [(usize, usize); 4] = [(0, 3), (9, 12), (2, 10), (13, 5)];

        /// Canonical two-qubit gate vs single-pair batch call, for same-word
        /// and cross-word pairs on a multi-word tableau.
        macro_rules! pair_gate_matches_batch {
            ($test:ident, $gate:ident, $gate_many:ident) => {
                #[test]
                fn $test() {
                    let n = 16;
                    for &(c, t) in &PAIRS {
                        let mut tab_canon = TTab::new(n);
                        seed(&mut tab_canon);
                        tab_canon.$gate(c, t);

                        let mut tab_batch = TTab::new(n);
                        seed(&mut tab_batch);
                        tab_batch.$gate_many(&[(c, t)]);

                        assert_eq!(
                            snapshot(&tab_canon),
                            snapshot(&tab_batch),
                            "canonical {} disagrees with batch {} on pair ({c}, {t})",
                            stringify!($gate),
                            stringify!($gate_many),
                        );
                    }
                }
            };
        }

        pair_gate_matches_batch!(canonical_cnot_matches_batch_multiword, cnot, cnot_many);
        pair_gate_matches_batch!(canonical_cz_matches_batch_multiword, cz, cz_many);
        pair_gate_matches_batch!(canonical_cy_matches_batch_multiword, cy, cy_many);

        /// Independent of the batch path: known forward-propagation rules for
        /// a single-qubit gate applied to a qubit living in the SECOND storage
        /// word. This checks the high-index mask and shift arithmetic against
        /// hand-derived expectations.
        ///
        /// Row layout (see `Tableau::new`): rows `0..n` are destabilizers
        /// `X_q` (row `q`), rows `n..2n` are stabilizers `Z_q` (row `n + q`).
        #[test]
        fn sqrt_x_on_second_word_qubit_known_transform() {
            // √X sends Z → -Y (x stays, z gains x, phase +2) and leaves X
            // unchanged.
            let n = 16;
            let q = 11; // word 1, bit 3
            let bit = 1u8 << (q % 8);
            let wi = q / 8; // 1

            let mut tab = TTab::new(n);
            tab.sqrt_x(q);

            // Destabilizer row (X_q): unchanged, no phase.
            let dz = &tab.data[q];
            assert_eq!(dz.word.xbits.data.as_raw_slice()[wi] & bit, bit);
            assert_eq!(dz.word.zbits.data.as_raw_slice()[wi] & bit, 0);
            assert_eq!(dz.phase, 0);

            // Stabilizer row (Z_q → -Y_q): x set, z set, phase +2.
            let st = &tab.data[n + q];
            assert_eq!(st.word.xbits.data.as_raw_slice()[wi] & bit, bit);
            assert_eq!(st.word.zbits.data.as_raw_slice()[wi] & bit, bit);
            assert_eq!(st.phase, 2);
        }

        /// CNOT with control and target in DIFFERENT storage words, checked
        /// against the known stabilizer-propagation truth table. Control = 2
        /// (word 0), target = 10 (word 1). Fresh tableau, so the relevant rows
        /// are the single-qubit destabilizers/stabilizers of qubits 2 and 10.
        ///
        /// Row layout (see `Tableau::new`): destabilizer `X_q` is row `q`,
        /// stabilizer `Z_q` is row `n + q`.
        #[test]
        fn cnot_cross_word_known_transform() {
            let n = 16;
            let (c, t) = (2usize, 10usize);
            let mut tab = TTab::new(n);
            tab.cnot(c, t);

            let cbit = 1u8 << (c % 8);
            let tbit = 1u8 << (t % 8);
            let (cw, tw) = (c / 8, t / 8);

            // Read (xc, zc, xt, zt, phase) for a given row index.
            let read = |row: usize| {
                let pw = &tab.data[row];
                let x = pw.word.xbits.data.as_raw_slice();
                let z = pw.word.zbits.data.as_raw_slice();
                (
                    x[cw] & cbit != 0,
                    z[cw] & cbit != 0,
                    x[tw] & tbit != 0,
                    z[tw] & tbit != 0,
                    pw.phase,
                )
            };

            // X_c (destabilizer of control) → X_c X_t, no phase.
            assert_eq!(read(c), (true, false, true, false, 0));
            // Z_c (stabilizer of control) → Z_c, unchanged.
            assert_eq!(read(n + c), (false, true, false, false, 0));
            // X_t (destabilizer of target) → X_t, unchanged.
            assert_eq!(read(t), (false, false, true, false, 0));
            // Z_t (stabilizer of target) → Z_c Z_t, no phase.
            assert_eq!(read(n + t), (false, true, false, true, 0));
        }
    }

    // G-061: a repeated qubit index in a batched Clifford owes conjugation by
    // `Gᵏ`, i.e. the per-index loop — `X 0 0` is legal Stim meaning
    // apply-to-each-target-in-order, and the fused mask carries one bit per
    // *distinct* site. `Gᵏ` is the identity only for the involutory gates
    // (`S² = Z`, `(√X)² = X`, `(√Y)² = Y`), so neither one bit nor an
    // XOR-cancelled bit is right for every family.
    mod duplicate_index_tests {
        use super::*;
        use ppvm_pauli_sum::config::fxhash::ByteF64;

        // 2 u8 words → qubits 0..8 in word 0, qubits 8..16 in word 1.
        type TC = ByteF64<2>;
        type TTab = Tableau<TC>;
        type GenTab = GeneralizedTableau<TC>;

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

        /// A frame with all four bit patterns (`I`/`X`/`Y`/`Z`) live in the
        /// rows, so a wrong number of gate applications cannot cancel out by
        /// accident.
        fn prepared(n: usize) -> TTab {
            let mut tab = TTab::new(n);
            for q in 0..n {
                tab.h(q);
                if q % 2 == 0 {
                    tab.s(q);
                }
            }
            for q in 1..n {
                tab.cnot(q - 1, q);
            }
            tab
        }

        /// The batch on `$indices` must equal `$one` applied once per
        /// occurrence.
        macro_rules! check {
            ($many:ident, $one:ident, $indices:expr) => {{
                let indices: &[usize] = $indices;
                let mut batched = prepared(16);
                let mut looped = prepared(16);
                batched.$many(indices);
                for &q in indices {
                    looped.$one(q);
                }
                assert_eq!(
                    snapshot(&batched),
                    snapshot(&looped),
                    "{} on repeated indices {:?} must equal the {} loop",
                    stringify!($many),
                    indices,
                    stringify!($one),
                );
            }};
        }

        macro_rules! check_all {
            ($indices:expr) => {{
                check!(x_many, x, $indices);
                check!(y_many, y, $indices);
                check!(z_many, z, $indices);
                check!(h_many, h, $indices);
                check!(s_many, s, $indices);
                check!(s_dag_many, s_dag, $indices);
                check!(sqrt_x_many, sqrt_x, $indices);
                check!(sqrt_x_dag_many, sqrt_x_dag, $indices);
                check!(sqrt_y_many, sqrt_y, $indices);
                check!(sqrt_y_dag_many, sqrt_y_dag, $indices);
            }};
        }

        #[test]
        fn duplicate_indices_batch_equals_per_index_loop() {
            check_all!(&[0, 0]);
            check_all!(&[1, 1, 1]);
            check_all!(&[0, 2, 2, 3]);
        }

        /// Independent of the loop: `G²` is a Clifford in its own right, so
        /// `g_many(&[q, q])` must equal conjugation by `G²`.
        #[test]
        fn duplicate_indices_equal_conjugation_by_g_squared() {
            // Involutory: G² = I, so the batch must leave the frame untouched.
            for (name, apply) in [
                ("x_many", TTab::x_many as fn(&mut TTab, &[usize])),
                ("y_many", TTab::y_many),
                ("z_many", TTab::z_many),
                ("h_many", TTab::h_many),
            ] {
                let mut batched = prepared(16);
                let untouched = prepared(16);
                apply(&mut batched, &[1, 1]);
                assert_eq!(
                    snapshot(&batched),
                    snapshot(&untouched),
                    "{name}(&[1, 1]) must be the identity (G² = I)"
                );
            }

            // S² = Z, (S†)² = Z, (√X)² = X, (√X†)² = X, (√Y)² = Y, (√Y†)² = Y.
            for (name, apply, square) in [
                (
                    "s_many",
                    TTab::s_many as fn(&mut TTab, &[usize]),
                    TTab::z as fn(&mut TTab, usize),
                ),
                ("s_dag_many", TTab::s_dag_many, TTab::z),
                ("sqrt_x_many", TTab::sqrt_x_many, TTab::x),
                ("sqrt_x_dag_many", TTab::sqrt_x_dag_many, TTab::x),
                ("sqrt_y_many", TTab::sqrt_y_many, TTab::y),
                ("sqrt_y_dag_many", TTab::sqrt_y_dag_many, TTab::y),
            ] {
                let mut batched = prepared(16);
                let mut squared = prepared(16);
                apply(&mut batched, &[1, 1]);
                square(&mut squared, 1);
                assert_eq!(
                    snapshot(&batched),
                    snapshot(&squared),
                    "{name}(&[1, 1]) must equal conjugation by G²"
                );
            }
        }

        /// The same defect across storage words: the duplicate sits in word 1
        /// while another index sits in word 0.
        #[test]
        fn duplicate_indices_batch_equals_loop_multiword() {
            check!(h_many, h, &[3, 10, 10]);
            check!(s_many, s, &[3, 10, 10]);
            check!(sqrt_y_many, sqrt_y, &[3, 10, 10]);
        }

        /// The `GeneralizedTableau` wrappers filter lost qubits and forward the
        /// rest, duplicates included, so they inherit the same contract.
        #[test]
        fn duplicate_indices_generalized_tableau_equals_loop() {
            let indices: &[usize] = &[0, 2, 2];
            for lost in [None, Some(0usize)] {
                let mut batched: GenTab = GeneralizedTableau::new(16, 1e-12);
                let mut looped: GenTab = GeneralizedTableau::new(16, 1e-12);
                for tab in [&mut batched, &mut looped] {
                    tab.tableau = prepared(16);
                    if let Some(q) = lost {
                        tab.is_lost[q] = true;
                    }
                }
                batched.sqrt_x_many(indices);
                for &q in indices {
                    looped.sqrt_x(q);
                }
                assert_eq!(
                    snapshot(&batched.tableau),
                    snapshot(&looped.tableau),
                    "sqrt_x_many(&[0, 2, 2]) with lost = {lost:?}"
                );
            }
        }

        /// The user-visible symptom from G-061's `.stim` reproduction:
        /// `X 0 0` on `|0…0⟩` is the identity, so the sampled bit must be 0.
        #[test]
        fn batched_x_twice_on_ground_state_measures_zero() {
            let mut tab: GenTab = GeneralizedTableau::new(1, 1e-12);
            tab.x_many(&[0, 0]);
            assert_eq!(tab.measure(0), Some(false));
        }
    }
}
