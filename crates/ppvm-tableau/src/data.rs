// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::{fmt::Debug, marker::PhantomData};

// Only the rayon branch-coefficient helpers coalesce into a map now; the apply
// path and the default-feature branch path both use flat-Vec sort-merge/relabel.
#[cfg(feature = "rayon")]
use fxhash::FxHashMap as HashMap;

use bitvec::array::BitArray;
use bitvec::view::{BitView, BitViewSized};
use num::PrimInt;

use crate::prelude::*;
use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

type PhasedPauliWordNoHash<A, H> = PhasedPauliWord<A, H, PauliWord<A, H, false>>;

/// A `2n`-row stabilizer / destabilizer tableau.
///
/// Rows `0..n` hold the destabilizers; rows `n..2n` hold the
/// stabilizers. Each row is a [`PhasedPauliWord`] tracking both its
/// `X`/`Z` bits and a phase in `{±1, ±i}`. Implements every
/// Clifford-only operation natively (Hadamard, phase, CNOT, CZ, etc.).
///
/// # Examples
///
/// ```
/// use ppvm_pauli_sum::config::fxhash::ByteF64;
/// use ppvm_traits::traits::Clifford;
/// use ppvm_tableau::data::Tableau;
///
/// let mut tab: Tableau<ByteF64<1>> = Tableau::new(2);
/// tab.h(0);
/// tab.cnot(0, 1);
/// assert_eq!(tab.n_qubits, 2);
/// assert_eq!(tab.stabilizers().len(), 2);
/// ```
#[derive(Clone, Debug)]
pub struct Tableau<T: Config> {
    /// Number of qubits.
    pub n_qubits: usize,
    /// Destabilizer / Stabilizer tableau
    /// * Entries 0..n are the destabilizers
    /// * Entries n..2n are the stabilizers
    pub data: Vec<PhasedPauliWordNoHash<T::Storage, T::BuildHasher>>,
    pub(crate) rng: SmallRng,
}

impl<T: Config> Tableau<T> {
    fn new_data(n_qubits: usize) -> Vec<PhasedPauliWordNoHash<T::Storage, T::BuildHasher>> {
        // Initialize tableau for 0 state
        let mut data: Vec<PhasedPauliWordNoHash<T::Storage, T::BuildHasher>> =
            Vec::with_capacity(2 * n_qubits);
        let pw_cache = PhasedPauliWordNoHash::<T::Storage, T::BuildHasher>::new(n_qubits);
        for i in 0..n_qubits {
            // destabilizer
            let mut pw = pw_cache;
            pw.set(i, Pauli::X);
            data.push(pw);
        }
        for i in 0..n_qubits {
            // stabilizer
            let mut pw = pw_cache;
            pw.set(i, Pauli::Z);
            data.push(pw);
        }
        data
    }

    /// Construct a fresh tableau initialised to `|0…0⟩`.
    pub fn new(n_qubits: usize) -> Self {
        let data = Tableau::<T>::new_data(n_qubits);
        Self {
            n_qubits,
            data,
            rng: rand::make_rng(),
        }
    }

    /// Same as [`Tableau::new`], but seed the RNG deterministically.
    pub fn new_with_seed(n_qubits: usize, seed: u64) -> Self {
        let mut t = Self::new(n_qubits);
        t.rng = SmallRng::seed_from_u64(seed);
        t
    }

    pub fn reset_all(&mut self) {
        let data = Tableau::<T>::new_data(self.n_qubits);
        self.data = data;
    }

    /// View of the stabilizer rows (the upper half of the tableau).
    #[inline]
    pub fn stabilizers(&self) -> &[PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &self.data[self.n_qubits..]
    }

    /// Mutable view of the stabilizer rows.
    #[inline]
    pub fn stabilizers_mut(&mut self) -> &mut [PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &mut self.data[self.n_qubits..]
    }

    /// View of the destabilizer rows (the lower half of the tableau).
    #[inline]
    pub fn destabilizers(&self) -> &[PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &self.data[..self.n_qubits]
    }

    /// Mutable view of the destabilizer rows.
    #[inline]
    pub fn destabilizers_mut(
        &mut self,
    ) -> &mut [PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &mut self.data[..self.n_qubits]
    }

    // some helper functions for measurement impl
    pub(crate) fn find_z_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        // Find first stabilizer that anticommutes with Z_addr0
        self.stabilizers()
            .iter()
            .position(|stab| stab.word.anticommutes_at(addr0, (false, true)))
    }

    pub(crate) fn get_deterministic_outcome(&self, addr0: usize) -> bool
    where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        // find the outcome: either Z_addr0 or -Z_addr0 is a stabilizer
        // the stabilizer can be computed as the product of all stabilizers S_i
        // whose corresponding destabilizer D_i anticommutes with Z_addr0 (has X at addr0).
        // We have to actually multiply Paulis to also account for products of +i/-i
        let destabilizers = self.destabilizers();
        let stabilizers = self.stabilizers();
        let n = self.n_qubits;
        let mut result = PhasedPauliWordNoHash::<T::Storage, T::BuildHasher>::new(n);
        for (i, destab) in destabilizers.iter().enumerate() {
            if destab.word.xbits[addr0] {
                result *= &stabilizers[i];
            }
        }

        // phase >= 2 means -Z eigenvalue → outcome |1⟩ (true)
        debug_assert!(
            result.phase == 0 || result.phase == 2,
            "Measurement result cannot be imaginary!"
        );
        result.phase >= 2
    }

    pub(crate) fn update_tableau_according_to_outcome(
        &mut self,
        addr0: usize,
        q_idx: usize,
        outcome: bool,
    ) where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        let n = self.n_qubits;
        let (destabilizers, stabilizers) = self.data.split_at_mut(n);

        // Copy g_q once before the loop
        let g_q = stabilizers[q_idx];

        // Check if there are other stabilizers that anticommute with Z_addr0
        // If so, replace with g_j = g_j * g_q
        for i in 0..n {
            if i == q_idx {
                continue;
            }
            if stabilizers[i].word.xbits[addr0] {
                // Stabilizer i also anticommutes, so multiply by g_q to eliminate
                stabilizers[i] *= &g_q;
            }
            if destabilizers[i].word.xbits[addr0] {
                destabilizers[i] *= &g_q;
            }
        }

        // Update destabilizer q to be the old stabilizer q (before replacement)
        destabilizers[q_idx] = g_q;

        // Finally, replace g_q by ±Z
        let stab_q = &mut stabilizers[q_idx];
        stab_q.word.xbits = BitArray::ZERO;
        stab_q.word.zbits = BitArray::ZERO;
        stab_q.word.zbits.set(addr0, true);
        stab_q.phase = if outcome { 2 } else { 0 };
    }

    /// Apply CZ to N pairs with constant offset: (base+i, base+offset+i) for i in 0..count.
    /// All pairs must be in the same u64 word. This replaces N individual CZ calls
    /// with a single word-level shift+XOR operation per row.
    ///
    /// # Panics
    /// Debug-asserts that all bits are within the same word.
    #[inline]
    pub fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize)
    where
        <<T::Storage as BitView>::Store as TryFrom<usize>>::Error: Debug,
        <T::Storage as BitView>::Store: PrimInt + TryFrom<usize>,
    {
        if count == 0 {
            return;
        }
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let base_bit = base % bits_per_word;
        let word_idx = base / bits_per_word;

        debug_assert_eq!(
            (base + offset + count - 1) / bits_per_word,
            word_idx,
            "All CZ pairs must be in the same word"
        );

        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let count_mask = if count >= bits_per_word {
            !zero
        } else {
            (one << count) - one
        };
        let mask_c = count_mask << base_bit;
        let mask_t = count_mask << (base_bit + offset);

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();
            let x = xp[word_idx];
            let z = zp[word_idx];

            // Phase computation (must use original z before update)
            let xc = (x >> base_bit) & count_mask;
            let xt = (x >> (base_bit + offset)) & count_mask;
            let zc = (z >> base_bit) & count_mask;
            let zt = (z >> (base_bit + offset)) & count_mask;
            let phase_bits = xc & xt & (zc ^ zt);
            pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;

            // Z update: z[c] ^= x[t], z[t] ^= x[c]
            let z_delta = ((x >> offset) & mask_c) | ((x << offset) & mask_t);
            zp[word_idx] = z ^ z_delta;
        });
    }

    /// Apply CZ to N pairs with constant offset across two different words.
    /// Controls at (word_c, base_bit_c+i) and targets at (word_t, base_bit_t+i) for i in 0..count.
    /// word_c and word_t must be different.
    #[inline]
    pub fn cz_block_pairs_cross_word(
        &mut self,
        word_c: usize,
        base_bit_c: usize,
        word_t: usize,
        base_bit_t: usize,
        count: usize,
    ) where
        <T::Storage as BitView>::Store: PrimInt,
    {
        if count == 0 {
            return;
        }
        let one = <T::Storage as BitView>::Store::one();
        let zero = <T::Storage as BitView>::Store::zero();
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;

        debug_assert!(base_bit_c + count <= bits_per_word);
        debug_assert!(base_bit_t + count <= bits_per_word);
        debug_assert_ne!(word_c, word_t);

        let count_mask = if count >= bits_per_word {
            !zero
        } else {
            (one << count) - one
        };

        self.data.iter_mut().for_each(|pw| {
            let xp = pw.word.xbits.data.as_raw_mut_slice();
            let zp = pw.word.zbits.data.as_raw_mut_slice();

            // Extract aligned bits (shifted to 0..count-1)
            let xc = (xp[word_c] >> base_bit_c) & count_mask;
            let xt = (xp[word_t] >> base_bit_t) & count_mask;
            let zc = (zp[word_c] >> base_bit_c) & count_mask;
            let zt = (zp[word_t] >> base_bit_t) & count_mask;

            // Phase: x[c] & x[t] & (z[c] ^ z[t]) per pair
            let phase_bits = xc & xt & (zc ^ zt);
            pw.phase ^= ((phase_bits.count_ones() & 1) as u8) << 1;

            // z[c] ^= x[t]: place target x-bits at control positions
            zp[word_c] = zp[word_c] ^ (xt << base_bit_c);
            // z[t] ^= x[c]: place control x-bits at target positions
            zp[word_t] = zp[word_t] ^ (xc << base_bit_t);
        });
    }
}

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },  // +1
    Complex64 { re: 0.0, im: 1.0 },  // +i
    Complex64 { re: -1.0, im: 0.0 }, // -1
    Complex64 { re: 0.0, im: -1.0 }, // -i
];

/// Symplectic inner product of two tableau index values — the count of
/// shared set bits, used in stabilizer phase calculations.
#[inline]
pub fn symplectic_inner<I>(alpha: I, beta: I) -> u32
where
    I: TableauIndex,
{
    (alpha & beta).count_ones()
}

/// Compute the phase contribution from destabilizer anticommutation and odd-phase masks.
/// This is a pure function of its arguments (no self access needed), extracted to enable
/// use in parallel contexts where borrowing self is not possible.
#[inline]
pub(crate) fn compute_phase_with_mask_static<I: TableauIndex>(
    destab_anticomm_bits: I,
    basis_index: I,
    stab_anticomm_bits: I,
    odd_phase_mask: I,
) -> u8 {
    let mut phase = (2 * symplectic_inner(destab_anticomm_bits, basis_index) as u8) % 4;
    let active = basis_index & stab_anticomm_bits;
    let parity = (active & odd_phase_mask).count_ones() % 2;
    phase = (phase + 2 * parity as u8) % 4;
    phase
}

/// Minimum number of coefficients before engaging rayon parallelism.
/// Below this threshold, the sequential path is always used even with the `rayon` feature.
/// Benchmarked: at 8K coefficients rayon has ~24% overhead; at 32K it's 35% faster.
/// Set to 16384 to avoid regressions while capturing the large-coefficient wins.
#[cfg(feature = "rayon")]
pub(crate) const RAYON_COEFF_THRESHOLD: usize = 16384;

/// Parameters describing how a Pauli branch splits each coefficient.
///
/// Bundles the four anticommutation/phase-mask fields produced by
/// [`compute_decomposition`](GeneralizedTableau::compute_decomposition) with the
/// two scaling factors applied to the branch and non-branch contributions. They
/// always travel together, so passing them as one value keeps the coefficient
/// accumulation helpers under clippy's argument-count limit without an
/// `#[allow]`.
///
/// Only the rayon coefficient helpers consume this; the default-feature build
/// uses the inlined sort-merge path in
/// [`branch_with_coefficients`](GeneralizedTableau::branch_with_coefficients).
#[cfg(feature = "rayon")]
#[derive(Clone, Copy)]
pub(crate) struct BranchParams<I, CoeffType> {
    pub stab_anticomm_bits: I,
    pub destab_anticomm_bits: I,
    pub odd_phase_mask: I,
    pub phase_decomp: u8,
    pub coefficient_factor: Complex<CoeffType>,
    pub branch_factor: Complex<CoeffType>,
}

/// Sequential accumulation of branch coefficients.
#[cfg(feature = "rayon")]
fn branch_coefficients_seq<I, CoeffType>(
    items: impl IntoIterator<Item = (Complex<CoeffType>, I)>,
    capacity: usize,
    params: BranchParams<I, CoeffType>,
) -> HashMap<I, Complex<CoeffType>>
where
    I: TableauIndex,
    CoeffType: One + Zero + Clone + num::Num,
    Complex<CoeffType>:
        std::ops::Mul<Output = Complex<CoeffType>> + std::ops::AddAssign + From<Complex64> + Copy,
{
    let BranchParams {
        stab_anticomm_bits,
        destab_anticomm_bits,
        odd_phase_mask,
        phase_decomp,
        coefficient_factor,
        branch_factor,
    } = params;
    let mut map: HashMap<I, Complex<CoeffType>> =
        HashMap::with_capacity_and_hasher(capacity, Default::default());
    for (coeff, idx) in items {
        debug_assert!(
            !(coeff.re == CoeffType::zero() && coeff.im == CoeffType::zero()),
            "Coefficient should not be zero"
        );
        let branch_index = idx ^ stab_anticomm_bits;
        let branch_phase_contribution = compute_phase_with_mask_static(
            destab_anticomm_bits,
            idx,
            stab_anticomm_bits,
            odd_phase_mask,
        );
        let branch_phase = (branch_phase_contribution + phase_decomp) % 4;
        let phase_factor: Complex<CoeffType> =
            COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();
        let branch_coefficient = phase_factor * coeff * branch_factor;
        let nonbranch_coefficient = coeff * coefficient_factor;
        *map.entry(branch_index).or_insert(Complex::zero()) += branch_coefficient;
        *map.entry(idx).or_insert(Complex::zero()) += nonbranch_coefficient;
    }
    map
}

/// Accumulate branch coefficients. When the coefficient count exceeds
/// `RAYON_COEFF_THRESHOLD`, uses parallel map/collect into a Vec followed
/// by sequential accumulation. Below the threshold, falls back to sequential.
#[cfg(feature = "rayon")]
fn branch_coefficients_parallel<I, CoeffType>(
    items: &[(Complex<CoeffType>, I)],
    params: BranchParams<I, CoeffType>,
) -> HashMap<I, Complex<CoeffType>>
where
    I: TableauIndex + Send + Sync,
    CoeffType: One + Zero + Clone + Send + Sync + num::Num,
    Complex<CoeffType>:
        std::ops::Mul<Output = Complex<CoeffType>> + std::ops::AddAssign + From<Complex64> + Copy,
{
    let BranchParams {
        stab_anticomm_bits,
        destab_anticomm_bits,
        odd_phase_mask,
        phase_decomp,
        coefficient_factor,
        branch_factor,
    } = params;
    if items.len() >= RAYON_COEFF_THRESHOLD {
        use rayon::prelude::*;

        // Parallel phase: compute all (branch_idx, branch_coeff, idx, nonbranch_coeff) tuples.
        // This is pure math with no shared mutable state.
        let pairs: Vec<(I, Complex<CoeffType>, I, Complex<CoeffType>)> = items
            .par_iter()
            .map(|&(coeff, idx)| {
                let branch_index = idx ^ stab_anticomm_bits;
                let branch_phase_contribution = compute_phase_with_mask_static(
                    destab_anticomm_bits,
                    idx,
                    stab_anticomm_bits,
                    odd_phase_mask,
                );
                let branch_phase = (branch_phase_contribution + phase_decomp) % 4;
                let phase_factor: Complex<CoeffType> =
                    COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();
                (
                    branch_index,
                    phase_factor * coeff * branch_factor,
                    idx,
                    coeff * coefficient_factor,
                )
            })
            .collect();

        // Sequential phase: accumulate into a pre-sized HashMap.
        // HashMap inserts dominate the cost but benefit from cache locality.
        let mut map: HashMap<I, Complex<CoeffType>> =
            HashMap::with_capacity_and_hasher(2 * pairs.len(), Default::default());
        for (branch_idx, branch_coeff, idx, nonbranch_coeff) in pairs {
            *map.entry(branch_idx).or_insert(Complex::zero()) += branch_coeff;
            *map.entry(idx).or_insert(Complex::zero()) += nonbranch_coeff;
        }
        return map;
    }

    branch_coefficients_seq(items.iter().copied(), 2 * items.len(), params)
}

/// Sequential relabel of apply coefficients.
///
/// Pauli application sends every branch to `branch_index = idx ^ stab_anticomm_bits`.
/// XOR by a fixed constant is a bijection, so distinct input indices always map
/// to distinct branch indices: unlike the T-gate branch split (which emits two
/// streams that genuinely collide), the apply path produces no index collisions
/// at all. A per-index coalesce can therefore never merge two entries — the
/// `entry()`-keyed map was pure overhead (hash every key + table allocation) for
/// what is a straight relabel. We instead build a flat `Vec` in one sequential,
/// prefetch-friendly pass and let the caller apply the magnitude cutoff. (The
/// returned keys are unique by the bijection above; the `Vec` backing relies on
/// that, exactly as the old map did implicitly.)
fn apply_coefficients_seq<I, CoeffType>(
    items: impl IntoIterator<Item = (Complex<CoeffType>, I)>,
    capacity: usize,
    stab_anticomm_bits: I,
    destab_anticomm_bits: I,
    odd_phase_mask: I,
    phase_decomp: u8,
) -> Vec<(I, Complex<CoeffType>)>
where
    I: TableauIndex,
    CoeffType: One + Zero + Clone + num::Num,
    Complex<CoeffType>:
        std::ops::Mul<Output = Complex<CoeffType>> + std::ops::AddAssign + From<Complex64> + Copy,
{
    let mut out: Vec<(I, Complex<CoeffType>)> = Vec::with_capacity(capacity);
    for (coeff, idx) in items {
        debug_assert!(
            !(coeff.re == CoeffType::zero() && coeff.im == CoeffType::zero()),
            "Coefficient should not be zero"
        );
        let branch_index = idx ^ stab_anticomm_bits;
        let branch_phase_contribution = compute_phase_with_mask_static(
            destab_anticomm_bits,
            idx,
            stab_anticomm_bits,
            odd_phase_mask,
        );
        let branch_phase = (branch_phase_contribution + phase_decomp) % 4;
        let phase_factor: Complex<CoeffType> =
            COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();
        out.push((branch_index, phase_factor * coeff));
    }
    out
}

/// Relabel coefficients for pauli application. When the coefficient count
/// exceeds `RAYON_COEFF_THRESHOLD`, the per-element relabel runs as a parallel
/// map; below the threshold it falls back to the sequential relabel.
///
/// Because the relabel `idx ^ stab_anticomm_bits` is a bijection (see
/// [`apply_coefficients_seq`]), the parallel map's output already has unique
/// keys — there is nothing to coalesce, so the result is collected straight into
/// a flat `Vec` with no sequential fold afterwards.
#[cfg(feature = "rayon")]
fn apply_coefficients_parallel<I, CoeffType>(
    items: &[(Complex<CoeffType>, I)],
    stab_anticomm_bits: I,
    destab_anticomm_bits: I,
    odd_phase_mask: I,
    phase_decomp: u8,
) -> Vec<(I, Complex<CoeffType>)>
where
    I: TableauIndex + Send + Sync,
    CoeffType: One + Zero + Clone + Send + Sync + num::Num,
    Complex<CoeffType>:
        std::ops::Mul<Output = Complex<CoeffType>> + std::ops::AddAssign + From<Complex64> + Copy,
{
    // See `branch_coefficients_parallel`: avoid nesting rayon inside shot-level
    // parallelism; the main-thread (single-shot) path is unaffected.
    if items.len() >= RAYON_COEFF_THRESHOLD && rayon::current_thread_index().is_none() {
        use rayon::prelude::*;

        return items
            .par_iter()
            .map(|&(coeff, idx)| {
                let branch_index = idx ^ stab_anticomm_bits;
                let branch_phase_contribution = compute_phase_with_mask_static(
                    destab_anticomm_bits,
                    idx,
                    stab_anticomm_bits,
                    odd_phase_mask,
                );
                let branch_phase = (branch_phase_contribution + phase_decomp) % 4;
                let phase_factor: Complex<CoeffType> =
                    COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();
                (branch_index, phase_factor * coeff)
            })
            .collect();
    }

    apply_coefficients_seq(
        items.iter().copied(),
        items.len(),
        stab_anticomm_bits,
        destab_anticomm_bits,
        odd_phase_mask,
        phase_decomp,
    )
}

/// A [`Tableau`] extended with sparse coefficient tracking to handle
/// non-Clifford gates.
///
/// Non-Clifford gates (T, rotations) split a single tableau into a sum
/// of weighted branches indexed by bitstrings. `GeneralizedTableau`
/// stores those weights in a [`SparseVector`] keyed by an
/// [`IndexType`](TableauIndex). Choose:
/// * `IndexType = usize` for up to 64 qubits,
/// * `IndexType = u128` for up to 128,
/// * `IndexType = bnum::types::U256` and friends for the very wide
///   regime.
///
/// Per-qubit loss is tracked in [`is_lost`](GeneralizedTableau::is_lost);
/// gates respect it automatically.
///
/// # Examples
///
/// Prepare a Bell pair and sample one shot. With a fixed seed the two
/// measurements are perfectly correlated on every shot:
///
/// ```
/// use ppvm_pauli_sum::config::fxhash::ByteF64;
/// use ppvm_traits::traits::{Clifford, LossyMeasure};
/// use ppvm_tableau::data::GeneralizedTableau;
///
/// let mut tab: GeneralizedTableau<ByteF64<1>> =
///     GeneralizedTableau::new_with_seed(2, 1e-12, 0);
/// tab.h(0);
/// tab.cnot(0, 1);
///
/// let r0 = LossyMeasure::measure(&mut tab, 0);
/// let r1 = LossyMeasure::measure(&mut tab, 1);
/// assert_eq!(r0, r1);
/// ```
///
/// Non-Clifford gates work through the same interface — apply a `T` gate
/// followed by `T†` and the state is unchanged:
///
/// ```
/// use ppvm_pauli_sum::config::fxhash::ByteF64;
/// use ppvm_traits::traits::{Clifford, TGate};
/// use ppvm_tableau::data::GeneralizedTableau;
///
/// let mut tab: GeneralizedTableau<ByteF64<1>> =
///     GeneralizedTableau::new_with_seed(1, 1e-12, 0);
/// tab.h(0);
/// tab.t(0);
/// tab.t_dag(0);
/// // T followed by T† is the identity; the |+⟩ state is restored.
/// ```
#[derive(Clone)]
pub struct GeneralizedTableau<
    T: Config,
    IndexType = usize,
    SparseVectorType: SparseVector<Complex<T::Coeff>, IndexType> = Vec<(Complex64, IndexType)>,
> {
    /// Underlying Clifford tableau.
    pub tableau: Tableau<T>,
    /// Sparse coefficient vector indexed by bitstrings.
    pub coefficients: SparseVectorType,
    /// Per-qubit loss flags.
    pub is_lost: Vec<bool>,
    /// Coefficient-magnitude threshold below which branches are dropped.
    pub coefficient_threshold: T::Coeff,
    /// Ordered log of every measurement performed (mirrors stim's record).
    pub measurement_record: Vec<Option<bool>>,
    _index_phantom: PhantomData<IndexType>,
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
where
    T::Coeff: One + Zero + Clone + num::Num,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex,
{
    /// Construct a generalized tableau in the `|0…0⟩` state.
    ///
    /// Branches whose coefficient magnitude falls below
    /// `coefficient_threshold` are dropped during gate application.
    pub fn new(n_qubits: usize, coefficient_threshold: T::Coeff) -> Self {
        let mut coefficients = C::new();
        let complex_one = Complex {
            re: T::Coeff::one(),
            im: T::Coeff::zero(),
        };
        coefficients.unsafe_insert(I::zero(), complex_one);
        Self {
            tableau: Tableau::new(n_qubits),
            coefficients,
            is_lost: vec![false; n_qubits],
            coefficient_threshold,
            measurement_record: Vec::new(),
            _index_phantom: PhantomData,
        }
    }

    /// Same as [`GeneralizedTableau::new`], but seed the RNG deterministically.
    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: T::Coeff, seed: u64) -> Self {
        let mut s = Self::new(n_qubits, coefficient_threshold);
        s.tableau.rng = SmallRng::seed_from_u64(seed);
        s
    }

    pub fn reset_all(&mut self) {
        self.tableau.reset_all();

        let mut coefficients = C::new();
        let complex_one = Complex {
            re: T::Coeff::one(),
            im: T::Coeff::zero(),
        };
        coefficients.unsafe_insert(I::zero(), complex_one);
        self.coefficients = coefficients;
        for l in self.is_lost.iter_mut() {
            *l &= false;
        }
        self.measurement_record.clear();
    }

    /// Clone the quantum state but reinitialize the RNG, producing an independent simulation
    /// branch. If `seed` is `Some`, the new RNG is seeded deterministically; if `None`, it is
    /// seeded from OS entropy.
    pub fn fork(&self, seed: Option<u64>) -> Self {
        let mut cloned = self.clone();
        cloned.tableau.rng = match seed {
            Some(s) => SmallRng::seed_from_u64(s),
            None => rand::make_rng(),
        };
        cloned
    }

    /// Number of qubits.
    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits
    }

    /// All measurement outcomes recorded so far, in order.
    pub fn current_measurement_record(&self) -> &[Option<bool>] {
        &self.measurement_record
    }

    /// Append an externally defined measurement result to the record.
    ///
    /// Used by Stim instructions such as `MPAD`, which append measurement
    /// record bits without measuring a qubit.
    pub fn append_measurement_record(&mut self, result: Option<bool>) {
        self.measurement_record.push(result);
    }

    /// Replace the most recent measurement record entry.
    ///
    /// Used by noisy measurement paths where the quantum state follows the
    /// true outcome but the public record should hold the reported bit.
    pub fn overwrite_last_measurement_record(&mut self, result: Option<bool>) {
        if let Some(last) = self.measurement_record.last_mut() {
            *last = result;
        }
    }

    /// Apply CZ to N pairs with constant offset: (base+i, base+offset+i) for i in 0..count.
    /// Falls back to individual CZ calls if any qubit in the range is lost.
    pub fn cz_block_pairs(&mut self, base: usize, offset: usize, count: usize)
    where
        <<T::Storage as BitView>::Store as TryFrom<usize>>::Error: Debug,
        <T::Storage as BitView>::Store: PrimInt + TryFrom<usize>,
    {
        // Check if any qubit in the range is lost
        let any_lost =
            (0..count).any(|i| self.is_lost[base + i] || self.is_lost[base + offset + i]);
        if !any_lost {
            self.tableau.cz_block_pairs(base, offset, count);
        } else {
            // Fallback to individual CZ calls
            for i in 0..count {
                let c = base + i;
                let t = base + offset + i;
                if !self.is_lost[c] && !self.is_lost[t] {
                    Clifford::cz(&mut self.tableau, c, t);
                }
            }
        }
    }

    /// Apply CZ to N cross-word pairs. Controls at word_c, targets at word_t.
    /// Falls back to individual CZ calls if any qubit is lost.
    pub fn cz_block_pairs_cross_word(
        &mut self,
        word_c: usize,
        base_bit_c: usize,
        word_t: usize,
        base_bit_t: usize,
        count: usize,
    ) where
        <T::Storage as BitView>::Store: PrimInt,
    {
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let any_lost = (0..count).any(|i| {
            let c = word_c * bits_per_word + base_bit_c + i;
            let t = word_t * bits_per_word + base_bit_t + i;
            self.is_lost[c] || self.is_lost[t]
        });
        if !any_lost {
            self.tableau
                .cz_block_pairs_cross_word(word_c, base_bit_c, word_t, base_bit_t, count);
        } else {
            for i in 0..count {
                let c = word_c * bits_per_word + base_bit_c + i;
                let t = word_t * bits_per_word + base_bit_t + i;
                if !self.is_lost[c] && !self.is_lost[t] {
                    Clifford::cz(&mut self.tableau, c, t);
                }
            }
        }
    }

    /// Apply CZ to `count` pairs with a constant offset, given in qubit-index
    /// terms: `(control_base + i, target_base + i)` for `i in 0..count`.
    ///
    /// This is the high-level entry point for a fused block of CZs: it splits
    /// the run at storage-word boundaries internally and dispatches each
    /// segment to [`Self::cz_block_pairs`] (control and target in the same
    /// word) or [`Self::cz_block_pairs_cross_word`] (straddling two words), so
    /// callers never need to reason about the `u64` packing. CZ is symmetric,
    /// so the two bases may be passed in either order.
    pub fn cz_block(&mut self, control_base: usize, target_base: usize, count: usize)
    where
        <<T::Storage as BitView>::Store as TryFrom<usize>>::Error: Debug,
        <T::Storage as BitView>::Store: PrimInt + TryFrom<usize>,
    {
        if count == 0 {
            return;
        }
        // cz_block_pairs needs a non-negative offset; CZ is symmetric, so order
        // the two bases.
        let (lo, hi) = if control_base <= target_base {
            (control_base, target_base)
        } else {
            (target_base, control_base)
        };
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let mut i = 0;
        while i < count {
            let (c, t) = (lo + i, hi + i);
            let (wc, bc) = (c / bits_per_word, c % bits_per_word);
            let (wt, bt) = (t / bits_per_word, t % bits_per_word);
            // Longest run before either index crosses into the next word.
            let run = (bits_per_word - bc).min(bits_per_word - bt).min(count - i);
            if wc == wt {
                self.cz_block_pairs(c, t - c, run);
            } else {
                self.cz_block_pairs_cross_word(wc, bc, wt, bt, run);
            }
            i += run;
        }
    }

    // helper functions

    /// Compute the decomposition of a pauli into stabilizer destabilizer products
    /// Any Pauli can be written as P_addr0 = phase * prod(d_k ^ gamma_k) * prod(s_l ^ lambda_l)
    /// where: gamma_k == 1 iff {P_addr0, s_k} = 0
    /// lambda_l == 1 iff {P_addr0, d_l} = 0
    /// Lemma 5. from T. J. Yoder (2012)
    /// NOTE: this is O(n^2)
    ///
    /// The function returns `(phase, stab_anticomm_bits, destab_anticomm_bits)`, where
    /// `stab_anticomm_bits[k] = 1` iff P_addr0 anticommutes with stabilizer s_k, and
    /// `destab_anticomm_bits[l] = 1` iff P_addr0 anticommutes with destabilizer d_l.
    /// Note that stab_anticomm_bits is equal to the shift of the index when branching
    /// (`beta` in Eq(4) of the SOFT paper).
    pub fn compute_decomposition(&self, addr0: usize, pauli: Pauli) -> (u8, I, I)
    where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        let n = self.n_qubits();

        // the actual decomposition, which we need to track the phase
        let mut p_word = PhasedPauliWordNoHash::<T::Storage, T::BuildHasher>::new(n);
        p_word.set(addr0, pauli);

        // the bit strings defining the contributions
        let mut destab_anticomm_bits = I::zero();
        let mut stab_anticomm_bits = I::zero();

        debug_assert_ne!(pauli, Pauli::I);
        let pauli_bits = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
            _ => unreachable!("Pauli L cannot occur in tableau"),
        };

        let stabilizers = self.tableau.stabilizers();
        let destabilizers = self.tableau.destabilizers();
        let one = I::one();

        for (i, stab) in stabilizers.iter().enumerate() {
            if !destabilizers[i].word.anticommutes_at(addr0, pauli_bits) {
                // commutes
                continue;
            }

            // contributes, so set the corresponding bit in destab_anticomm_bits to 1
            destab_anticomm_bits |= one << i;

            // destabilizer anti-commutes, so the stabilizer contributes
            // the stabilizer is its own inverse up to its phase
            // to avoid inverting the stabilizer, we just multiply by it
            // and then divide by its phase squared
            p_word *= stab;
            p_word.add_phase(8 - 2 * stab.phase);
        }

        // NOTE: destabilizers also commute with one another in a valid tableau
        // since the form a basis together with stabilizers
        for (i, destab) in destabilizers.iter().enumerate() {
            if !stabilizers[i].word.anticommutes_at(addr0, pauli_bits) {
                // commutes
                continue;
            }

            // contributes, so set the corresponding bit in stab_anticomm_bits to 1
            stab_anticomm_bits |= one << i;

            // stabilizer anti-commutes, so the destabilizer contributes
            p_word *= destab;
            p_word.add_phase(8 - 2 * destab.phase);
        }

        (p_word.phase, stab_anticomm_bits, destab_anticomm_bits)
    }

    /// Multi-qubit generalization of [`compute_decomposition`]: conjugate an
    /// arbitrary `PauliWord` through the tableau and return the same triple
    /// `(phase, stab_anticomm_bits, destab_anticomm_bits)`.
    ///
    /// Algorithm: call [`compute_decomposition`] for each non-identity qubit
    /// in the input, then multiply the resulting single-qubit conjugates in
    /// canonical-basis form `i^φ X^x Z^z`. Pauli multiplication picks up a
    /// `(-1)^{popcount(z_running & x_new)}` cross-phase from
    /// `Z^z_a X^x_b = (-1)^{z_a · x_b} X^x_b Z^z_a`.
    pub(crate) fn compute_decomposition_word<W: PauliWordTrait>(&self, word: &W) -> (u8, I, I)
    where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        let mut phase = 0u8;
        let mut stab_anticomm = I::zero();
        let mut destab_anticomm = I::zero();
        for q in 0..self.n_qubits() {
            let p_q = word.get(q);
            if p_q == Pauli::I {
                continue;
            }
            let (q_phase, q_stab, q_destab) = self.compute_decomposition(q, p_q);
            let cross = 2 * (symplectic_inner(destab_anticomm, q_stab) as u8 % 2);
            phase = (phase + q_phase + cross) % 4;
            stab_anticomm = stab_anticomm ^ q_stab;
            destab_anticomm = destab_anticomm ^ q_destab;
        }
        (phase, stab_anticomm, destab_anticomm)
    }

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    #[cfg(test)]
    pub(crate) fn compute_phase(
        &self,
        destab_anticomm_bits: I,
        basis_index: I,
        stab_anticomm_bits: I,
    ) -> u8 {
        // phase convention: 0: +1, 1: +i, 2: -1, 3: -i
        let one = I::one();
        let zero = I::zero();

        // contribution 1: each destabilizer D_i with basis_index[i]=1 that anticommutes
        // with P (destab_anticomm_bits[i]=1) contributes a -1 sign; this is the symplectic inner product
        let mut phase = (2 * symplectic_inner(destab_anticomm_bits, basis_index) as u8) % 4;

        // contribution 2: destabilizers that appear twice (basis_index[i]=1 and stab_anticomm_bits[i]=1)
        // contribute an extra -1 if their phase is imaginary
        let active = basis_index & stab_anticomm_bits;
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if active & (one << i) == zero {
                continue;
            }
            if destab.phase % 2 != 0 {
                phase = (phase + 2) % 4;
            }
        }

        phase
    }

    /// Build a bitmask where bit i is set if destabilizer i has odd phase (phase % 2 != 0).
    pub fn odd_phase_destabilizer_mask(&self) -> I {
        let mut mask = I::zero();
        let one = I::one();
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if destab.phase % 2 != 0 {
                mask |= one << i;
            }
        }
        mask
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
where
    T::Coeff: One + Zero + Clone + Send + Sync + num::Num + PartialOrd,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync,
{
    pub(crate) fn branch_with_coefficients(
        &mut self,
        addr0: usize,
        pauli: Pauli,
        coefficient_factor: Complex<T::Coeff>,
        branch_factor: Complex<T::Coeff>,
    ) where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, pauli);

        let odd_phase_mask = self.odd_phase_destabilizer_mask();
        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let n_coefficients = old_coefficients.len();

        let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();

        // When rayon is enabled and above the threshold, use parallel map then apply cutoff
        // directly. Early-return so old_coefficients stays valid for the sequential path below.
        #[cfg(feature = "rayon")]
        if n_coefficients >= RAYON_COEFF_THRESHOLD {
            let items: Vec<_> = old_coefficients.into_iter().collect();
            let map = branch_coefficients_parallel(
                &items,
                BranchParams {
                    stab_anticomm_bits,
                    destab_anticomm_bits,
                    odd_phase_mask,
                    phase_decomp,
                    coefficient_factor,
                    branch_factor,
                },
            );
            self.coefficients.reserve(map.len());
            for (idx, coeff) in map {
                if coeff.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(idx, coeff);
                }
            }
            return;
        }

        // Sequential sort-merge path: build nonbranch (nb), branch-values (brv), and a
        // packed branch-key stream, sort if needed, then 2-way merge directly into
        // self.coefficients with inline cutoff — no intermediate output Vec.
        //
        // Fast path: pack each (branch_key, build_pos) as a single u64
        // `(key << 16) | pos` and sort the u64 array. Half the data movement vs the
        // generic (I, u32) elements, and the most-optimised sort_unstable path.
        // Preconditions for packable: n_coefficients ≤ 0xFFFF (pos fits 16 bits) and
        // every branch key fits in 47 bits. Both hold for cultivation_d5 (≤13 active
        // bits, small coefficient counts). The fallback reproduces the prior behaviour
        // exactly for wide index types or large keys.
        let mut nb: Vec<(I, Complex<T::Coeff>)> = Vec::with_capacity(n_coefficients);
        let mut brv: Vec<Complex<T::Coeff>> = Vec::with_capacity(n_coefficients);
        let mut packed: Vec<u64> = Vec::with_capacity(n_coefficients);
        let mut packable = n_coefficients <= 0xFFFF;
        let mut nb_sorted = true;
        let mut prev: Option<I> = None;
        for (pos, (coeff, idx)) in (0_u32..).zip(old_coefficients) {
            let branch_index = idx ^ stab_anticomm_bits;
            let bpc = compute_phase_with_mask_static(
                destab_anticomm_bits,
                idx,
                stab_anticomm_bits,
                odd_phase_mask,
            );
            let branch_phase = (bpc + phase_decomp) % 4;
            let pf: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();
            brv.push(pf * coeff * branch_factor);
            match <u64 as num::NumCast>::from(branch_index) {
                Some(k) if k < (1u64 << 47) => packed.push((k << 16) | (pos as u64)),
                _ => {
                    packable = false;
                    packed.push(pos as u64);
                }
            }
            nb.push((idx, coeff * coefficient_factor));
            if let Some(p) = prev
                && idx < p
            {
                nb_sorted = false;
            }
            prev = Some(idx);
        }

        self.coefficients.reserve(nb.len() + brv.len());
        let mut i = 0;
        if packable {
            if !nb_sorted {
                nb.sort_unstable_by_key(|a| a.0);
            }
            packed.sort_unstable();
            // Decode the 47-bit key from the high bits of a packed entry using byte-by-byte
            // construction from I::from(u8) — avoids num::NumCast which panics for bnum types.
            let decode_key = |w: u64| -> I {
                let k = w >> 16; // k < 2^47; 6 bytes suffice
                let mut v = I::zero();
                for b in 0..6usize {
                    let byte = ((k >> (b * 8)) & 0xFF) as u8;
                    v |= <I as From<u8>>::from(byte) << (b * 8);
                }
                v
            };
            let mut j = 0;
            while i < nb.len() && j < packed.len() {
                let bp = (packed[j] & 0xFFFF) as usize;
                let bk = decode_key(packed[j]);
                match nb[i].0.cmp(&bk) {
                    std::cmp::Ordering::Less => {
                        if nb[i].1.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
                        }
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        let v = brv[bp];
                        if v.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(bk, v);
                        }
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        let mut sv = nb[i].1;
                        sv += brv[bp];
                        if sv.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, sv);
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
            while j < packed.len() {
                let bp = (packed[j] & 0xFFFF) as usize;
                let bk = decode_key(packed[j]);
                let v = brv[bp];
                if v.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(bk, v);
                }
                j += 1;
            }
        } else {
            // Fallback for wide index types, large keys (≥ 2^47), or many coefficients
            // (> 65535). nb is still in build order here; reconstruct brk from it before
            // sorting nb, so build-position p correctly indexes brv[p].
            let mut brk: Vec<(I, u32)> = (0_u32..)
                .zip(nb.iter())
                .map(|(p, &(idx, _))| (idx ^ stab_anticomm_bits, p))
                .collect();
            if !nb_sorted {
                nb.sort_unstable_by_key(|a| a.0);
            }
            brk.sort_unstable_by_key(|a| a.0);
            let mut j = 0;
            while i < nb.len() && j < brk.len() {
                let (bk, bp) = brk[j];
                match nb[i].0.cmp(&bk) {
                    std::cmp::Ordering::Less => {
                        if nb[i].1.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
                        }
                        i += 1;
                    }
                    std::cmp::Ordering::Greater => {
                        let v = brv[bp as usize];
                        if v.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(bk, v);
                        }
                        j += 1;
                    }
                    std::cmp::Ordering::Equal => {
                        let mut sv = nb[i].1;
                        sv += brv[bp as usize];
                        if sv.norm_sqr() > cutoff_sq {
                            self.coefficients.unsafe_insert(nb[i].0, sv);
                        }
                        i += 1;
                        j += 1;
                    }
                }
            }
            while j < brk.len() {
                let (bk, bp) = brk[j];
                let v = brv[bp as usize];
                if v.norm_sqr() > cutoff_sq {
                    self.coefficients.unsafe_insert(bk, v);
                }
                j += 1;
            }
        }
        while i < nb.len() {
            if nb[i].1.norm_sqr() > cutoff_sq {
                self.coefficients.unsafe_insert(nb[i].0, nb[i].1);
            }
            i += 1;
        }
    }

    pub(crate) fn compute_coefficients_after_pauli_apply(
        &self,
        coefficients: &mut C,
        addr0: usize,
        pauli: Pauli,
    ) where
        <<T as Config>::Storage as BitView>::Store: PrimInt,
    {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, pauli);

        let odd_phase_mask = self.odd_phase_destabilizer_mask();
        let n_coefficients = coefficients.len();
        let old_coefficients = std::mem::replace(coefficients, C::new());

        #[cfg(feature = "rayon")]
        let new_coefficients = if n_coefficients >= RAYON_COEFF_THRESHOLD {
            let items: Vec<_> = old_coefficients.into_iter().collect();
            apply_coefficients_parallel(
                &items,
                stab_anticomm_bits,
                destab_anticomm_bits,
                odd_phase_mask,
                phase_decomp,
            )
        } else {
            apply_coefficients_seq(
                old_coefficients,
                n_coefficients,
                stab_anticomm_bits,
                destab_anticomm_bits,
                odd_phase_mask,
                phase_decomp,
            )
        };

        #[cfg(not(feature = "rayon"))]
        let new_coefficients = apply_coefficients_seq(
            old_coefficients,
            n_coefficients,
            stab_anticomm_bits,
            destab_anticomm_bits,
            odd_phase_mask,
            phase_decomp,
        );

        let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();
        coefficients.reserve(new_coefficients.len());
        for (idx, coeff) in new_coefficients {
            if coeff.norm_sqr() > cutoff_sq {
                coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::BUint;
    use ppvm_pauli_sum::config::fxhash::ByteF64;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;
    type TestTableauBUint = GeneralizedTableau<TestConfig, BUint<1>>;

    #[test]
    fn test_compute_phase_z_2_single_qubit_plus_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);

        // After H: stabilizer = +X, destabilizer = +Z
        // stab_anticomm_bits = 1 (stabilizer has xbit[0]=true)
        // both phases should be 0
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            tab.compute_decomposition(0, Pauli::Z);
        let phase0 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 0, stab_anticomm_bits);
        assert_eq!(phase0, 0);
        let phase1 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 1, stab_anticomm_bits);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            tab.compute_decomposition(0, Pauli::Z);
        let phase0 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 0, stab_anticomm_bits);
        assert_eq!(phase0, 0);
        let phase1 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 1, stab_anticomm_bits);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_mx_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.z(0);

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            tab.compute_decomposition(0, Pauli::Z);
        let phase0 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 0, stab_anticomm_bits);
        assert_eq!(phase0, 0);
        let phase1 = phase_decomp + tab.compute_phase(destab_anticomm_bits, 1, stab_anticomm_bits);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer_2() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);
        tab.tableau.h(0);

        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            tab.compute_decomposition(0, Pauli::Z);
        let phase0 =
            (phase_decomp + tab.compute_phase(destab_anticomm_bits, 0, stab_anticomm_bits)) % 4;
        assert_eq!(phase0, 1);
        let phase1 =
            (phase_decomp + tab.compute_phase(destab_anticomm_bits, 1, stab_anticomm_bits)) % 4;
        assert_eq!(phase1, 3);
    }

    #[test]
    fn test_index_type() {
        let mut tab: TestTableauBUint = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
    }

    /// Snapshot all rows (xbits, zbits, phase) for comparison.
    fn snapshot_tableau<C: Config>(tab: &Tableau<C>) -> Vec<(C::Storage, C::Storage, u8)> {
        tab.data
            .iter()
            .map(|pw| (pw.word.xbits.data, pw.word.zbits.data, pw.phase))
            .collect()
    }

    #[test]
    fn test_cz_block_pairs_matches_individual() {
        // Test CZ pairs (0,4), (1,5), (2,6), (3,7) — offset=4, count=4
        type TTab = Tableau<ByteF64<1>>;
        let n = 8;
        let base = 0;
        let offset = 4;
        let count = 4;

        let mut tab1 = TTab::new(n);
        Clifford::h(&mut tab1, 0);
        Clifford::h(&mut tab1, 3);
        Clifford::s(&mut tab1, 1);
        let mut tab2 = tab1.clone();

        // Individual
        for i in 0..count {
            Clifford::cz(&mut tab1, base + i, base + offset + i);
        }

        // Batch
        tab2.cz_block_pairs(base, offset, count);

        assert_eq!(snapshot_tableau(&tab1), snapshot_tableau(&tab2));
    }

    #[test]
    fn test_cz_block_pairs_offset_17() {
        // Simulate MSD-like CZ: (0,17), (1,18), ..., (16,33) — all in one u64 word
        use ppvm_pauli_sum::config::fx64hash::Byte8F64;
        type LargeTab = Tableau<Byte8F64<2>>;
        let n = 34;
        let mut tab1 = LargeTab::new(n);
        // Apply some gates to create non-trivial state
        for i in 0..n {
            Clifford::h(&mut tab1, i);
        }
        let mut tab2 = tab1.clone();

        // Individual
        for i in 0..17 {
            Clifford::cz(&mut tab1, i, 17 + i);
        }

        // Batch
        tab2.cz_block_pairs(0, 17, 17);

        assert_eq!(snapshot_tableau(&tab1), snapshot_tableau(&tab2));
    }

    #[test]
    fn test_cz_block_pairs_nonzero_base() {
        // Test CZ pairs starting from a non-zero base: (10,27), (11,28), ..., (14,31)
        // All within one u64 word (bits 0-63)
        use ppvm_pauli_sum::config::fx64hash::Byte8F64;
        type LargeTab = Tableau<Byte8F64<2>>;
        let n = 32;
        let base = 10;
        let offset = 17;
        let count = 5;

        let mut tab1 = LargeTab::new(n);
        for i in 0..n {
            Clifford::h(&mut tab1, i);
        }
        Clifford::s(&mut tab1, 12);
        Clifford::s(&mut tab1, 28);
        let mut tab2 = tab1.clone();

        for i in 0..count {
            Clifford::cz(&mut tab1, base + i, base + offset + i);
        }

        tab2.cz_block_pairs(base, offset, count);

        assert_eq!(snapshot_tableau(&tab1), snapshot_tableau(&tab2));
    }

    #[test]
    fn test_cz_block_pairs_single_pair() {
        // Degenerate case: count=1 should be same as one CZ
        type TTab = Tableau<ByteF64<1>>;
        let n = 8;
        let mut tab1 = TTab::new(n);
        Clifford::h(&mut tab1, 2);
        Clifford::s(&mut tab1, 5);
        let mut tab2 = tab1.clone();

        Clifford::cz(&mut tab1, 2, 5);
        tab2.cz_block_pairs(2, 3, 1);

        assert_eq!(snapshot_tableau(&tab1), snapshot_tableau(&tab2));
    }

    #[test]
    fn test_cz_block_pairs_zero_count() {
        // count=0 should be a no-op
        type TTab = Tableau<ByteF64<1>>;
        let n = 8;
        let mut tab1 = TTab::new(n);
        Clifford::h(&mut tab1, 0);
        let before = snapshot_tableau(&tab1);
        tab1.cz_block_pairs(0, 4, 0);
        assert_eq!(before, snapshot_tableau(&tab1));
    }

    #[test]
    fn test_generalized_tableau_cz_block_pairs() {
        // Test through GeneralizedTableau wrapper
        use ppvm_pauli_sum::config::fx64hash::Byte8F64;
        type GTab = GeneralizedTableau<Byte8F64<2>>;
        let n = 34;
        let mut tab1: GTab = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            Clifford::h(&mut tab1.tableau, i);
        }
        let mut tab2 = tab1.clone();

        // Individual via Clifford trait
        for i in 0..17 {
            Clifford::cz(&mut tab1, i, 17 + i);
        }

        // Batch
        tab2.cz_block_pairs(0, 17, 17);

        assert_eq!(
            snapshot_tableau(&tab1.tableau),
            snapshot_tableau(&tab2.tableau)
        );
    }

    #[test]
    fn test_generalized_tableau_cz_block_pairs_with_loss() {
        // When a qubit is lost, should fall back to individual CZ (skipping lost ones)
        type GTab = GeneralizedTableau<ByteF64<1>>;
        let n = 8;
        let mut tab1: GTab = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            Clifford::h(&mut tab1.tableau, i);
        }
        tab1.is_lost[2] = true; // Mark qubit 2 as lost
        let mut tab2 = tab1.clone();

        // Individual, skipping lost qubits
        for i in 0..4 {
            let c = i;
            let t = 4 + i;
            if !tab1.is_lost[c] && !tab1.is_lost[t] {
                Clifford::cz(&mut tab1.tableau, c, t);
            }
        }

        // Batch (should fall back internally)
        tab2.cz_block_pairs(0, 4, 4);

        assert_eq!(
            snapshot_tableau(&tab1.tableau),
            snapshot_tableau(&tab2.tableau)
        );
    }

    #[test]
    fn test_cz_block_matches_individual_across_word_boundary() {
        // cz_block must split a run that straddles the u64 boundary into the
        // right within-word + cross-word segments. control_base=34,
        // target_base=51, count=17 reproduces the MSD ql[2]xql[3] sweep:
        // (34,51)..(46,63) in word 0, then (47,64)..(50,67) cross-word.
        use ppvm_pauli_sum::config::fx64hash::Byte8F64;
        type GTab = GeneralizedTableau<Byte8F64<2>>;
        let n = 85;
        let mut tab1: GTab = GeneralizedTableau::new(n, 1e-12);
        for i in 0..n {
            Clifford::h(&mut tab1.tableau, i);
        }
        let mut tab2 = tab1.clone();

        let (control_base, target_base, count) = (34, 51, 17);
        for i in 0..count {
            Clifford::cz(&mut tab1, control_base + i, target_base + i);
        }
        tab2.cz_block(control_base, target_base, count);

        assert_eq!(
            snapshot_tableau(&tab1.tableau),
            snapshot_tableau(&tab2.tableau)
        );

        // Reversed bases (CZ is symmetric) must give the same result.
        let mut tab3 = GeneralizedTableau::<Byte8F64<2>>::new(n, 1e-12);
        for i in 0..n {
            Clifford::h(&mut tab3.tableau, i);
        }
        tab3.cz_block(target_base, control_base, count);
        assert_eq!(
            snapshot_tableau(&tab1.tableau),
            snapshot_tableau(&tab3.tableau)
        );
    // ─── reset_all ────────────────────────────────────────────────────

    /// `GeneralizedTableau::reset_all` restores the full state to a fresh
    /// `|0…0⟩` tableau: identical stabilizer/destabilizer rows and a single
    /// identity coefficient, even after non-Clifford branching.
    #[test]
    fn reset_all_restores_fresh_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(3, 1e-12);
        let fresh: TestTableau = GeneralizedTableau::new(3, 1e-12);

        tab.h(0);
        tab.cnot(0, 1);
        tab.ry(2, 0.7); // non-Clifford: branches the coefficient vector
        assert!(
            tab.coefficients.iter().count() > 1,
            "rotation should branch the coefficient vector"
        );

        tab.reset_all();

        assert_eq!(
            snapshot_tableau(&tab.tableau),
            snapshot_tableau(&fresh.tableau)
        );
        let coeffs: Vec<_> = tab.coefficients.iter().copied().collect();
        let fresh_coeffs: Vec<_> = fresh.coefficients.iter().copied().collect();
        assert_eq!(coeffs, fresh_coeffs);
    }

    /// A full reset clears the measurement record. Regression guard: an earlier
    /// version left it intact, so `current_measurement_record` returned stale
    /// outcomes after a reset.
    #[test]
    fn reset_all_clears_measurement_record() {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.append_measurement_record(Some(true));
        tab.append_measurement_record(None);
        assert_eq!(tab.current_measurement_record().len(), 2);

        tab.reset_all();

        assert!(tab.current_measurement_record().is_empty());
    }

    /// A full reset clears per-qubit loss flags.
    #[test]
    fn reset_all_clears_loss_flags() {
        let mut tab: TestTableau = GeneralizedTableau::new(3, 1e-12);
        tab.is_lost[0] = true;
        tab.is_lost[2] = true;

        tab.reset_all();

        assert!(tab.is_lost.iter().all(|&lost| !lost));
    }

    /// `Tableau::reset_all` restores the fresh identity tableau rows.
    #[test]
    fn tableau_reset_all_restores_fresh_rows() {
        let mut tab: Tableau<TestConfig> = Tableau::new(4);
        let fresh: Tableau<TestConfig> = Tableau::new(4);

        tab.h(0);
        tab.s(1);
        tab.h(3);

        tab.reset_all();

        assert_eq!(snapshot_tableau(&tab), snapshot_tableau(&fresh));
    }
}
