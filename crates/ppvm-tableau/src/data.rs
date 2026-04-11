use std::{fmt::Debug, marker::PhantomData};

use fxhash::FxHashMap as HashMap;

use bitvec::array::BitArray;
use bitvec::view::BitView;
use num::PrimInt;

use crate::prelude::*;
use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

type PhasedPauliWordNoHash<A, H> = PhasedPauliWord<A, H, PauliWord<A, H, false>>;

#[derive(Clone, Debug)]
pub struct Tableau<T: Config> {
    pub n_qubits: usize,
    /// Destabilizer / Stabilizer tableau
    /// * Entries 0..n are the destabilizers
    /// * Entries n..2n are the stabilizers
    pub data: Vec<PhasedPauliWordNoHash<T::Storage, T::BuildHasher>>,
    pub(crate) rng: SmallRng,
}

impl<T: Config> Tableau<T> {
    pub fn new(n_qubits: usize) -> Self {
        // Initialize tableau for 0 state
        let mut data: Vec<PhasedPauliWordNoHash<T::Storage, T::BuildHasher>> =
            Vec::with_capacity(2 * n_qubits);
        let pw_cache = PhasedPauliWordNoHash::<T::Storage, T::BuildHasher>::new(n_qubits);
        for i in 0..n_qubits {
            // destabilizer
            let mut pw = pw_cache.clone();
            pw.set(i, Pauli::X);
            data.push(pw);
        }
        for i in 0..n_qubits {
            // stabilizer
            let mut pw = pw_cache.clone();
            pw.set(i, Pauli::Z);
            data.push(pw);
        }

        Self {
            n_qubits,
            data,
            rng: rand::make_rng(),
        }
    }

    pub fn new_with_seed(n_qubits: usize, seed: u64) -> Self {
        let mut t = Self::new(n_qubits);
        t.rng = SmallRng::seed_from_u64(seed);
        t
    }

    #[inline]
    pub fn stabilizers(&self) -> &[PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &self.data[self.n_qubits..]
    }

    #[inline]
    pub fn stabilizers_mut(&mut self) -> &mut [PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &mut self.data[self.n_qubits..]
    }

    #[inline]
    pub fn destabilizers(&self) -> &[PhasedPauliWordNoHash<T::Storage, T::BuildHasher>] {
        &self.data[..self.n_qubits]
    }

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

        // Clone g_q once before the loop
        let g_q = stabilizers[q_idx].clone();

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
}

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },  // +1
    Complex64 { re: 0.0, im: 1.0 },  // +i
    Complex64 { re: -1.0, im: 0.0 }, // -1
    Complex64 { re: 0.0, im: -1.0 }, // -i
];

#[inline]
pub fn symplectic_inner<I>(alpha: I, beta: I) -> u32
where
    I: TableauIndex,
{
    (alpha & beta).count_ones()
}

// TODO: builder
#[derive(Clone)]
pub struct GeneralizedTableau<
    T: Config,
    IndexType = usize,
    SparseVectorType: SparseVector<Complex<T::Coeff>, IndexType> = Vec<(Complex64, IndexType)>,
> {
    pub tableau: Tableau<T>,
    pub coefficients: SparseVectorType,
    pub is_lost: Vec<bool>,
    pub coefficient_threshold: T::Coeff,
    _index_phantom: PhantomData<IndexType>,
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
where
    T::Coeff: One + Zero + Clone,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + std::ops::AddAssign
        + From<Complex64>
        + ComplexFloat,
    I: TableauIndex,
{
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
            _index_phantom: PhantomData,
        }
    }

    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: T::Coeff, seed: u64) -> Self {
        let mut s = Self::new(n_qubits, coefficient_threshold);
        s.tableau.rng = SmallRng::seed_from_u64(seed);
        s
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

    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits
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
    pub(crate) fn compute_decomposition(&self, addr0: usize, pauli: Pauli) -> (u8, I, I)
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

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    #[allow(dead_code)]
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
    pub(crate) fn odd_phase_destabilizer_mask(&self) -> I {
        let mut mask = I::zero();
        let one = I::one();
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if destab.phase % 2 != 0 {
                mask |= one << i;
            }
        }
        mask
    }

    /// Like `compute_phase`, but uses a precomputed odd-phase bitmask instead of
    /// looping over all destabilizers. The mask should be obtained from
    /// `odd_phase_destabilizer_mask()`.
    pub(crate) fn compute_phase_with_mask(
        &self,
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

    /// Keep only coefficients that correspond to the correct eigenvalue of a Z measurement.
    /// Applying the projector to a basis state, we have three phases:
    /// 1. The actual measurement outcome (k)
    /// 2. The sign from whether +Z or -Z is a stabilizer (m) - can get that from the decomposition
    /// 3. Contribution from commuting Z_addr0 through the destabilizers (xi)
    ///    Only coefficients where m*k*xi == 1 are kept, equivalently written as (xi * k) == m
    pub(crate) fn trim_coefficients_for_measurement(
        &mut self,
        destab_anticomm_bits: I,
        z_sign: bool,
        outcome: bool,
    ) {
        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let old_len = old_coefficients.len();
        for (coeff, alpha) in old_coefficients.into_iter() {
            let mut phase = false; // false: +1 eigenspace of Z, true: -1 eigenspace

            // get the phase from the anti-commutation with the product over all destabilizers
            let parity = symplectic_inner(alpha, destab_anticomm_bits) % 2 != 0;
            phase ^= parity;

            // (xi * k) == m
            if (phase ^ z_sign) == outcome {
                self.coefficients.unsafe_insert(alpha, coeff);
            }
        }

        // renormalize only if coefficients were actually trimmed
        if self.coefficients.len() < old_len {
            self.coefficients.normalize();
        }
    }

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
        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::default();
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ stab_anticomm_bits; // stab_anticomm_bits is the index shift

            // get the phase contributions from duplicate destabilizers
            // and anti-commuting through destabilizers
            let branch_phase_contribution = self.compute_phase_with_mask(
                destab_anticomm_bits,
                idx,
                stab_anticomm_bits,
                odd_phase_mask,
            );

            // the total phase is the product of the above with the decomposition phase
            let branch_phase = (branch_phase_contribution + phase_decomp) % 4;

            let phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            let branch_coefficient = phase_factor * coeff * branch_factor;
            let nonbranch_coefficient = coeff * coefficient_factor;

            *new_coefficients
                .entry(branch_index)
                .or_insert(Complex::zero()) += branch_coefficient;
            *new_coefficients.entry(idx).or_insert(Complex::zero()) += nonbranch_coefficient;
        }

        let cutoff = Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        }
        .abs();
        for (idx, coeff) in new_coefficients {
            if coeff.abs() > cutoff {
                self.coefficients.unsafe_insert(idx, coeff);
            }
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
        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::default();
        let old_coefficients = std::mem::replace(coefficients, C::new());
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ stab_anticomm_bits; // stab_anticomm_bits is the index shift

            // get the phase contributions from duplicate destabilizers
            // and anti-commuting through destabilizers
            let branch_phase_contribution = self.compute_phase_with_mask(
                destab_anticomm_bits,
                idx,
                stab_anticomm_bits,
                odd_phase_mask,
            );
            let branch_phase = (branch_phase_contribution + phase_decomp) % 4;

            let phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            let branch_coefficient = phase_factor * coeff;

            *new_coefficients
                .entry(branch_index)
                .or_insert(Complex::zero()) += branch_coefficient;
        }

        let cutoff = Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        }
        .abs();

        for (idx, coeff) in new_coefficients {
            if coeff.abs() > cutoff {
                coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bnum::BUint;
    use ppvm_runtime::config::fxhash::ByteF64;

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
}
