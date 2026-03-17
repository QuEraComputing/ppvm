use std::{collections::HashMap, fmt::Debug, marker::PhantomData};

use super::sparsevec::SparseVector;
use crate::traits::PauliWordTrait;
use crate::{char::Pauli, config::Config};
use crate::{phase::PhasedPauliWord, tableau::traits::TableauIndex};
use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};
use rand::SeedableRng;
use rand::rngs::SmallRng;

#[derive(Clone, Debug)]
pub struct Tableau<T: Config> {
    pub n_qubits: usize,
    /// Destabilizer / Stabilizer tableau
    /// * Entries 0..n are the destabilizers
    /// * Entries n..2n are the stabilizers
    pub data: Vec<PhasedPauliWord<T::Storage, T::BuildHasher>>,
    pub(crate) rng: SmallRng,
}

impl<T: Config> Tableau<T> {
    pub fn new(n_qubits: usize) -> Self {
        // Initialize tableau for 0 state
        let mut data: Vec<PhasedPauliWord<T::Storage, T::BuildHasher>> =
            Vec::with_capacity(2 * n_qubits);
        let pw_cache = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(n_qubits);
        for i in 0..n_qubits {
            // destabilizer
            let mut pw = pw_cache.clone();
            pw.set(i, crate::char::Pauli::X);
            data.push(pw);
        }
        for i in 0..n_qubits {
            // stabilizer
            let mut pw = pw_cache.clone();
            pw.set(i, crate::char::Pauli::Z);
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
    pub fn stabilizers(&self) -> &[PhasedPauliWord<T::Storage, T::BuildHasher>] {
        &self.data[self.n_qubits..]
    }

    #[inline]
    pub fn stabilizers_mut(&mut self) -> &mut [PhasedPauliWord<T::Storage, T::BuildHasher>] {
        &mut self.data[self.n_qubits..]
    }

    #[inline]
    pub fn destabilizers(&self) -> &[PhasedPauliWord<T::Storage, T::BuildHasher>] {
        &self.data[..self.n_qubits]
    }

    #[inline]
    pub fn destabilizers_mut(&mut self) -> &mut [PhasedPauliWord<T::Storage, T::BuildHasher>] {
        &mut self.data[..self.n_qubits]
    }

    // some helper functions for measurement impl
    pub(crate) fn find_z_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        // Find first stabilizer that anticommutes with Z_addr0
        self.stabilizers()
            .iter()
            .position(|stab| stab.word.anticommutes_at(addr0, (false, true)))
    }

    pub(crate) fn get_deterministic_outcome(&self, addr0: usize) -> bool {
        // find the outcome: either Z_addr0 or -Z_addr0 is a stabilizer
        // the stabilizer can be computed as the product of all stabilizers S_i
        // whose corresponding destabilizer D_i anticommutes with Z_addr0 (has X at addr0).
        // We have to actually multiply Paulis to also account for products of +i/-i
        let destabilizers = self.destabilizers();
        let stabilizers = self.stabilizers();
        let n = self.n_qubits;
        let mut result = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(n);
        for (i, destab) in destabilizers.iter().enumerate() {
            if destab.word.xbits[addr0] {
                result *= stabilizers[i].clone();
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
    ) {
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
                stabilizers[i] *= g_q.clone();
            }
            if destabilizers[i].word.xbits[addr0] {
                destabilizers[i] *= g_q.clone();
            }
        }

        // Update destabilizer q to be the old stabilizer q (before replacement)
        destabilizers[q_idx] = g_q;

        // Finally, replace g_q by ±Z
        let stab_q = &mut stabilizers[q_idx];
        for i in 0..stab_q.n_qubits() {
            stab_q.word.xbits.set(i, false);
            stab_q.word.zbits.set(i, i == addr0);
        }
        stab_q.phase = if outcome { 2 } else { 0 };
    }
}

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64 { re: 1.0, im: 0.0 },  // +1
    Complex64 { re: 0.0, im: 1.0 },  // +i
    Complex64 { re: -1.0, im: 0.0 }, // -1
    Complex64 { re: 0.0, im: -1.0 }, // -i
];

pub fn symplectic_inner<I>(alpha: I, beta: I, n_qubits: usize) -> u32
where
    I: TableauIndex,
{
    let one = I::from(1u8);
    let zero = I::from(0u8);
    let mut parity = 0u32;
    for i in 0..n_qubits {
        if (alpha & beta) & (one << i) != zero {
            parity ^= 1;
        }
    }
    parity
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
        coefficients.unsafe_insert(I::from(0u8), complex_one);
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
    /// The function returns `(phase, gamma, lambda)`, where `gamma = (gamma_1, ..., gamma_n) as I`
    /// and `lambda = (lambda_1, ..., lambda_n) as I`. Note that gamma is equal to the shift of
    /// the index when branching (`beta` in Eq(4) of the SOFT paper).
    pub(crate) fn compute_decomposition(&self, addr0: usize, pauli: Pauli) -> (u8, I, I) {
        let n = self.n_qubits();

        // the actual decomposition, which we need to track the phase
        let mut p_word = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(n);
        p_word.set(addr0, pauli);

        // the bit strings defining the contributions
        let mut lambda = I::from(0u8);
        let mut gamma = I::from(0u8);

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
        let one = I::from(1u8);

        for (i, stab) in stabilizers.iter().enumerate() {
            if !destabilizers[i].word.anticommutes_at(addr0, pauli_bits) {
                // commutes
                continue;
            }

            // contributes, so set the corresponding bit in the lambda bit string
            // to 1
            lambda |= one << i;

            // destabilizer anti-commutes, so the stabilizer contributes
            // FIXME: don't need to clone here, just divide by phase twice
            let mut stab_inv = stab.clone();
            stab_inv.phase = (4 - stab.phase) % 4;
            p_word *= stab_inv;
        }

        // NOTE: destabilizers also commute with one another in a valid tableau
        // since the form a basis together with stabilizers
        for (i, destab) in destabilizers.iter().enumerate() {
            if !stabilizers[i].word.anticommutes_at(addr0, pauli_bits) {
                // commutes
                continue;
            }

            // contributes, so set the corresponding entry in the gamma bit string
            // to 1
            gamma |= one << i;

            // stabilizer anti-commutes, so the destabilizer contributes
            let mut destab_inv = destab.clone();
            destab_inv.phase = (4 - destab.phase) % 4;
            p_word *= destab_inv;
        }

        (p_word.phase, gamma, lambda)
    }

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    pub(crate) fn compute_phase(&self, lambda: I, basis_index: I, index_shift: I) -> u8 {
        // phase convention: 0: +1, 1: +i, 2: -1, 3: -i
        let one = I::from(1u8);
        let zero = I::from(0u8);
        let n = self.n_qubits();

        // contribution 1: each destabilizer D_i with basis_index[i]=1 that anticommutes
        // with P (lambda[i]=1) contributes a -1 sign; this is the symplectic inner product
        let mut phase = (2 * symplectic_inner(lambda, basis_index, n) as u8) % 4;

        // contribution 2: destabilizers that appear twice (basis_index[i]=1 and index_shift[i]=1)
        // contribute an extra -1 if their phase is imaginary
        let active = basis_index & index_shift;
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

    /// Keep only coefficients that correspond to the correct eigenvalue of a Z measurement.
    ///
    /// `outcome` controls which eigenspace to project onto:
    /// - `false` (outcome 0): keep terms where `phase == false` (commutes with Z, +1 eigenspace).
    ///   Use this when the tableau has already been updated to ±Z, because the new reference
    ///   state is already the correct eigenstate — `!phase` is always right in that frame.
    /// - `true` (outcome 1): keep terms where `phase == true` (anticommutes with Z, -1 eigenspace).
    ///   Use this when the tableau was *not* updated, so the reference state is unchanged and
    ///   we must explicitly select the -1 eigenspace.
    pub(crate) fn trim_coefficients_for_measurement(
        &mut self,
        addr0: usize,
        z_sign: bool,
        outcome: bool,
    ) {
        // TODO: more efficient update of coefficients in-place
        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let destabilizers = self.tableau.destabilizers();
        let n = self.n_qubits();
        let one = I::from(1u8);
        let zero = I::from(0u8);
        for (coeff, alpha) in old_coefficients.into_iter() {
            let mut phase = false; // false: +1 eigenspace of Z, true: -1 eigenspace

            // get the phase from the anti-commutation with the product over all destabilizers
            for i in 0..n {
                if alpha & (one << i) == zero {
                    // this index doesn't pick D_i
                    continue;
                }
                phase ^= destabilizers[i].word.xbits[addr0];
            }

            // (xi * k) == m
            if (phase ^ z_sign) == outcome {
                self.coefficients.add_or_insert(alpha, coeff);
            }
        }

        // renormalize
        self.coefficients.normalize();
    }

    pub(crate) fn branch_with_coefficients(
        &mut self,
        addr0: usize,
        pauli: Pauli,
        coefficient_factor: Complex<T::Coeff>,
        branch_factor: Complex<T::Coeff>,
    ) {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, index_shift, lambda) = self.compute_decomposition(addr0, pauli);

        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::new();
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ index_shift;

            // get the phase contributions from duplicate destabilizers
            // and anti-commuting through destabilizers
            let branch_phase_contribution = self.compute_phase(lambda, idx, index_shift);
            let branch_phase = (branch_phase_contribution + phase_decomp) % 4;

            let phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            let branch_coefficient = phase_factor * coeff.clone() * branch_factor.clone();
            let nonbranch_coefficient = coeff * coefficient_factor.clone();

            *new_coefficients
                .entry(branch_index)
                .or_insert(Complex::zero()) += branch_coefficient;
            *new_coefficients.entry(idx).or_insert(Complex::zero()) += nonbranch_coefficient;
        }

        let cutoff = Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        };
        for (idx, coeff) in new_coefficients {
            if coeff.abs() > cutoff.abs() {
                self.coefficients.unsafe_insert(idx, coeff);
            }
        }
    }

    pub(crate) fn compute_coefficients_after_pauli_apply(
        &self,
        coefficients: &mut C,
        addr0: usize,
        pauli: Pauli,
    ) {
        if self.is_lost[addr0] {
            return;
        }

        let (phase_decomp, index_shift, lambda) = self.compute_decomposition(addr0, pauli);

        let mut new_coefficients: HashMap<I, Complex<T::Coeff>> = HashMap::new();
        let old_coefficients = std::mem::replace(coefficients, C::new());
        for (coeff, idx) in old_coefficients.into_iter() {
            debug_assert!(
                !(coeff.re == T::Coeff::zero() && coeff.im == T::Coeff::zero()),
                "Coefficient should not be zero"
            );

            let branch_index = idx ^ index_shift;

            // get the phase contributions from duplicate destabilizers
            // and anti-commuting through destabilizers
            let branch_phase_contribution = self.compute_phase(lambda, idx, index_shift);
            let branch_phase = (branch_phase_contribution + phase_decomp) % 4;

            let phase_factor: Complex<T::Coeff> =
                COMPLEX_PHASE_CONVERSION[branch_phase as usize].into();

            let branch_coefficient = phase_factor * coeff.clone();

            *new_coefficients
                .entry(branch_index)
                .or_insert(Complex::zero()) += branch_coefficient;
        }

        let cutoff = Complex {
            re: self.coefficient_threshold.clone(),
            im: T::Coeff::zero(),
        };

        for (idx, coeff) in new_coefficients {
            if coeff.abs() > cutoff.abs() {
                coefficients.unsafe_insert(idx, coeff);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use bnum::BUint;

    use super::*;
    use crate::config::fxhash::ByteF64;
    use crate::traits::Clifford;

    type TestConfig = ByteF64<1>;
    type TestTableau = GeneralizedTableau<TestConfig>;
    type TestTableauBUint = GeneralizedTableau<TestConfig, BUint<1>>;

    #[test]
    fn test_compute_phase_z_2_single_qubit_plus_state() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);

        // After H: stabilizer = +X, destabilizer = +Z
        // shift = 1 (stabilizer has xbit[0]=true)
        // both phases should be 0
        let (decomp, gamma, lambda) = tab.compute_decomposition(0, Pauli::Z);
        let phase0 = decomp + tab.compute_phase(lambda, 0, gamma);
        assert_eq!(phase0, 0);
        let phase1 = decomp + tab.compute_phase(lambda, 1, gamma);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);

        let (decomp, gamma, lambda) = tab.compute_decomposition(0, Pauli::Z);
        let phase0 = decomp + tab.compute_phase(lambda, 0, gamma);
        assert_eq!(phase0, 0);
        let phase1 = decomp + tab.compute_phase(lambda, 1, gamma);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_mx_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.z(0);

        let (decomp, gamma, lambda) = tab.compute_decomposition(0, Pauli::Z);
        let phase0 = decomp + tab.compute_phase(lambda, 0, gamma);
        assert_eq!(phase0, 0);
        let phase1 = decomp + tab.compute_phase(lambda, 1, gamma);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer_2() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);
        tab.tableau.h(0);

        let (decomp, gamma, lambda) = tab.compute_decomposition(0, Pauli::Z);
        let phase0 = (decomp + tab.compute_phase(lambda, 0, gamma)) % 4;
        assert_eq!(phase0, 1);
        let phase1 = (decomp + tab.compute_phase(lambda, 1, gamma)) % 4;
        assert_eq!(phase1, 3);
    }

    #[test]
    fn test_index_type() {
        let mut tab: TestTableauBUint = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
    }
}
