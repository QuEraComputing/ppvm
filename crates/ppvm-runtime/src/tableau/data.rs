use std::{
    collections::HashMap,
    fmt::Debug,
    marker::PhantomData,
    ops::{BitAnd, Shl},
};

use super::sparsevec::SparseVector;
use crate::{char::Pauli, config::Config};
use crate::{phase::PhasedPauliWord, tableau::traits::TableauIndex};
use num::{
    One, Zero,
    complex::{Complex, Complex64, ComplexFloat},
};

#[derive(Clone, Debug)]
pub struct Tableau<T: Config> {
    pub n_qubits: usize,
    /// Destabilizer / Stabilizer tableau
    /// * Entries 0..n are the destabilizers
    /// * Entries n..2n are the stabilizers
    pub data: Vec<PhasedPauliWord<T::Storage, T::BuildHasher>>,
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

        Self { n_qubits, data }
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
        // the stabilizer can be computed as the product of all destabilizers
        // it anticommutes with; we do this and then check the phase to determine if it's Z or -Z
        // NOTE: we can just skip building the actual Pauli string since we only need the phase
        let destabilizers = self.destabilizers();
        let stabilizers = self.stabilizers();
        let mut phase = 0;
        for (i, destab) in destabilizers.iter().enumerate() {
            if destab.word.xbits[addr0] {
                phase = (phase + stabilizers[i].phase) % 4;
            }
        }

        // phase >= 2 means -Z eigenvalue → outcome |1⟩ (true)
        phase >= 2
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
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
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
            coefficients: coefficients,
            is_lost: vec![false; n_qubits],
            coefficient_threshold,
            _index_phantom: PhantomData,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits
    }

    // helper functions

    pub(crate) fn compute_decomposition_phase(&self, addr0: usize, pauli: Pauli) -> u8 {
        // NOTE: this is O(n ^ 2); can we improve it since we only need the phase?

        // P_addr0 = phase * prod(d_k ^ gamma_k) * prod(s_l ^ lambda_l)
        // where: gamma_k == 1 iff {P_addr0, s_k} = 0
        // lambda_l == 1 iff {P_addr0, d_l} = 0
        // Lemma 5. from T. J. Yoder (2012)
        // now, we just need to invert the expression to compute the phase
        let n = self.n_qubits();
        let mut p_word = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(n);
        p_word.set(addr0, pauli);

        debug_assert_ne!(pauli, Pauli::I);
        let pauli_bits = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
        };

        let stabilizers = self.tableau.stabilizers();
        let destabilizers = self.tableau.destabilizers();

        for (i, stab) in stabilizers.iter().enumerate() {
            if !destabilizers[i].word.anticommutes_at(addr0, pauli_bits) {
                // commutes
                continue;
            }

            // destabilizer anti-commutes, so the stabilizer contributes
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

            // stabilizer anti-commutes, so the destabilizer contributes
            let mut destab_inv = destab.clone();
            destab_inv.phase = (4 - destab.phase) % 4;
            p_word *= destab_inv;
        }

        p_word.phase
    }

    /// Compute the index shift when applying a Pauli
    pub(crate) fn compute_shift(&self, addr0: usize, pauli: (bool, bool)) -> I {
        // NOTE: we use LSB ordering
        let mut shift = I::from(0u8);
        let one = I::from(1u8);
        debug_assert!(pauli.0 || pauli.1); // should never be called on Pauli::I
        for (i, stab) in self.tableau.stabilizers().iter().enumerate() {
            if stab.word.anticommutes_at(addr0, pauli) {
                shift |= one << i;
            }
        }
        shift
    }

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    pub(crate) fn compute_phase(
        &self,
        addr0: usize,
        pauli: (bool, bool),
        basis_index: I,
        index_shift: I,
    ) -> u8 {
        // phase convention: 0: +1, 1: +i, 2: -1, 3: -i
        let mut phase = 0u8;
        let one = I::from(1u8);
        let zero = I::from(0u8);
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if basis_index & (one << i) == zero {
                // NOTE: LSB ordering; has to be consistent with shift computation
                continue;
            }

            if destab.word.anticommutes_at(addr0, pauli) {
                // We have an xbit set, so we anticommute, leading to a -1 sign
                phase = (phase + 2) % 4;
            }

            if index_shift & (I::from(1u8) << i) == I::from(0u8) {
                continue;
            }

            // this particular destabilizer occurs twice: once from the P decomposition
            // this is given by the index shift, since the corresponding bit in the shift
            // is only 1 if P anti-commutes with the stabilizer, meaning its decomposition
            // features the destabilizer here
            if destab.phase % 2 != 0 {
                // phase of the destabilizer is ~i, so its square gives another -1
                phase = (phase + 2) % 4;
            }
        }
        phase
    }

    /// Keep only coefficients that correspond to the correct eigenvalue of
    /// a Z measurement
    /// NOTE: this is called AFTER the tableau has been updated
    pub(crate) fn trim_coefficients_for_measurement(&mut self, addr0: usize) {
        // TODO: more efficient update of coefficients in-place
        let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
        let destabilizers = self.tableau.destabilizers();
        let n = self.n_qubits();
        let one = I::from(1u8);
        let zero = I::from(0u8);
        for (coeff, alpha) in old_coefficients.into_iter() {
            let mut phase = false; // false: 1, true: -1

            // get the phase from the anti-commutation with the product over all destabilizers
            for i in 0..n {
                if alpha & (one << i) == zero {
                    // this index doesn't pick D_i
                    continue;
                }
                phase ^= destabilizers[i].word.xbits[addr0];
            }

            // NOTE: if the term accumulates a phase, then the projector
            // (I + P) |b_alpha> ~ (I - P) |psi_s>, where P is +Z or -Z
            // since P is a stabilizer in the updated tableau, any term
            // where a negative phase is accumulated zeros out
            if !phase {
                // keep term
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

        let pauli_booleans = match pauli {
            Pauli::I => (false, false),
            Pauli::X => (true, false),
            Pauli::Y => (true, true),
            Pauli::Z => (false, true),
        };

        let index_shift = self.compute_shift(addr0, pauli_booleans);
        let phase_decomp = self.compute_decomposition_phase(addr0, pauli);

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
            let branch_phase_contribution =
                self.compute_phase(addr0, pauli_booleans, idx, index_shift);
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
        let phase0 = tab.compute_decomposition_phase(0, Pauli::Z)
            + tab.compute_phase(0, (false, true), 0, 1);
        assert_eq!(phase0, 0);
        let phase1 = tab.compute_decomposition_phase(0, Pauli::Z)
            + tab.compute_phase(0, (false, true), 1, 1);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);

        let shift = tab.compute_shift(0, (false, true));
        let decomp = tab.compute_decomposition_phase(0, Pauli::Z);
        let phase0 = decomp + tab.compute_phase(0, (false, true), 0, shift);
        assert_eq!(phase0, 0);
        let phase1 = decomp + tab.compute_phase(0, (false, true), 1, shift);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_mx_stabilizer() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.z(0);

        let shift = tab.compute_shift(0, (false, true));
        let decomp = tab.compute_decomposition_phase(0, Pauli::Z);
        let phase0 = decomp + tab.compute_phase(0, (false, true), 0, shift);
        assert_eq!(phase0, 0);
        let phase1 = decomp + tab.compute_phase(0, (false, true), 1, shift);
        assert_eq!(phase1, 0);
    }

    #[test]
    fn test_compute_phase_z_2_y_stabilizer_2() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
        tab.tableau.s(0);
        tab.tableau.h(0);

        let shift = tab.compute_shift(0, (false, true));
        let decomp = tab.compute_decomposition_phase(0, Pauli::Z);
        let phase0 = (decomp + tab.compute_phase(0, (false, true), 0, shift)) % 4;
        assert_eq!(phase0, 1);
        let phase1 = (decomp + tab.compute_phase(0, (false, true), 1, shift)) % 4;
        assert_eq!(phase1, 3);
    }

    #[test]
    fn test_index_type() {
        let mut tab: TestTableauBUint = GeneralizedTableau::new(1, 1e-12);
        tab.tableau.h(0);
    }
}
