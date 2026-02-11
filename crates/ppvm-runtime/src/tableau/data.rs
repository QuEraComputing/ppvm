use super::sparsevec::SparseVector;
use crate::config::Config;
use crate::phase::PhasedPauliWord;
use num::{
    One, Zero,
    complex::{Complex, Complex64},
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
    pub(crate) fn find_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        // Find first stabilizer that anticommutes with Z_addr0
        let mut q = None;
        for (i, stab) in self.stabilizers().iter().enumerate() {
            if stab.word.xbits[addr0] {
                // X or Y anticommutes with Z
                q = Some(i);
                break;
            }
        }
        q
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

// TODO: builder
pub struct GeneralizedTableau<T: Config, C: SparseVector<Complex<T::Coeff>>> {
    pub tableau: Tableau<T>,
    pub coefficients: C,
    pub is_lost: Vec<bool>,
    pub coefficient_threshold: T::Coeff,
}

impl<T: Config, C: SparseVector<Complex<T::Coeff>>> GeneralizedTableau<T, C>
where
    T::Coeff: One + Zero + Clone,
{
    pub fn new(n_qubits: usize, coefficient_threshold: T::Coeff) -> Self {
        let mut coefficients = C::new();
        let complex_one = Complex {
            re: T::Coeff::one(),
            im: T::Coeff::zero(),
        };
        coefficients.unsafe_insert(0, complex_one);
        Self {
            tableau: Tableau::new(n_qubits),
            coefficients: coefficients,
            is_lost: vec![false; n_qubits],
            coefficient_threshold,
        }
    }

    pub fn n_qubits(&self) -> usize {
        self.tableau.n_qubits
    }

    // helper functions

    /// Compute the index shift when applying a Z Pauli
    pub(crate) fn compute_shift_z(&self, addr0: usize) -> usize {
        // NOTE: we use LSB ordering
        let mut shift = 0usize;
        for (i, stab) in self.tableau.stabilizers().iter().enumerate() {
            shift |= (stab.word.xbits[addr0] as usize) << i;
        }
        shift
    }

    /// every basis index is a bit string alpha defining the basis state
    /// the phase when applying a Pauli is the product of all destabilizer phases
    /// and the phase contributions from the commutation relations
    /// we need to check every destabilizer where the basis index has a 1 bit.
    pub(crate) fn compute_phase_z(&self, addr0: usize, basis_index: usize) -> u8 {
        // phase convention: 0: +1, 1: +i, 2: -1, 3: -i
        let mut phase = 0u8;
        for (i, destab) in self.tableau.destabilizers().iter().enumerate() {
            if basis_index & (1 << i) == 0 {
                // NOTE: LSB ordering; has to be consistent with shift computation
                continue;
            }

            let has_x = destab.word.xbits[addr0];
            let has_z = destab.word.zbits[addr0];

            // need to account for destabilizer phase
            phase = (phase + destab.phase) % 4;

            if has_x && has_z {
                // Y operator contributes a phase of -i
                phase = (phase + 3) % 4;
            } else if has_x {
                // X operator contributes a phase of -1
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
        for (coeff, alpha) in old_coefficients.into_iter() {
            let mut phase = false; // false: 1, true: -1

            // get the phase from the anti-commutation with the product over all destabilizers
            for i in 0..n {
                if alpha & (1 << i) == 0 {
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
}
