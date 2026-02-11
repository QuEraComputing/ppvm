use super::data::{GeneralizedTableau, Tableau};
use super::traits::{GeneralizedTableauTGate, Measure};
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};

impl<const N: usize, T: Config> Measure for Tableau<N, T> {
    fn find_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        // Find first stabilizer that anticommutes with Z_addr0
        let mut q = None;
        for (i, stab) in self.stabilizers.iter().enumerate() {
            if stab.word.xbits[addr0] {
                // X or Y anticommutes with Z
                q = Some(i);
                break;
            }
        }
        q
    }

    fn update_tableau_according_to_outcome(&mut self, addr0: usize, q_idx: usize, outcome: bool) {
        // Check if there are other stabilizers that anticommute with Z_addr0
        // If so, replace with g_j = g_j * g_q
        for i in 0..N {
            if i == q_idx {
                continue;
            }
            if self.stabilizers[i].word.xbits[addr0] {
                // Stabilizer i also anticommutes, so multiply by g_q to eliminate
                let g_q = self.stabilizers[q_idx].clone();
                self.stabilizers[i] *= g_q;
            }
            if self.destabilizers[i].word.xbits[addr0] {
                let g_q = self.stabilizers[q_idx].clone();
                self.destabilizers[i] *= g_q;
            }
        }

        // Update destabilizer q to be the old stabilizer q (before replacement)
        self.destabilizers[q_idx] = self.stabilizers[q_idx].clone();

        // Finally, replace g_q by \pm Z
        for i in 0..self.stabilizers[q_idx].n_qubits() {
            // set the q_idx stabilizer to the Pauli string IIZIII...I
            self.stabilizers[q_idx].word.xbits.set(i, false);
            self.stabilizers[q_idx].word.zbits.set(i, i == addr0);
        }

        // Set phase depending on outcome
        self.stabilizers[q_idx].phase = if outcome { 2 } else { 0 };
    }

    /// Measure qubit `addr0` in Z basis
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_anticommuting_stabilizer(addr0);
        match q {
            Some(q_idx) => {
                // Case a: random measurement outcome
                // At least one stabilizer anticommutes with Z_addr0

                // Generate random measurement outcome (50/50)
                let outcome = rand::random::<bool>();

                self.update_tableau_according_to_outcome(addr0, q_idx, outcome);

                outcome
            }
            None => {
                // Case b: deterministic measurement outcome

                // find the outcome: either Z_addr0 or -Z_addr0 is a stabilizer
                // the stabilizer can be computed as the product of all destabilizers
                // it anticommutes with; we do this and then check the phase to determine if it's Z or -Z
                // NOTE: we can just skip building the actual Pauli string since we only need the phase
                let mut phase = 0;
                for (i, destab) in self.destabilizers.iter().enumerate() {
                    if destab.word.xbits[addr0] {
                        phase = (phase + self.stabilizers[i].phase) % 4;
                    }
                }

                // phase >= 2 means -Z eigenvalue → outcome |1⟩ (true)
                phase >= 2
            }
        }
    }
}

const COMPLEX_PHASE_CONVERSION: [Complex64; 4] = [
    Complex64::new(1.0, 0.0),  // +1
    Complex64::new(0.0, 1.0),  // +i
    Complex64::new(-1.0, 0.0), // -1
    Complex64::new(0.0, -1.0), // -i
];

impl<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> Measure
    for GeneralizedTableau<N, T, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>> + std::fmt::Debug,
    T::Coeff: One + Zero + Clone + num::Num + ToPrimitive,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
{
    fn find_anticommuting_stabilizer(&self, addr0: usize) -> Option<usize> {
        self.tableau.find_anticommuting_stabilizer(addr0)
    }

    fn update_tableau_according_to_outcome(&mut self, addr0: usize, q_idx: usize, outcome: bool) {
        self.tableau
            .update_tableau_according_to_outcome(addr0, q_idx, outcome);
    }

    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_anticommuting_stabilizer(addr0);

        match q {
            Some(q_idx) => {
                // Case a: Random outcome

                // evaluate the action of Z on the state
                // i.e. shift + phase
                let shift = self.compute_shift_z(addr0);
                let mut z_overlap = Complex64::from(0.0);
                // Compute the probabilities by computing the overlap <psi|Z|psi>
                // which is proportional to sum(alpha) conj(v_alpha) * v_(alpha + shift) * xi_(alpha)
                for (coeff, idx) in self.coefficients.clone().into_iter() {
                    let branch_index = idx ^ shift;
                    // TODO: double-check the phase, this might need to be computed with the branch_index
                    let phase = self.compute_phase_z(addr0, idx);
                    let complex_phase: Complex<T::Coeff> =
                        COMPLEX_PHASE_CONVERSION[phase as usize].into();
                    // let eigenvalue = phase >= 2;
                    let coeff_branch = self.coefficients.get(&branch_index);
                    let overlap = complex_phase * coeff.conj() * coeff_branch;
                    z_overlap.re += overlap.re.to_f64().unwrap_or(0.0);
                    z_overlap.im += overlap.im.to_f64().unwrap_or(0.0);
                    // let prob = (complex_phase * coeff.conj() * coeff_branch)
                    //     .re
                    //     .to_f64()
                    //     .unwrap_or(0.0);
                    // if eigenvalue {
                    //     prob_minus += prob;
                    // } else {
                    //     prob_plus += prob;
                    // }
                }

                debug_assert!(
                    z_overlap.im.abs() < 1e-6,
                    "Overlap should be real, got {}",
                    z_overlap
                );

                let prob_0 = 0.5 + 0.5 * z_overlap.re;
                let prob_1 = 0.5 - 0.5 * z_overlap.re;

                debug_assert!(
                    (prob_0 + prob_1 - 1.0).abs() < 1e-6,
                    "Probabilities should sum to 1, got {} + {} = {}",
                    prob_0,
                    prob_1,
                    prob_0 + prob_1
                );

                let outcome = rand::random::<f64>() < prob_1;

                // update the coefficients so only the ones with the correct outcome are kept
                let mut new_coefficients = C::new();
                for (coeff, idx) in self.coefficients.clone().into_iter() {
                    let branch_index = idx ^ shift;
                    let phase = self.compute_phase_z(addr0, idx);
                    let mut complex64_phase = 0.5 * COMPLEX_PHASE_CONVERSION[phase as usize];
                    if outcome {
                        complex64_phase *= -1.0;
                    }
                    let complex_phase: Complex<T::Coeff> = complex64_phase.into();

                    let value = complex_phase * coeff;
                    new_coefficients.add_or_insert(branch_index, value);
                }

                println!("{:?}", new_coefficients);

                for (_coeff, idx) in self.coefficients.clone().into_iter() {
                    self.coefficients
                        .mul_element_by(idx, Complex64::from(0.5).into());
                }

                for (new_coeff, idx) in new_coefficients.clone().into_iter() {
                    self.coefficients.add_or_insert(idx, new_coeff);
                }

                self.update_tableau_according_to_outcome(addr0, q_idx, outcome);

                outcome
            }

            None => {
                // Case b: deterministic outcome

                // TODO: more efficient update of coefficients in-place
                let old_coefficients = std::mem::replace(&mut self.coefficients, C::new());
                for (coeff, alpha) in old_coefficients.into_iter() {
                    let mut phase = false; // false: 1, true: -1

                    // get the phase from the anti-commutation with the product over all destabilizers
                    for i in 0..N {
                        if alpha & (1 << i) == 0 {
                            // this index doesn't pick D_i
                            continue;
                        }
                        phase ^= self.tableau.destabilizers[i].word.xbits[addr0];
                    }

                    if !phase {
                        // keep term
                        self.coefficients.add_or_insert(alpha, coeff);
                    } // else drop it, since it would flip the sign in (1 + P)|b_alpha> regardless of whether P is +Z or -Z
                }

                // renormalize
                self.coefficients.normalize();

                // TODO: for this bit, move the logic into a new function and share it for tableau
                // find the outcome of the deterministic measurement
                let mut phase = 0;
                for (i, destab) in self.tableau.destabilizers.iter().enumerate() {
                    if destab.word.xbits[addr0] {
                        phase = (phase + self.tableau.stabilizers[i].phase) % 4;
                    }
                }

                // phase >= 2 means -Z eigenvalue → outcome |1⟩ (true)
                let outcome = phase >= 2;
                outcome
            }
        }
    }
}
