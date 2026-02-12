use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};

impl<T: Config> Measure for Tableau<T> {
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

                self.get_deterministic_outcome(addr0)
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

impl<T: Config, C: SparseVector<Complex<T::Coeff>>> Measure for GeneralizedTableau<T, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>> + std::fmt::Debug,
    T::Coeff: One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
{
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.tableau.find_anticommuting_stabilizer(addr0);

        let outcome = match q {
            Some(q_idx) => {
                // Case a: Random outcome
                // NOTE: this may also be a deterministic outcome, with branches present

                // evaluate the action of Z on the state
                // i.e. shift + phase
                let shift = self.compute_shift_z(addr0);
                let mut z_overlap = Complex64::from(0.0);

                // NOTE: for the overlap, need to get the phase
                // of Z in the stabilizer state
                let mut z_phase = 0u8;
                // for stab in self.tableau.stabilizers().iter() {
                //     if stab.word.xbits[addr0] {
                //         let has_z = stab.word.zbits[addr0];

                //         // TODO: check whether we need to account for stabilizer z_phase
                //         // z_phase = (z_phase + stab.z_phase) % 4;

                //         if has_z {
                //             // Y operator contributes a z_phase of -i
                //             z_phase = (z_phase + 3) % 4;
                //         } else {
                //             // X operator contributes a z_phase of -1
                //             z_phase = (z_phase + 2) % 4;
                //         }
                //     }
                // }
                // Compute the probabilities by computing the overlap <psi|Z|psi>
                // which is proportional to sum(alpha) conj(v_alpha) * v_(alpha + shift) * xi_(alpha)
                for (coeff, idx) in self.coefficients.clone().into_iter() {
                    let branch_index = idx ^ shift;
                    // TODO: double-check the phase, this might need to be computed with the branch_index
                    let phase = (z_phase + self.compute_phase_z_2(addr0, idx)) % 4;
                    let complex_phase: Complex<T::Coeff> =
                        COMPLEX_PHASE_CONVERSION[phase as usize].into();
                    let coeff_branch = self.coefficients.get(&branch_index);
                    let overlap = complex_phase.conj() * coeff.conj() * coeff_branch;
                    z_overlap.re += overlap.re.to_f64().unwrap_or(0.0);
                    z_overlap.im += overlap.im.to_f64().unwrap_or(0.0);
                }

                debug_assert!(
                    z_overlap.im.abs() < 1e-6,
                    "Overlap should be real, got {}",
                    z_overlap
                );

                // TODO: directly compute one of these probs above and skip the other
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
                self.tableau
                    .update_tableau_according_to_outcome(addr0, q_idx, outcome);

                outcome
            }

            None => {
                // Case b: deterministic outcome
                self.tableau.get_deterministic_outcome(addr0)
            }
        };

        self.trim_coefficients_for_measurement(addr0);

        outcome
    }
}
