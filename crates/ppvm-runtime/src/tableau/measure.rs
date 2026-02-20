use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{BitAnd, BitOrAssign, BitXor, Shl};

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

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Measure for GeneralizedTableau<T, I, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One + Zero + Clone + num::Num + ToPrimitive + std::fmt::Debug,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    I: PartialEq
        + Eq
        + Hash
        + Copy
        + From<u8>
        + Shl<usize>
        + BitOrAssign<<I as Shl<usize>>::Output>
        + BitAnd<<I as Shl<usize>>::Output, Output = I>
        + BitXor<Output = I>,
    <I as BitAnd<<I as Shl<usize>>::Output>>::Output: PartialEq<I>,
{
    fn measure(&mut self, addr0: usize) -> bool {
        // NOTE: regardless of whether Z is a stabilizer, we need to compute
        // the probabilities, since the coefficients may make a Z stabilizer
        // state random, or a seemingly random one deterministic
        // the probabilities should just account for that

        // evaluate the action of Z on the state
        // i.e. shift + phase
        let shift = self.compute_shift(addr0, (false, true));
        let mut z_overlap = Complex64::from(0.0);

        // TODO: this is O(n^2), but we know the probabilities are always real
        // however, whether the decomposition phase is imaginary or not tells us
        // whether we need to pick the real or imaginary part of the overlap
        // we still might be able to optimize here
        let phase_decomp = self.compute_z_decomposition_phase(addr0);

        // build a temporary lookup table for faster lookup in the loop
        let coeff_map: HashMap<I, Complex<T::Coeff>> = self
            .coefficients
            .clone()
            .into_iter()
            .map(|(v, i)| (i, v))
            .collect();
        // Compute the probabilities by computing the overlap <psi|Z|psi>
        // which is proportional to sum(alpha) conj(v_alpha) * v_(alpha + shift) * xi_(alpha)
        // NOTE: this could probably be optimized
        for (&idx, coeff) in &coeff_map {
            let branch_index = idx ^ shift;
            let phase = (phase_decomp + self.compute_phase_z(addr0, idx, shift)) % 4;
            let complex_phase: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[phase as usize].into();
            let coeff_branch = coeff_map
                .get(&branch_index)
                .cloned()
                .unwrap_or(Complex::zero());
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

        // Now, we may need to update the tableau if Z is not a stabilizer
        let q = self.tableau.find_anticommuting_stabilizer(addr0);

        match q {
            Some(q_idx) => {
                // Case a: Z is not a stabilizer
                // tableau needs to be updated
                self.tableau
                    .update_tableau_according_to_outcome(addr0, q_idx, outcome);
            }

            None => {
                // Case b: tableau does not change since Z is a stabilizer; no-op
                // NOTE: the outcome here may still be random and need not match
                // the deterministic outcome from the tableau
            }
        };

        self.trim_coefficients_for_measurement(addr0);

        outcome
    }
}
