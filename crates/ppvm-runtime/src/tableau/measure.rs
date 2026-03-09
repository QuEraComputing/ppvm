use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::tableau::sparsevec::SparseVector;
use crate::tableau::traits::TableauIndex;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use std::collections::HashMap;
use std::fmt::Debug;

impl<T: Config> Measure for Tableau<T> {
    /// Measure qubit `addr0` in Z basis
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_z_anticommuting_stabilizer(addr0);
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
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    I: TableauIndex + Debug,
{
    fn measure(&mut self, addr0: usize) -> bool {
        if self.is_lost[addr0] {
            return false;
        }
        // NOTE: regardless of whether Z is a stabilizer, we need to compute
        // the probabilities, since the coefficients may make a Z stabilizer
        // state random, or a seemingly random one deterministic
        // the probabilities should just account for that

        // TODO: we can optimize this by looking at which states get eliminated
        // first and then computing the probabilities as the norm from there
        // this skips the O(n ^ 2) evaluation of <Z>

        // evaluate the action of Z on the state
        // i.e. shift + phase
        let shift = self.compute_shift(addr0, (false, true));
        let mut z_overlap = Complex64::from(0.0);

        // TODO: this is O(n^2), but we know the probabilities are always real
        // however, whether the decomposition phase is imaginary or not tells us
        // whether we need to pick the real or imaginary part of the overlap
        // we still might be able to optimize here
        let phase_decomp = self.compute_decomposition_phase(addr0, crate::char::Pauli::Z);

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
            let phase = (phase_decomp + self.compute_phase(addr0, (false, true), idx, shift)) % 4;
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

        let q = self.tableau.find_z_anticommuting_stabilizer(addr0);

        match q {
            Some(q_idx) => {
                debug_assert_ne!(shift, I::from(0u8), "Shift 0, but Z is not a stabilizer!");
                // Case a: Z is not a stabilizer

                // In this case, we cannot simply trim the coefficients (though some
                // might be smaller than the threshold)

                // coefficient algorithm from T.J. Yoder, adapted for state vectors
                // see Algorithm 2 in https://www.scottaaronson.com/showcase2/report/ted-yoder.pdf

                // get k: bit string with a single 1 entry at the position
                // of the first 1 in shift
                let mut k = I::from(0u8);
                let one = I::from(1u8);
                let zero = I::from(0u8);
                for i in 0..self.n_qubits() {
                    if shift & (one << i) != zero {
                        k = one << i;
                        break;
                    }
                }

                // Find the stabilizer index of decomposing Z_addr0 into
                // stabilizers and destabilizers
                // TODO: combine this with getting the decomposition phase
                let mut c = I::from(0u8);
                let destabilizers = self.tableau.destabilizers();
                for (i, destab) in destabilizers.iter().enumerate() {
                    if destab.word.xbits[addr0] {
                        // anti-commuting destabilizer
                        // meaning the stabilizer contributes to the decomp
                        c |= one << i;
                    }
                }

                let alpha = if outcome {
                    (phase_decomp + 2) % 4
                } else {
                    phase_decomp
                };

                // TODO: hashmap for assigning new coefficients
                let mut new_coefficients = C::new();
                for (idx, coeff) in &coeff_map {
                    let mut x = idx.clone();
                    let mut q: Complex<T::Coeff> = Complex::one();
                    if (*idx & k) != zero {
                        // q = phase_decomp * (-1).pow(symplectic_inner(*idx, c)) * q;
                        let symp_inner = {
                            let mut parity = 0u32;
                            for i in 0..self.n_qubits() {
                                if (*idx & c) & (one << i) != zero {
                                    parity ^= 1;
                                }
                            }
                            parity
                        };
                        let phase_idx =
                            ((alpha as i32 + if symp_inner % 2 == 1 { 2 } else { 0 }) % 4) as usize;
                        q = COMPLEX_PHASE_CONVERSION[phase_idx].into();
                        x = *idx ^ shift;
                    }
                    let half: Complex<T::Coeff> = Complex64::new(0.5, 0.0).into();
                    let new_coeff = q * *coeff * half;
                    new_coefficients.add_or_insert(x, new_coeff);
                }

                new_coefficients.normalize();

                self.coefficients = new_coefficients;

                // update the tableau, coefficients can be updated independently
                self.tableau
                    .update_tableau_according_to_outcome(addr0, q_idx, outcome);
            }
            None => {
                debug_assert_eq!(shift, I::from(0u8), "Shift !=0 but Z is a stabilizer!");
                // Case b: +Z or -Z already is a stabilizer; we just need
                // to trim the coefficients accordingly; tableau remains unchanged

                // Applying the projector to a basis state, we have three phases:
                // 1. The actual measurement outcome (k)
                // 2. The sign from whether +Z or -Z is a stabilizer (m)
                // 3. Contribution from commuting Z_addr0 through the destabilizers (xi)
                // Only coefficients where m*k*xi == 1 are kept

                // 2. get the sign
                let z_sign = self.tableau.get_deterministic_outcome(addr0);

                // 3. check the anticommutation -- combine with coefficient update
                self.trim_coefficients_for_measurement(addr0, outcome, z_sign);
            }
        };

        outcome
    }
}
