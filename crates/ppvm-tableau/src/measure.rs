use super::data::symplectic_inner;
use crate::prelude::*;
use bitvec::view::BitView;
use fxhash::FxHashMap as HashMap;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use rand::RngExt;
use std::fmt::Debug;

impl<T: Config> Measure for Tableau<T>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    /// Measure qubit `addr0` in Z basis
    fn measure(&mut self, addr0: usize) -> bool {
        let q = self.find_z_anticommuting_stabilizer(addr0);
        match q {
            Some(q_idx) => {
                // Case a: random measurement outcome
                // At least one stabilizer anticommutes with Z_addr0

                // Generate random measurement outcome (50/50)
                let outcome = self.rng.random::<bool>();

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

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> LossyMeasure
    for GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
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
    fn measure(&mut self, addr0: usize) -> Option<bool> {
        if self.is_lost[addr0] {
            return None;
        }
        // NOTE: regardless of whether Z is a stabilizer, we need to compute
        // the probabilities, since the coefficients may make a Z stabilizer
        // state random, or a seemingly random one deterministic
        // the probabilities should just account for that

        // evaluate the action of Z on the state
        // i.e. shift + phase
        let mut z_overlap = Complex64::from(0.0);

        // TODO: this is O(n^2), but we know the probabilities are always real
        // however, whether the decomposition phase is imaginary or not tells us
        // whether we need to pick the real or imaginary part of the overlap
        // we still might be able to optimize here
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(addr0, Pauli::Z);

        // build a temporary lookup table for faster lookup in the loop
        let mut coeff_map: HashMap<I, Complex<T::Coeff>> = self
            .coefficients
            .clone()
            .into_iter()
            .map(|(v, i)| (i, v))
            .collect();
        // Compute the probabilities by computing the overlap <psi|Z|psi>
        // which is proportional to sum(alpha) conj(v_alpha) * v_(alpha + shift) * xi_(alpha)
        // NOTE: this could probably be optimized
        let odd_phase_mask = self.odd_phase_destabilizer_mask();
        for (&idx, coeff) in &coeff_map {
            let branch_index = idx ^ stab_anticomm_bits; // stab_anticomm_bits is the index shift
            let phase = (phase_decomp
                + self.compute_phase_with_mask(
                    destab_anticomm_bits,
                    idx,
                    stab_anticomm_bits,
                    odd_phase_mask,
                ))
                % 4;
            let complex_phase: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[phase as usize].into();
            let Some(coeff_branch) = coeff_map.get(&branch_index).copied() else {
                continue;
            };
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

        let outcome = self.tableau.rng.random::<f64>() < prob_1;

        if stab_anticomm_bits != I::zero() {
            // Case a: Z is not a stabilizer

            // first anti-commuting stabilizer is just the first nonzero bit in stab_anticomm_bits
            let q_idx = stab_anticomm_bits.trailing_zeros() as usize;

            // In this case, we cannot simply trim the coefficients (though some
            // might be smaller than the threshold)

            // coefficient algorithm from T.J. Yoder, adapted for state vectors
            // see Algorithm 2 in https://www.scottaaronson.com/showcase2/report/ted-yoder.pdf

            // k = lowest set bit of stab_anticomm_bits (same position as q_idx)
            let one = I::one();
            let zero = I::zero();
            let k = one << q_idx;

            let alpha = if outcome {
                (phase_decomp + 2) % 4
            } else {
                phase_decomp
            };

            // Drain B entries (k-bit=1) into their A counterparts (k-bit=0) in-place,
            // avoiding O(M^2) add_or_insert on Vec
            let b_keys: Vec<I> = coeff_map
                .keys()
                .filter(|idx| (**idx & k) != zero)
                .cloned()
                .collect();
            for idx in b_keys {
                let coeff = coeff_map.remove(&idx).unwrap();
                // q = phase_decomp * (-1).pow(symplectic_inner(*idx, destab_anticomm_bits)) * q;
                let symp_inner = symplectic_inner(idx, destab_anticomm_bits);
                let phase_idx =
                    ((alpha as i32 + if symp_inner % 2 == 1 { 2 } else { 0 }) % 4) as usize;
                let q: Complex<T::Coeff> = COMPLEX_PHASE_CONVERSION[phase_idx].into();
                *coeff_map
                    .entry(idx ^ stab_anticomm_bits) // stab_anticomm_bits is the index shift
                    .or_insert(Complex::zero()) += q * coeff;
            }

            // Keep entries where |c|/norm > threshold.
            let norm_sqr = coeff_map
                .values()
                .fold(T::Coeff::zero(), |acc, c: &Complex<T::Coeff>| {
                    acc + c.norm_sqr()
                });

            let cutoff_sq = self.coefficient_threshold.clone() * self.coefficient_threshold.clone();
            let threshold = cutoff_sq.to_f64().unwrap_or(0.0) * norm_sqr.to_f64().unwrap_or(0.0);
            self.coefficients = C::new();
            for (idx, coeff) in coeff_map {
                if coeff.norm_sqr() > threshold {
                    self.coefficients.unsafe_insert(idx, coeff);
                }
            }

            // normalize again, since dropping coefficients may reduce the norm
            self.coefficients.normalize();

            // update the tableau, coefficients can be updated independently
            self.tableau
                .update_tableau_according_to_outcome(addr0, q_idx, outcome);
        } else {
            // Case b: +Z or -Z already is a stabilizer; we just need
            // to trim the coefficients accordingly; tableau remains unchanged

            // Applying the projector to a basis state, we have three phases:
            // 1. The actual measurement outcome (k)
            // 2. The sign from whether +Z or -Z is a stabilizer (m) - can get that from the decomposition
            // 3. Contribution from commuting Z_addr0 through the destabilizers (xi)
            // Only coefficients where m*k*xi == 1 are kept

            // 2. get the sign
            debug_assert!(
                phase_decomp == 0 || phase_decomp == 2,
                "Measurement result cannot be imaginary!"
            );
            let z_sign = phase_decomp == 2;

            // 3. check the anticommutation -- combine with coefficient update
            self.trim_coefficients_for_measurement(destab_anticomm_bits, outcome, z_sign);
        }
        Some(outcome)
    }
}
