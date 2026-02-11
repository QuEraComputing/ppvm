use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
use crate::phase::PhasedPauliWord;
use crate::tableau::sparsevec::SparseVector;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, Zero};

// NOTE: this trait impl was 100% vibe-coded, please don't judge
impl<const N: usize, T: Config> Measure for Tableau<N, T> {
    /// Measure qubit `addr0` in the computational (Z) basis.
    ///
    /// Returns the measurement outcome: false for |0⟩, true for |1⟩.
    ///
    /// This implements the standard stabilizer measurement algorithm:
    /// - If any stabilizer anticommutes with Z_i (has X or Y), the outcome is random
    /// - Otherwise, the outcome is deterministic based on stabilizer phases
    fn measure(&mut self, addr0: usize) -> bool {
        // Step 1: Find first stabilizer that anticommutes with Z_addr0
        // (i.e., has X or Y at position addr0)
        let mut p = None;
        for (i, stab) in self.stabilizers.iter().enumerate() {
            if stab.word.xbits[addr0] {
                // Has X or Y (X bit is set)
                p = Some(i);
                break;
            }
        }

        match p {
            Some(p_idx) => {
                // RANDOM MEASUREMENT CASE:
                // At least one stabilizer anticommutes with Z_addr0

                // Generate random measurement outcome
                let outcome = rand::random::<bool>();

                // Perform Gaussian elimination: multiply other anticommuting stabilizers by g_p
                // This ensures only g_p has X or Y at addr0
                for i in 0..N {
                    if i != p_idx && self.stabilizers[i].word.xbits[addr0] {
                        // Stabilizer i also anticommutes, so multiply by g_p to eliminate
                        let g_p = self.stabilizers[p_idx].clone();
                        self.stabilizers[i] *= g_p;
                    }
                }

                // Also update destabilizers that anticommute with Z_addr0
                // (skip p_idx since it will be overwritten)
                for i in 0..N {
                    if i != p_idx && self.destabilizers[i].word.xbits[addr0] {
                        let g_p = self.stabilizers[p_idx].clone();
                        self.destabilizers[i] *= g_p;
                    }
                }

                // Step (b): Copy old stabilizer to destabilizer
                // (must happen before overwriting the stabilizer)
                self.destabilizers[p_idx] = self.stabilizers[p_idx].clone();

                // Step (c): Replace stabilizer p with ±Z_addr0 based on outcome
                // outcome = false (|0⟩) → +Z, outcome = true (|1⟩) → -Z
                for i in 0..self.stabilizers[p_idx].n_qubits() {
                    self.stabilizers[p_idx].word.xbits.set(i, false);
                    self.stabilizers[p_idx].word.zbits.set(i, i == addr0);
                }
                self.stabilizers[p_idx].phase = if outcome { 2 } else { 0 };

                outcome
            }
            None => {
                // DETERMINISTIC MEASUREMENT CASE (Aaronson-Gottesman Case II):
                // All stabilizers commute with Z_addr0 (no X or Y at addr0).
                // Express Z_addr0 as a product of stabilizer generators by
                // examining which destabilizers anticommute with Z_addr0.
                // For each destabilizer with X at addr0, multiply the
                // corresponding stabilizer into a scratch accumulator.
                let mut scratch = PhasedPauliWord::<T::Storage, T::BuildHasher>::new(N);
                for (i, destab) in self.destabilizers.iter().enumerate() {
                    if destab.word.xbits[addr0] {
                        scratch *= self.stabilizers[i].clone();
                    }
                }

                // Phase encoding: 0 → +1, 2 → -1
                // phase >= 2 means -Z eigenvalue → outcome |1⟩ (true)
                scratch.phase >= 2
            }
        }
    }
}

impl<const N: usize, T: Config, C: SparseVector<Complex<T::Coeff>>> Measure
    for GeneralizedTableau<N, T, C>
where
    T: Config,
    C: SparseVector<Complex<T::Coeff>>,
    T::Coeff: One + Zero + Clone + num::Num,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat,
    <Complex<T::Coeff> as ComplexFloat>::Real: num::ToPrimitive,
{
    /// Measure qubit `addr0` in Z basis for a generalized stabilizer state.
    ///
    /// The state is represented as |ψ⟩ = Σ_i α_i |φ_i⟩ where each |φ_i⟩ is a stabilizer state.
    /// This implements measurement by:
    /// 1. Computing the Z eigenvalue (+1 or -1) for each stabilizer state
    /// 2. Computing Born probabilities for each outcome
    /// 3. Sampling according to Born rule
    /// 4. Keeping only coefficients with matching eigenvalue
    /// 5. Renormalizing the state
    fn measure(&mut self, addr0: usize) -> bool {
        // Step 1: For each basis state (indexed by idx), compute eigenvalue of Z_addr0
        // The eigenvalue is determined by the destabilizers
        let eigenvalues: Vec<(bool, usize)> = self
            .coefficients
            .clone()
            .into_iter()
            .map(|(_, idx)| {
                // Track sign: false = +1 eigenvalue, true = -1 eigenvalue
                let mut sign = false;

                // Check each destabilizer
                for (i, destab) in self.tableau.destabilizers.iter().enumerate() {
                    let bit_i_set = (idx & (1 << i)) != 0;

                    // Case 1: Destabilizer has Z at addr0 (like ZI) → directly determines eigenvalue
                    // If bit i is set, we're in the -1 eigenstate, so Z eigenvalue is -1
                    if !destab.word.xbits[addr0] && destab.word.zbits[addr0] {
                        if bit_i_set {
                            sign = !sign; // -1 eigenstate of Z means -1 eigenvalue
                        }
                    }
                    // Case 2: Destabilizer has X or Y at addr0 → Z anticommutes with it
                    // This contributes a phase flip when the bit is set
                    else if destab.word.xbits[addr0] && bit_i_set {
                        sign = !sign;
                    }
                }

                // Return (is_minus_one_eigenvalue, index)
                (sign, idx)
            })
            .collect();

        // Step 2: Compute Born probabilities for each outcome
        let mut prob_plus_one = 0.0_f64;
        let mut prob_minus_one = 0.0_f64;

        for ((coeff, _), (is_minus_one, _)) in
            self.coefficients.clone().into_iter().zip(&eigenvalues)
        {
            // |α_i|² = |re|² + |im|²
            use num::ToPrimitive;
            let re_sq = (coeff.re() * coeff.re()).to_f64().unwrap_or(0.0);
            let im_sq = (coeff.im() * coeff.im()).to_f64().unwrap_or(0.0);
            let prob = re_sq + im_sq;
            if *is_minus_one {
                prob_minus_one += prob;
            } else {
                prob_plus_one += prob;
            }
        }

        // Step 3: Sample outcome according to Born rule
        let total_prob = prob_plus_one + prob_minus_one;
        let outcome = if total_prob == 0.0 {
            // Degenerate case: assume |0⟩ outcome
            false
        } else {
            // Generate random number in [0, 1]
            let r: f64 = rand::random::<f64>();
            // If r < P(+1) / P(total), measure +1 (outcome=false)
            // Otherwise measure -1 (outcome=true)
            let threshold = prob_plus_one / total_prob;
            r >= threshold
        };

        // Step 4: Keep only coefficients with matching eigenvalue
        self.coefficients.retain(|(_, idx)| {
            let eigenvalue_is_minus_one = eigenvalues
                .iter()
                .find(|(_, i)| i == idx)
                .map(|(sign, _)| *sign)
                .unwrap_or(false);
            eigenvalue_is_minus_one == outcome
        });

        // Step 5: Update the underlying tableau to reflect the measurement
        self.tableau.measure(addr0);

        // Step 6: Renormalize the remaining state
        if !self.coefficients.is_empty() {
            self.coefficients.normalize();
        }

        outcome
    }
}
