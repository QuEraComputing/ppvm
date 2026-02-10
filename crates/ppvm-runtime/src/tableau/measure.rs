use super::data::{GeneralizedTableau, Tableau};
use super::traits::Measure;
use crate::config::Config;
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
                for i in 0..N {
                    if i != p_idx && self.destabilizers[i].word.xbits[addr0] {
                        let g_p = self.stabilizers[p_idx].clone();
                        self.destabilizers[i] *= g_p;
                    }
                }

                // Replace stabilizer p with ±Z_addr0 based on outcome
                // outcome = false (|0⟩) → +Z, outcome = true (|1⟩) → -Z
                self.stabilizers[p_idx].word.xbits.set(addr0, false);
                self.stabilizers[p_idx].word.zbits.set(addr0, true);
                for i in 0..self.stabilizers[p_idx].n_qubits() {
                    if i != addr0 {
                        self.stabilizers[p_idx].word.xbits.set(i, false);
                        self.stabilizers[p_idx].word.zbits.set(i, false);
                    }
                }
                // Phase: 0 for +Z (outcome=false), 2 for -Z (outcome=true)
                self.stabilizers[p_idx].phase = if outcome { 2 } else { 0 };

                // Update corresponding destabilizer to X_addr0
                self.destabilizers[p_idx].word.xbits.set(addr0, true);
                self.destabilizers[p_idx].word.zbits.set(addr0, false);
                for i in 0..self.destabilizers[p_idx].n_qubits() {
                    if i != addr0 {
                        self.destabilizers[p_idx].word.xbits.set(i, false);
                        self.destabilizers[p_idx].word.zbits.set(i, false);
                    }
                }
                self.destabilizers[p_idx].phase = 0;

                outcome
            }
            None => {
                // DETERMINISTIC MEASUREMENT CASE:
                // All stabilizers commute with Z_addr0 (no X or Y at addr0)
                // The outcome is determined by the product of stabilizers with Z at addr0

                let mut phase = 0u8;
                for stab in self.stabilizers.iter() {
                    if stab.word.zbits[addr0] {
                        // This stabilizer has Z (or Y, but Y would have X bit set)
                        phase = (phase + stab.phase) % 4;
                    }
                }

                // Phase encoding: 0 → +1, 1 → +i, 2 → -1, 3 → -i
                // For Z eigenstates: phase 0,1 → |0⟩ (false), phase 2,3 → |1⟩ (true)
                // This is because -Z|1⟩ = -|1⟩, so negative phase means |1⟩ state
                phase >= 2
            }
        }
    }
}
