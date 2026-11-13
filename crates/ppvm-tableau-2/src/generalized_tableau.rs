// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`GeneralizedTableau`]: the simulator, `|ψ⟩ = C |χ⟩`.

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use num_complex::Complex64;
use ppvm_traits_2::{
    ApplyGateRule, Clifford, CliffordBatch, Pauli, RotationOne, RotationTwo, TGate,
};

use crate::{InverseTableau, VirtualState};

/// The state `|ψ⟩ = C |χ⟩`: the Clifford `C` as an [`InverseTableau`], and the
/// virtual state `|χ⟩` in the representation `V`.
///
/// Physical gates `C ← G·C` touch only the tableau. Measurements ask `V` for
/// the outcome weights, sample here, and let `V` collapse.
#[derive(Clone, Debug)]
pub struct GeneralizedTableau<V: VirtualState> {
    tableau: InverseTableau,
    state: V,
    rng: SmallRng,
}

impl<V: VirtualState> GeneralizedTableau<V> {
    /// `|0…0⟩` on `n` qubits, with the RNG seeded from OS entropy.
    pub fn new(n: usize) -> Self {
        Self::with_rng(n, rand::make_rng())
    }

    /// Same as [`new`](Self::new), but seed the RNG deterministically.
    pub fn new_with_seed(n: usize, seed: u64) -> Self {
        Self::with_rng(n, SmallRng::seed_from_u64(seed))
    }

    fn with_rng(n: usize, rng: SmallRng) -> Self {
        Self {
            tableau: InverseTableau::new(n),
            state: V::new(n),
            rng,
        }
    }

    pub fn num_qubits(&self) -> usize {
        self.tableau.num_qubits()
    }

    pub fn tableau(&self) -> &InverseTableau {
        &self.tableau
    }

    pub fn state(&self) -> &V {
        &self.state
    }

    /// Measure `Z_q` and collapse; returns the outcome (`true` for `|1⟩`).
    pub fn measure(&mut self, q: usize) -> bool {
        let mut view = self.tableau.virtual_view();
        let [w0, w1] = self.state.z_weights(&mut view, q);
        let m = self.rng.random::<f64>() * (w0 + w1) < w1;
        self.state.collapse(&mut view, q, m);
        m
    }

    /// Reset qubit `q` to `|0⟩`: measure it, then flip it if it was `|1⟩`.
    pub fn reset(&mut self, q: usize) {
        if self.measure(q) {
            self.x(q);
        }
    }

    /// `|ψ⟩ ← (a·I + b·∏ Z_q) |ψ⟩` over `qubits`, which must be distinct.
    fn rotate_z(&mut self, qubits: &[usize], a: Complex64, b: Complex64) {
        let mut view = self.tableau.virtual_view();
        self.state.rotate(&mut view, qubits, a, b);
    }

    /// `exp(−i θ/2 ∏ Z_q) = cos(θ/2)·I − i sin(θ/2)·∏ Z_q`.
    fn rz_product(&mut self, qubits: &[usize], theta: f64) {
        let (s, c) = (theta / 2.0).sin_cos();
        self.rotate_z(qubits, Complex64::new(c, 0.0), Complex64::new(0.0, -s));
    }

    /// Map `axis` on `q` to `Z`: apply `W` with `W·axis·W† = Z`.
    fn axis_to_z(&mut self, axis: Pauli, q: usize) {
        match axis {
            Pauli::X => self.h(q),
            Pauli::Y => {
                self.s_dag(q);
                self.h(q);
            }
            Pauli::Z | Pauli::I => {}
        }
    }

    /// Undo [`axis_to_z`](Self::axis_to_z): apply `W†`.
    fn z_to_axis(&mut self, axis: Pauli, q: usize) {
        match axis {
            Pauli::X => self.h(q),
            Pauli::Y => {
                self.h(q);
                self.s(q);
            }
            Pauli::Z | Pauli::I => {}
        }
    }
}

/// `T = diag(1, ω) = ½(1 + ω)·I + ½(1 − ω)·Z` with `ω = e^{iπ/4}`, exact
/// including the global phase.
impl<V: VirtualState> TGate for GeneralizedTableau<V> {
    fn t(&mut self, q: usize) {
        let w = Complex64::from_polar(1.0, std::f64::consts::FRAC_PI_4);
        self.rotate_z(&[q], (1.0 + w) / 2.0, (1.0 - w) / 2.0);
    }

    fn t_dag(&mut self, q: usize) {
        let w = Complex64::from_polar(1.0, -std::f64::consts::FRAC_PI_4);
        self.rotate_z(&[q], (1.0 + w) / 2.0, (1.0 - w) / 2.0);
    }
}

/// `exp(−i θ/2 P)`: `P` is mapped to `Z` by a physical Clifford.
impl<V: VirtualState> RotationOne<f64> for GeneralizedTableau<V> {
    fn rotate_1(&mut self, axis: Pauli, q: usize, theta: f64) {
        if axis == Pauli::I {
            return;
        }
        self.axis_to_z(axis, q);
        self.rz_product(&[q], theta);
        self.z_to_axis(axis, q);
    }
}

/// `exp(−i θ/2 P_a ⊗ P_b)`: each factor is mapped to `Z` by a physical Clifford.
impl<V: VirtualState> RotationTwo<f64> for GeneralizedTableau<V> {
    fn rotate_2(&mut self, axis_a: Pauli, axis_b: Pauli, a: usize, b: usize, theta: f64) {
        assert_ne!(a, b, "two-qubit rotation on one qubit");
        let qubits: Vec<usize> = [(axis_a, a), (axis_b, b)]
            .into_iter()
            .filter(|&(axis, _)| axis != Pauli::I)
            .map(|(_, q)| q)
            .collect();
        self.axis_to_z(axis_a, a);
        self.axis_to_z(axis_b, b);
        self.rz_product(&qubits, theta);
        self.z_to_axis(axis_a, a);
        self.z_to_axis(axis_b, b);
    }
}

/// Physical gates `C ← G·C` go to the tableau; `|χ⟩` is untouched.
impl<V: VirtualState> ApplyGateRule for GeneralizedTableau<V> {
    type B = bool;

    fn apply1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.tableau.apply1(q, rule);
    }

    fn apply2(&mut self, a: usize, b: usize, rule: impl Fn([bool; 4]) -> ([bool; 4], bool)) {
        self.tableau.apply2(a, b, rule);
    }

    fn apply_sign1(&mut self, q: usize, rule: impl Fn(bool, bool) -> ([bool; 2], bool)) {
        self.tableau.apply_sign1(q, rule);
    }
}

impl<V: VirtualState> CliffordBatch for GeneralizedTableau<V> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ActiveDormant, DenseState};

    type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

    /// Bell pair: the two outcomes always agree, and both occur.
    #[test]
    fn bell_pair_outcomes_agree() {
        let mut ones = 0;
        for seed in 0..200 {
            let mut s = Sim::new_with_seed(2, seed);
            s.h(0);
            s.cnot(0, 1);
            let (a, b) = (s.measure(0), s.measure(1));
            assert_eq!(a, b);
            ones += usize::from(a);
        }
        assert!((60..140).contains(&ones), "{ones} of 200");
    }

    /// A measured qubit stays in its outcome, and reset leaves `|0⟩`.
    #[test]
    fn measurement_repeats_and_reset_clears() {
        for seed in 0..50 {
            let mut s = Sim::new_with_seed(3, seed);
            s.h(0);
            s.cnot(0, 1);
            s.sqrt_x(2);
            s.cz(1, 2);
            for q in 0..3 {
                let m = s.measure(q);
                assert_eq!(s.measure(q), m);
            }
            for q in 0..3 {
                s.reset(q);
                assert!(!s.measure(q));
            }
        }
    }

    /// `X` then measure is deterministic: `Z|1⟩ = −|1⟩`.
    #[test]
    fn x_flips_the_outcome() {
        let mut s = Sim::new_with_seed(1, 0);
        s.x(0);
        assert!(s.measure(0));
    }
}
