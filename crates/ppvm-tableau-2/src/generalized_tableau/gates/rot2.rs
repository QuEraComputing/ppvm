// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Pauli, RotationTwo, RotationTwoBatch};

use super::super::GeneralizedTableau;
use crate::VirtualState;

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

impl<V: VirtualState> RotationTwoBatch<f64> for GeneralizedTableau<V> {}
