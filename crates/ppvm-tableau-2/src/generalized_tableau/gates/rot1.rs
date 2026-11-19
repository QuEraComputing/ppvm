// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Pauli, RotXY, RotationOne, RotationOneBatch};

use super::super::GeneralizedTableau;
use crate::VirtualState;

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

impl<V: VirtualState> RotationOneBatch<f64> for GeneralizedTableau<V> {}

/// `R(φ, θ) = RZ(φ)·RX(θ)·RZ(−φ)`, exact: `RZ(φ)` turns `X` into
/// `cos φ·X + sin φ·Y`. Gates apply in circuit order, `RZ(−φ)` first.
impl<V: VirtualState> RotXY<f64> for GeneralizedTableau<V> {
    fn r(&mut self, q: usize, axis_angle: f64, theta: f64) {
        self.rz(q, -axis_angle);
        self.rx(q, theta);
        self.rz(q, axis_angle);
    }
}
