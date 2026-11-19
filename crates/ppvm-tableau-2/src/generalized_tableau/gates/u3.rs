// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{RotationOne, U3Gate};

use super::super::GeneralizedTableau;
use crate::VirtualState;

/// `U3(θ, φ, λ) = RZ(φ)·RY(θ)·RZ(λ)` up to the global phase `e^{i(φ+λ)/2}`.
/// Gates apply in circuit order, `RZ(λ)` first.
impl<V: VirtualState> U3Gate<f64> for GeneralizedTableau<V> {
    fn u3(&mut self, q: usize, theta: f64, phi: f64, lambda: f64) {
        self.rz(q, lambda);
        self.ry(q, theta);
        self.rz(q, phi);
    }
}
