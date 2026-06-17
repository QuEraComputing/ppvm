// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::char::Pauli;
use crate::config::Config;

/// Single-qubit Pauli rotations `exp(-i θ/2 · P)`.
pub trait RotationOne<T: Config> {
    /// Rotate about `axis` (one of `X`, `Y`, `Z`) by angle `theta`.
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: T::Coeff);
    /// `RX(θ)` on each target qubit.
    fn rx(&mut self, targets: impl crate::traits::Targets, theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for q in targets.each() {
            self.rotate_1(Pauli::X, q, theta.clone())
        }
    }
    /// `RY(θ)` on each target qubit.
    fn ry(&mut self, targets: impl crate::traits::Targets, theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for q in targets.each() {
            self.rotate_1(Pauli::Y, q, theta.clone())
        }
    }
    /// `RZ(θ)` on each target qubit.
    fn rz(&mut self, targets: impl crate::traits::Targets, theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for q in targets.each() {
            self.rotate_1(Pauli::Z, q, theta.clone())
        }
    }
}
