// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::char::Pauli;
use crate::config::Config;

/// Single-qubit Pauli rotations `exp(-i θ/2 · P)`.
pub trait RotationOne<T: Config> {
    /// Rotate about `axis` (one of `X`, `Y`, `Z`) by angle `theta`.
    fn rotate_1(&mut self, axis: Pauli, addr0: usize, theta: T::Coeff);
    /// `RX(θ)` on one qubit.
    fn rx(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::X, addr0, theta.into())
    }
    /// `RY(θ)` on one qubit.
    fn ry(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Y, addr0, theta.into())
    }
    /// `RZ(θ)` on one qubit.
    fn rz(&mut self, addr0: usize, theta: impl Into<T::Coeff>) {
        self.rotate_1(Pauli::Z, addr0, theta.into())
    }

    /// Explicit batched `RX(θ)`.
    fn rx_many(&mut self, targets: &[usize], theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for &q in targets {
            self.rx(q, theta.clone())
        }
    }
    /// Explicit batched `RY(θ)`.
    fn ry_many(&mut self, targets: &[usize], theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for &q in targets {
            self.ry(q, theta.clone())
        }
    }
    /// Explicit batched `RZ(θ)`.
    fn rz_many(&mut self, targets: &[usize], theta: impl Into<T::Coeff>) {
        let theta = theta.into();
        for &q in targets {
            self.rz(q, theta.clone())
        }
    }
}

/// Rotation about an axis in the x/y plane:
/// `R(axis_angle, θ) = exp(-i θ/2 · (cos(axis_angle)·X + sin(axis_angle)·Y))`.
///
/// The in-plane axis is `X` rotated about `Z` by `axis_angle`, so
/// `R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle)`.
pub trait RotXY<T: Config> {
    /// `R(axis_angle, θ)` on qubit `addr0`.
    fn r(&mut self, addr0: usize, axis_angle: T::Coeff, theta: T::Coeff);
}
