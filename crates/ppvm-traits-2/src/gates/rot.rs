// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::{
    arithmetic::{Angle, Coefficient},
    pauli::Pauli,
};

/// Single-qubit rotations parameterized by an angle domain `A` that yields
/// amplitudes in coefficient domain `C`.
///
/// The angle defaults to the coefficient (`A = C`), recovering today's
/// `rx(theta: C)` while permitting a symbolic/parametric angle over an
/// `f64`-coefficient sum.
pub trait RotationOne<C: Coefficient, A: Angle<C> = C> {
    /// Rotate about `axis` (one of `X`, `Y`, `Z`) on `qubit` by `theta`.
    ///
    /// `Pauli::I` commutes with every term, so an `I` axis is a no-op.
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: A);

    /// Rotate about `X` on `qubit` by `theta`.
    fn rx(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::X, qubit, theta)
    }
    /// Rotate about `Y` on `qubit` by `theta`.
    fn ry(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::Y, qubit, theta)
    }
    /// Rotate about `Z` on `qubit` by `theta`.
    fn rz(&mut self, qubit: usize, theta: A) {
        self.rotate_1(Pauli::Z, qubit, theta)
    }

    /// Explicit batched `RX(θ)`.
    fn rx_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.rx(q, theta.clone())
        }
    }
    /// Explicit batched `RY(θ)`.
    fn ry_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.ry(q, theta.clone())
        }
    }
    /// Explicit batched `RZ(θ)`.
    fn rz_many(&mut self, targets: &[usize], theta: A)
    where
        A: Clone,
    {
        for &q in targets {
            self.rz(q, theta.clone())
        }
    }
}

/// Two-qubit rotations `exp(-i θ/2 · P_a ⊗ P_b)`.
pub trait RotationTwo<C: Coefficient, A: Angle<C> = C> {
    /// Rotate about the supplied Pauli axes.
    ///
    /// The generator is `axis_a ⊗ axis_b` on sites `a` and `b`.
    fn rotate_2(&mut self, axis_a: Pauli, axis_b: Pauli, a: usize, b: usize, theta: A);

    /// Rotate about X ⊗ X.
    fn rxx(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::X, Pauli::X, a, b, theta);
    }

    /// Rotate about X ⊗ Y.
    fn rxy(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::X, Pauli::Y, a, b, theta);
    }

    /// Rotate about X ⊗ Z.
    fn rxz(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::X, Pauli::Z, a, b, theta);
    }

    /// Rotate about Y ⊗ X.
    fn ryx(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Y, Pauli::X, a, b, theta);
    }

    /// Rotate about Y ⊗ Y.
    fn ryy(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Y, Pauli::Y, a, b, theta);
    }

    /// Rotate about Y ⊗ Z.
    fn ryz(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Y, Pauli::Z, a, b, theta);
    }

    /// Rotate about Z ⊗ X.
    fn rzx(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Z, Pauli::X, a, b, theta);
    }

    /// Rotate about Z ⊗ Y.
    fn rzy(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Z, Pauli::Y, a, b, theta);
    }

    /// Rotate about Z ⊗ Z.
    fn rzz(&mut self, a: usize, b: usize, theta: A) {
        self.rotate_2(Pauli::Z, Pauli::Z, a, b, theta);
    }

    /// Apply RXX to each pair in order.
    fn rxx_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rxx(a, b, theta.clone());
        }
    }

    /// Apply RXY to each pair in order.
    fn rxy_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rxy(a, b, theta.clone());
        }
    }

    /// Apply RXZ to each pair in order.
    fn rxz_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rxz(a, b, theta.clone());
        }
    }

    /// Apply RYX to each pair in order.
    fn ryx_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.ryx(a, b, theta.clone());
        }
    }

    /// Apply RYY to each pair in order.
    fn ryy_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.ryy(a, b, theta.clone());
        }
    }

    /// Apply RYZ to each pair in order.
    fn ryz_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.ryz(a, b, theta.clone());
        }
    }

    /// Apply RZX to each pair in order.
    fn rzx_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rzx(a, b, theta.clone());
        }
    }

    /// Apply RZY to each pair in order.
    fn rzy_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rzy(a, b, theta.clone());
        }
    }

    /// Apply RZZ to each pair in order.
    fn rzz_many(&mut self, pairs: &[(usize, usize)], theta: A)
    where
        A: Clone,
    {
        for &(a, b) in pairs {
            self.rzz(a, b, theta.clone());
        }
    }
}

/// Rotation about an axis in the x/y plane:
/// `R(axis_angle, θ) = exp(−i·θ/2·(cos(axis_angle)·X + sin(axis_angle)·Y))`.
///
/// The in-plane axis is `X` rotated about `Z` by `axis_angle`, so
/// `R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle)`
pub trait RotXY<C: Coefficient, A: Angle<C> = C> {
    /// `R(axis_angle, θ)` on `qubit`.
    fn r(&mut self, qubit: usize, axis_angle: A, theta: A);
}

/// Controlled `RX` rotation
pub trait CRx<C: Coefficient, A: Angle<C> = C> {
    /// Apply `CRX(θ)` with the given control and target.
    fn crx(&mut self, control: usize, target: usize, theta: A);
}

/// The general single-qubit `U3(θ, φ, λ)` gate
pub trait U3Gate<C: Coefficient, A: Angle<C> = C> {
    /// Apply `U3(θ, φ, λ)` to `qubit`.
    fn u3(&mut self, qubit: usize, theta: A, phi: A, lambda: A);
}

/// The non-Clifford `T` gate and its adjoint, `T = diag(1, e^{iπ/4})`.
pub trait TGate {
    /// Apply `T` (`diag(1, e^{iπ/4})`) to one qubit.
    fn t(&mut self, qubit: usize);
    /// Apply `T†` to one qubit.
    fn t_dag(&mut self, qubit: usize);

    /// Explicit batched `T`.
    fn t_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t(q);
        }
    }

    /// Explicit batched `T†`.
    fn t_dag_many(&mut self, targets: &[usize]) {
        for &q in targets {
            self.t_dag(q);
        }
    }
}
