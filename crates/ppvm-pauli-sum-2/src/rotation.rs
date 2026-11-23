// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Heisenberg rotations, with reduction and truncation left to the caller.
//! Two-site generators drop lost sites, leaving a rotation on the survivor.
//! CRx inherits this convention through Rx(theta/2) and RZX(-theta/2): a lost
//! control cancels the rotations; a lost target leaves Rz(-theta/2) on control.

use ppvm_traits_2::{
    Angle, CRx, Coefficient, Pauli, PauliBits, RotXY, RotationOne, RotationOneBatch, RotationTwo,
    RotationTwoBatch, Support, TGate, U3Gate, loss::LossInfo,
};

use crate::{BranchInPlace, Sum};

impl<S, P, A> RotationOne<S::Coeff, A> for Sum<S, P>
where
    S: BranchInPlace,
    S::Key: PauliBits + LossInfo,
    A: Angle<S::Coeff>,
{
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: A) {
        assert!(qubit < self.n_sites(), "qubit out of bounds");
        if axis == Pauli::I {
            return;
        }
        let (sin, cos) = theta.sin_cos();
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(qubit) {
                return None;
            }
            let x = key.x_bit(qubit);
            let z = key.z_bit(qubit);
            // Heisenberg conjugation: Rx sends Z to cos(theta) Z + sin(theta) Y.
            let (toggle_x, toggle_z, negative) = match axis {
                Pauli::X if z => (true, false, x),
                Pauli::Y if x != z => (true, true, z),
                Pauli::Z if x => (false, true, !z),
                _ => return None,
            };
            let branch = coeff.clone() * sin.mul_sign(if negative { -1 } else { 1 });
            *coeff *= &cos;
            Some((key.toggled_bits(qubit, toggle_x, toggle_z), branch))
        });
    }
}

impl<S, P, A> RotationOneBatch<S::Coeff, A> for Sum<S, P>
where
    S: BranchInPlace,
    S::Key: PauliBits + LossInfo,
    A: Angle<S::Coeff>,
{
}

impl<S, P, A> RotationTwo<S::Coeff, A> for Sum<S, P>
where
    S: BranchInPlace,
    S::Key: PauliBits + LossInfo,
    A: Angle<S::Coeff>,
{
    fn rotate_2(&mut self, axis_a: Pauli, axis_b: Pauli, a: usize, b: usize, theta: A) {
        assert!(
            a < self.n_sites() && b < self.n_sites(),
            "qubit out of bounds"
        );
        assert_ne!(a, b, "rotation requires distinct qubits");
        let (sin, cos) = theta.sin_cos();
        self.storage.branch_in_place(|key, coeff| {
            // A lost site drops out of the generator, preserving the reference's
            // single-qubit fallback on the surviving site.
            let axis_a = if key.is_lost(a) { Pauli::I } else { axis_a };
            let axis_b = if key.is_lost(b) { Pauli::I } else { axis_b };
            let (mask_a, anti_a, negative_a) = local_product(axis_a, key.x_bit(a), key.z_bit(a));
            let (mask_b, anti_b, negative_b) = local_product(axis_b, key.x_bit(b), key.z_bit(b));
            if anti_a == anti_b {
                return None;
            }
            // Exactly one local product anticommutes; the other has phase +1.
            let negative = if anti_a { negative_a } else { negative_b };
            let branch = coeff.clone() * sin.mul_sign(if negative { -1 } else { 1 });
            *coeff *= &cos;
            Some((key.toggled_bits2(a, mask_a, b, mask_b), branch))
        });
    }

    #[inline]
    fn rxx(&mut self, a: usize, b: usize, theta: A) {
        assert!(
            a < self.n_sites() && b < self.n_sites(),
            "qubit out of bounds"
        );
        assert_ne!(a, b, "rotation requires distinct qubits");
        let (sin, cos) = theta.sin_cos();
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(a) {
                return single_branch(key, coeff, Pauli::X, b, &sin, &cos);
            }
            if key.is_lost(b) {
                return single_branch(key, coeff, Pauli::X, a, &sin, &cos);
            }
            let za = key.z_bit(a);
            let zb = key.z_bit(b);
            if za == zb {
                return None;
            }
            let xa = key.x_bit(a);
            let xb = key.x_bit(b);
            let x_anti = if za { xa } else { xb };
            let sign = if x_anti { -1 } else { 1 };
            let branch = coeff.clone() * sin.mul_sign(sign);
            *coeff *= &cos;
            Some((
                key.toggled_bits2(a, [true, false], b, [true, false]),
                branch,
            ))
        });
    }

    #[inline]
    fn ryy(&mut self, a: usize, b: usize, theta: A) {
        assert!(
            a < self.n_sites() && b < self.n_sites(),
            "qubit out of bounds"
        );
        assert_ne!(a, b, "rotation requires distinct qubits");
        let (sin, cos) = theta.sin_cos();
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(a) {
                return single_branch(key, coeff, Pauli::Y, b, &sin, &cos);
            }
            if key.is_lost(b) {
                return single_branch(key, coeff, Pauli::Y, a, &sin, &cos);
            }
            let xa = key.x_bit(a);
            let za = key.z_bit(a);
            let xb = key.x_bit(b);
            let zb = key.z_bit(b);
            let anti_a = xa != za;
            if anti_a == (xb != zb) {
                return None;
            }
            let x_anti = if anti_a { xa } else { xb };
            let sign = if x_anti { 1 } else { -1 };
            let branch = coeff.clone() * sin.mul_sign(sign);
            *coeff *= &cos;
            Some((key.toggled_bits2(a, [true, true], b, [true, true]), branch))
        });
    }

    #[inline]
    fn rzz(&mut self, a: usize, b: usize, theta: A) {
        assert!(
            a < self.n_sites() && b < self.n_sites(),
            "qubit out of bounds"
        );
        assert_ne!(a, b, "rotation requires distinct qubits");
        let (sin, cos) = theta.sin_cos();
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(a) {
                return single_branch(key, coeff, Pauli::Z, b, &sin, &cos);
            }
            if key.is_lost(b) {
                return single_branch(key, coeff, Pauli::Z, a, &sin, &cos);
            }
            let xa = key.x_bit(a);
            let xb = key.x_bit(b);
            if xa == xb {
                return None;
            }
            let za = key.z_bit(a);
            let zb = key.z_bit(b);
            let z_anti = if xa { za } else { zb };
            let sign = if z_anti { 1 } else { -1 };
            let branch = coeff.clone() * sin.mul_sign(sign);
            *coeff *= &cos;
            Some((
                key.toggled_bits2(a, [false, true], b, [false, true]),
                branch,
            ))
        });
    }
}

/// Generator mask, anticommutation, and sign of i * axis * word (when anticommuting).
#[inline]
fn local_product(axis: Pauli, x: bool, z: bool) -> ([bool; 2], bool, bool) {
    match axis {
        Pauli::I => ([false, false], false, false),
        Pauli::X => ([true, false], z, x),
        Pauli::Y => ([true, true], x != z, z),
        Pauli::Z => ([false, true], x, !z),
    }
}

impl<S, P, A> RotationTwoBatch<S::Coeff, A> for Sum<S, P>
where
    S: BranchInPlace,
    S::Key: PauliBits + LossInfo,
    A: Angle<S::Coeff>,
{
}

impl<S: Support, P> TGate for Sum<S, P>
where
    Self: RotationOne<S::Coeff, f64>,
    f64: Angle<S::Coeff>,
{
    fn t(&mut self, qubit: usize) {
        self.rz(qubit, std::f64::consts::FRAC_PI_4);
    }

    fn t_dag(&mut self, qubit: usize) {
        self.rz(qubit, -std::f64::consts::FRAC_PI_4);
    }
}

impl<S: Support, P, A> U3Gate<S::Coeff, A> for Sum<S, P>
where
    Self: RotationOne<S::Coeff, A>,
    A: Angle<S::Coeff>,
{
    fn u3(&mut self, qubit: usize, theta: A, phi: A, lambda: A) {
        // U3 = Rz(phi) Ry(theta) Rz(lambda); propagate in Heisenberg order.
        self.rz(qubit, phi);
        self.ry(qubit, theta);
        self.rz(qubit, lambda);
    }
}

impl<S: Support, P, A> RotXY<S::Coeff, A> for Sum<S, P>
where
    Self: RotationOne<S::Coeff, A>,
    A: Angle<S::Coeff> + Clone + std::ops::Neg<Output = A>,
{
    fn r(&mut self, qubit: usize, axis_angle: A, theta: A) {
        self.rz(qubit, axis_angle.clone());
        self.rx(qubit, theta);
        self.rz(qubit, -axis_angle);
    }
}

impl<S: Support, P, A> CRx<S::Coeff, A> for Sum<S, P>
where
    Self: RotationOne<S::Coeff, A> + RotationTwo<S::Coeff, A>,
    A: Angle<S::Coeff> + Clone + std::ops::Mul<f64, Output = A> + std::ops::Neg<Output = A>,
{
    fn crx(&mut self, control: usize, target: usize, theta: A) {
        assert!(
            control < self.n_sites() && target < self.n_sites(),
            "qubit out of bounds"
        );
        assert_ne!(control, target, "rotation requires distinct qubits");
        // |1><1| ⊗ X = (IX - ZX)/2. The two generators commute.
        // Loss follows the constituent rotations' surviving-site convention.
        let half = theta * 0.5;
        self.rx(target, half.clone());
        self.rzx(control, target, -half);
    }
}

// Single-site fallback for a two-site rotation whose other site is lost.
#[inline]
fn single_branch<W: PauliBits + LossInfo + Clone, C: Coefficient>(
    key: &W,
    coeff: &mut C,
    axis: Pauli,
    qubit: usize,
    sin: &C,
    cos: &C,
) -> Option<(W, C)> {
    if key.is_lost(qubit) {
        return None;
    }
    let (mask, anti, negative) = local_product(axis, key.x_bit(qubit), key.z_bit(qubit));
    if !anti {
        return None;
    }
    let branch = coeff.clone() * sin.mul_sign(if negative { -1 } else { 1 });
    *coeff *= cos;
    Some((key.toggled_bits(qubit, mask[0], mask[1]), branch))
}
