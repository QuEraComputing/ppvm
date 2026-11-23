// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_traits_2::{Angle, Coefficient, Pauli, PauliBits, RotationOne, loss::LossInfo};

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
