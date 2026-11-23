// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Heisenberg amplitude damping: X, Y scale by sqrt(1 - gamma); Z maps to
//! (1 - gamma) Z + gamma I through the branching machinery. Lost sites are
//! skipped. Zero branches stay stored; reduction is caller-controlled.

use num::One;
use ppvm_traits_2::{AmplitudeDamping, PauliBits, loss::LossInfo};

use crate::{BranchInPlace, Sum};

impl<S: BranchInPlace, P> AmplitudeDamping<S::Coeff> for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: num::Float,
{
    fn amplitude_damping(&mut self, qubit: usize, gamma: S::Coeff) {
        self.check_qubit(qubit);
        let survive = S::Coeff::one() - gamma;
        let transverse = num::Float::sqrt(survive);
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(qubit) {
                return None;
            }
            if key.x_bit(qubit) {
                *coeff *= transverse;
                return None;
            }
            if !key.z_bit(qubit) {
                return None;
            }
            let branch = *coeff * gamma;
            *coeff *= survive;
            Some((key.toggled_bits(qubit, false, true), branch))
        });
    }
}
