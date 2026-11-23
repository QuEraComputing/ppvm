// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Z-basis projectors as `O ↦ P_m O P_m`, the same map in either picture.
//! On the projected site `I` and `Z` both become `(I ± Z)/2`; `X`, `Y`, and lost
//! sites become zero, which stays stored. Results are unnormalized, and
//! reduction and truncation stay caller-controlled.

use num::Zero;
use ppvm_traits_2::{Halvable, PauliBits, Projection, loss::LossInfo};

use crate::{BranchInPlace, Sum};

impl<S: BranchInPlace, P> Projection for Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: Halvable,
{
    fn p0(&mut self, qubit: usize) {
        self.project(qubit, false);
    }

    fn p1(&mut self, qubit: usize) {
        self.project(qubit, true);
    }
}

impl<S: BranchInPlace, P> Sum<S, P>
where
    S::Key: PauliBits + LossInfo,
    S::Coeff: Halvable,
{
    fn project(&mut self, qubit: usize, one: bool) {
        self.check_qubit(qubit);
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(qubit) || key.x_bit(qubit) {
                *coeff = S::Coeff::zero();
                return None;
            }
            // Keep half and send ±half to the I <-> Z partner: P0 maps I and Z to
            // (I + Z)/2; P1 maps I to (I - Z)/2 and Z to (Z - I)/2.
            *coeff = coeff.half();
            let partner = if one { -coeff.clone() } else { coeff.clone() };
            Some((key.toggled_bits(qubit, false, true), partner))
        });
    }
}
