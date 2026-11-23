// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Heisenberg loss maps matching the existing engine. Loss branches clear an
//! observable's loss flag; reset maps I/Z to an additional lost observable.
//! Zero coefficients stay stored, and truncation remains caller-controlled.

use crate::{BranchInPlace, Sum};
use num::{One, Zero};
use ppvm_traits_2::{
    CorrelatedLossChannel, LossChannel, PauliBits, ResetLossChannel,
    loss::{LossInfo, LossState},
};

impl<S: BranchInPlace, P> LossChannel<S::Coeff> for Sum<S, P>
where
    S::Key: LossState,
{
    fn loss_channel(&mut self, qubit: usize, p: S::Coeff) {
        self.check_qubit(qubit);
        let survive = S::Coeff::one() - p.clone();
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(qubit) {
                Some((key.loss_cleared(qubit), coeff.clone() * p.clone()))
            } else {
                *coeff *= &survive;
                None
            }
        });
    }
}

impl<S: BranchInPlace, P> CorrelatedLossChannel<S::Coeff> for Sum<S, P>
where
    S::Key: LossState,
{
    fn correlated_loss_channel(
        &mut self,
        a: usize,
        b: usize,
        [both, single, after_loss]: [S::Coeff; 3],
    ) {
        self.check_pair(a, b);
        let survive = S::Coeff::one() - (single.clone() + single.clone()) - both.clone();
        let one_lost = S::Coeff::one() - after_loss.clone();
        self.storage.branch_in_place(|key, coeff| {
            let mut branches = [None, None, None];
            match (key.is_lost(a), key.is_lost(b)) {
                (true, true) => {
                    let weight = coeff.clone() * after_loss.clone();
                    branches[0] = Some((key.loss_cleared(a), weight.clone()));
                    branches[1] = Some((key.loss_cleared(b), weight));
                    let mut cleared = key.loss_cleared(a);
                    cleared.clear_lost(b);
                    branches[2] = Some((cleared, coeff.clone() * both.clone()));
                }
                (lost_a, lost_b) if lost_a || lost_b => {
                    let qubit = if lost_a { a } else { b };
                    branches[0] = Some((key.loss_cleared(qubit), coeff.clone() * single.clone()));
                    *coeff *= &one_lost;
                }
                _ => *coeff *= &survive,
            }
            branches.into_iter().flatten()
        });
    }
}

impl<S: BranchInPlace, P> ResetLossChannel for Sum<S, P>
where
    S::Key: PauliBits + LossState,
{
    fn reset_loss_channel(&mut self, qubit: usize) {
        self.check_qubit(qubit);
        self.storage.branch_in_place(|key, coeff| {
            if key.is_lost(qubit) {
                *coeff = S::Coeff::zero();
                return None;
            }
            if key.x_bit(qubit) {
                return None;
            }
            let mut lost = key.clone();
            lost.set_lost(qubit);
            Some((lost, coeff.clone()))
        });
    }
}
