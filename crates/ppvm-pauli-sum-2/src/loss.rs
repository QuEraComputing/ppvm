// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Heisenberg loss maps matching `ppvm-tableau-2`. `I` projects onto a live
//! site and `L` onto a lost one. Loss branches clear an observable's loss flag;
//! reset adds an `L` copy of each `I`/`Z` term, since a lost qubit returns in
//! `|0>`. Zero coefficients stay stored; reduction is caller-controlled.

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
        let partner_survives = S::Coeff::one() - after_loss.clone();
        self.storage.branch_in_place(|key, coeff| {
            let mut branches = [None, None, None];
            match (key.is_lost(a), key.is_lost(b)) {
                (false, false) => *coeff *= &survive,
                // Reached from one already lost (p_LN) or from both live (p_LL).
                (true, true) => {
                    let weight = coeff.clone() * after_loss.clone();
                    let mut neither = key.loss_cleared(a);
                    neither.clear_lost(b);
                    branches = [
                        Some((key.loss_cleared(a), weight.clone())),
                        Some((key.loss_cleared(b), weight)),
                        Some((neither, coeff.clone() * both.clone())),
                    ];
                }
                // Reached from both live (p_LQ), or the partner survived p_LN.
                (lost_a, _) => {
                    let lost = if lost_a { a } else { b };
                    branches[0] = Some((key.loss_cleared(lost), coeff.clone() * single.clone()));
                    *coeff *= &partner_survives;
                }
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
