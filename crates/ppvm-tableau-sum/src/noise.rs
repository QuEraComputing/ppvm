// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use ppvm_traits::config::Config;
use ppvm_traits::traits::{
    Clifford, CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel, PauliError,
    ResetLossChannel, TwoQubitPauliError,
};
use rand::{RngExt, rngs::SmallRng};

use crate::{
    data::GeneralizedTableauSum,
    storage::{Branch, EntryStore, loss_mask, pauli_branch_phase_loss},
};

fn single_qubit_loss_branch<T, I, C>(
    addr0: usize,
    p: &T::Coeff,
    rng: &mut SmallRng,
    branches: &mut Vec<Branch<T, I, C>>,
    tab: &mut GeneralizedTableau<T, I, C>,
    p_sum: &mut T::Coeff,
    // The branch inherits its parent's cached fingerprint halves:
    // `(word_fingerprint, phase_loss_hash)`.
    (word_fp, phase_loss): (u64, u64),
) where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Send + Sync + Debug,
    C: SparseVector<Complex<T::Coeff>, I>,
{
    if tab.is_lost[addr0] {
        // Don't branch if it's already lost
        return;
    }

    let tab_seed = rng.random::<u64>();
    let mut tab_branch = tab.fork(Some(tab_seed));
    tab_branch.is_lost[addr0] = true;
    // is_lost flip leaves the Pauli words and phases unchanged, so
    // the branch reuses its parent's word-fingerprint and the only
    // change to the phase/loss hash is the lost qubit's mask.
    branches.push((
        tab_branch,
        p_sum.clone() * p.clone(),
        word_fp,
        phase_loss ^ loss_mask(addr0),
    ));
    *p_sum *= T::Coeff::one() - p.clone();
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> LossChannel<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn loss_channel(&mut self, addr0: usize, p: <T as Config>::Coeff) {
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            self.entries.len(),
        );
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                single_qubit_loss_branch(
                    addr0,
                    &p,
                    &mut self.rng,
                    &mut branches,
                    tab,
                    p_sum,
                    (word_fp, phase_loss),
                );
            });

        let needs_renormalize = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
        if needs_renormalize {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> Depolarizing<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn depolarize(&mut self, addr0: usize, p: T::Coeff) {
        let p_3 = p.clone() / 3.0.into();
        self.pauli_error(addr0, [p_3.clone(), p_3.clone(), p_3]);
    }
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> PauliError<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn pauli_error(&mut self, addr0: usize, p: [<T as Config>::Coeff; 3]) {
        let p_total: T::Coeff = p[0].clone() + p[1].clone() + p[2].clone();
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            3 * self.entries.len(),
        );

        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                if tab.is_lost[addr0] {
                    return;
                }

                let tab_seed_x = self.rng.random::<u64>();
                let tab_seed_y = self.rng.random::<u64>();
                let tab_seed_z = self.rng.random::<u64>();

                let mut tab_branch_x = tab.fork(Some(tab_seed_x));
                let mut tab_branch_y = tab.fork(Some(tab_seed_y));
                let mut tab_branch_z = tab.fork(Some(tab_seed_z));

                tab_branch_x.x(addr0);
                tab_branch_y.y(addr0);
                tab_branch_z.z(addr0);

                // X/Y/Z flip only phase bits, never the Pauli words, so all three
                // branches reuse the parent's word-fingerprint and derive their
                // phase/loss hash from the parent's by XORing the flipped rows.
                let hx = pauli_branch_phase_loss(tab, &tab_branch_x, phase_loss);
                let hy = pauli_branch_phase_loss(tab, &tab_branch_y, phase_loss);
                let hz = pauli_branch_phase_loss(tab, &tab_branch_z, phase_loss);
                branches.push((tab_branch_x, p_sum.clone() * p[0].clone(), word_fp, hx));
                branches.push((tab_branch_y, p_sum.clone() * p[1].clone(), word_fp, hy));
                branches.push((tab_branch_z, p_sum.clone() * p[2].clone(), word_fp, hz));

                *p_sum *= T::Coeff::one() - p_total.clone();
            });

        let needs_normalize = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
        if needs_normalize {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> TwoQubitPauliError<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [<T as Config>::Coeff; 15]) {
        let p_total: T::Coeff = p
            .iter()
            .fold(T::Coeff::zero(), |acc, prob| acc + prob.clone());
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            15 * self.entries.len(),
        );

        // The 15 non-identity two-qubit Paulis on (addr0, addr1), in the same
        // order as the probability array `p`: IX, IY, IZ, XI, XX, XY, XZ, YI,
        // YX, YY, YZ, ZI, ZX, ZY, ZZ. Encoding: 0 = I, 1 = X, 2 = Y, 3 = Z.
        //
        // `rustfmt::skip` keeps the rows grouped by the first Pauli (a readable
        // 4-wide grid); without it rustfmt repacks the tuples to fill the line
        // width and the grouping is lost.
        #[rustfmt::skip]
        const PAULI_PAIRS: [(u8, u8); 15] = [
            (0, 1), (0, 2), (0, 3),
            (1, 0), (1, 1), (1, 2), (1, 3),
            (2, 0), (2, 1), (2, 2), (2, 3),
            (3, 0), (3, 1), (3, 2), (3, 3),
        ];

        let apply = |t: &mut GeneralizedTableau<T, I, C>, op: u8, addr: usize| match op {
            1 => t.x(addr),
            2 => t.y(addr),
            3 => t.z(addr),
            _ => {}
        };

        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                if tab.is_lost[addr0] || tab.is_lost[addr1] {
                    return;
                }

                for (k, &(op0, op1)) in PAULI_PAIRS.iter().enumerate() {
                    let tab_seed = self.rng.random::<u64>();
                    let mut tab_branch = tab.fork(Some(tab_seed));
                    apply(&mut tab_branch, op0, addr0);
                    apply(&mut tab_branch, op1, addr1);
                    // X/Y/Z flips only phase bits, so the word-fingerprint is
                    // preserved and the phase/loss hash is derived incrementally.
                    let h = pauli_branch_phase_loss(tab, &tab_branch, phase_loss);
                    branches.push((tab_branch, p_sum.clone() * p[k].clone(), word_fp, h));
                }

                *p_sum *= T::Coeff::one() - p_total.clone();
            });

        let needs_normalize = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
        if needs_normalize {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> Depolarizing2<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn depolarize2(&mut self, addr0: usize, addr1: usize, p: <T as Config>::Coeff) {
        let ps: [T::Coeff; 15] = std::array::from_fn(|_| p.clone() / 15.0.into());
        self.two_qubit_pauli_error(addr0, addr1, ps);
    }
}

impl<
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
> CorrelatedLossChannel<T> for GeneralizedTableauSum<T, I, C, S>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn correlated_loss_channel(
        &mut self,
        addr0: usize,
        addr1: usize,
        p: [<T as Config>::Coeff; 3],
    ) {
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            3 * self.entries.len(),
        );
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                // if either is lost already, we just lose the other with probability p[2]
                if tab.is_lost[addr0] {
                    single_qubit_loss_branch(
                        addr1,
                        &p[2],
                        &mut self.rng,
                        &mut branches,
                        tab,
                        p_sum,
                        (word_fp, phase_loss),
                    );
                    return;
                } else if tab.is_lost[addr1] {
                    single_qubit_loss_branch(
                        addr0,
                        &p[2],
                        &mut self.rng,
                        &mut branches,
                        tab,
                        p_sum,
                        (word_fp, phase_loss),
                    );
                    return;
                }

                // if both are present, then we create 3 new branches:
                // losing both (p[0]), one, or the other qubit (p[1])

                let tab_seed_both = self.rng.random::<u64>();
                let mut tab_lose_both = tab.fork(Some(tab_seed_both));
                tab_lose_both.is_lost[addr0] = true;
                tab_lose_both.is_lost[addr1] = true;

                // is_lost flip leaves the Pauli words and phases unchanged, so
                // the branch reuses its parent's word-fingerprint and the only
                // change to the phase/loss hash is the lost qubit's mask.
                branches.push((
                    tab_lose_both,
                    p_sum.clone() * p[0].clone(),
                    word_fp,
                    phase_loss ^ loss_mask(addr0) ^ loss_mask(addr1),
                ));

                let tab_seed_0 = self.rng.random::<u64>();
                let mut tab_lose_0 = tab.fork(Some(tab_seed_0));
                tab_lose_0.is_lost[addr0] = true;
                branches.push((
                    tab_lose_0,
                    p_sum.clone() * (p[1].clone() / 2.0.into()),
                    word_fp,
                    phase_loss ^ loss_mask(addr0),
                ));

                let tab_seed_1 = self.rng.random::<u64>();
                let mut tab_lose_1 = tab.fork(Some(tab_seed_1));
                tab_lose_1.is_lost[addr1] = true;
                branches.push((
                    tab_lose_1,
                    p_sum.clone() * (p[1].clone() / 2.0.into()),
                    word_fp,
                    phase_loss ^ loss_mask(addr1),
                ));

                let p_total = p[0].clone() + p[1].clone();
                *p_sum *= T::Coeff::one() - p_total;
            });

        let needs_renormalize = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
        if needs_renormalize {
            self.normalize_probabilities();
        }
        self.truncate();
    }
}

impl<T, I, C, S> ResetLossChannel<T> for GeneralizedTableauSum<T, I, C, S>
where
    T: Config,
    I: TableauIndex + Send + Sync,
    C: SparseVector<Complex<T::Coeff>, I>,
    S: EntryStore<T, I, C>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + PartialOrd
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
{
    fn reset_loss_channel(&mut self, addr0: usize) {
        let delta = loss_mask(addr0);
        let mut branches = self.entries.drain_where(|tab| tab.is_lost[addr0]);
        for (tab, _, _, phase_loss) in branches.iter_mut() {
            tab.is_lost[addr0] = false;
            *phase_loss ^= delta;
        }
        // reset_loss preserves total probability mass and never drops entries
        // below sum_cutoff (merges only ever sum existing coefficients), so no
        // renormalize is needed regardless of insert_or_merge_batch's result.
        let _ = self
            .entries
            .insert_or_merge_batch(branches, &self.sum_cutoff);
    }
}
