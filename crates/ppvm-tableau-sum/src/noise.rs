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
    storage::{Branch, BranchMutation, EntryStore, RowMasks, loss_mask, pauli_branch_phase_loss},
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
        // Lazy branch materialization: describe each loss branch as a mutation
        // of its parent entry. The merge clones the parent only when the branch
        // survives as a NEW entry; merges/below-cutoff drops never clone.
        let mut branches =
            Vec::<(usize, BranchMutation, T::Coeff, u64, u64)>::with_capacity(self.entries.len());
        let mut idx = 0usize;
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                // Increment for EVERY entry, before the lost check, so
                // parent_idx aligns with for_each_mut_with_keys' order.
                let parent_idx = idx;
                idx += 1;
                if tab.is_lost[addr0] {
                    return;
                }
                branches.push((
                    parent_idx,
                    BranchMutation::Loss { q: addr0 },
                    p_sum.clone() * p.clone(),
                    word_fp,
                    phase_loss ^ loss_mask(addr0),
                ));
                *p_sum *= T::Coeff::one() - p.clone();
            });

        let needs_renormalize = self
            .entries
            .insert_or_merge_mutated_branches(branches, &self.sum_cutoff);
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
    fn depolarize1(&mut self, addr0: usize, p: T::Coeff) {
        let p_3 = p / 3.0.into();
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
        // Lazy branch materialization: describe each X/Y/Z branch as a Pauli
        // mutation of its parent. The phase/loss delta is computed by walking the
        // parent's column once (no clone) — X flips rows with z, Y with x^z, Z
        // with x — matching what `pauli_branch_phase_loss` would produce.
        let mut branches = Vec::<(usize, BranchMutation, T::Coeff, u64, u64)>::with_capacity(
            3 * self.entries.len(),
        );
        // Precompute the per-row sign masks once instead of recomputing the
        // splitmix `sign_mask` per row per entry in the hot loop below.
        let masks = RowMasks::new(self.n_qubits);
        let mut idx = 0usize;
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                let parent_idx = idx;
                idx += 1;
                if tab.is_lost[addr0] {
                    return;
                }

                let (mut dx, mut dy, mut dz) = (0u64, 0u64, 0u64);
                for (row, pw) in tab.tableau.data.iter().enumerate() {
                    let x: bool = pw.word.xbits[addr0];
                    let z: bool = pw.word.zbits[addr0];
                    let m = masks.sign[row];
                    if z {
                        dx ^= m;
                    }
                    if x ^ z {
                        dy ^= m;
                    }
                    if x {
                        dz ^= m;
                    }
                }

                branches.push((
                    parent_idx,
                    BranchMutation::Pauli { op: 1, addr0 },
                    p_sum.clone() * p[0].clone(),
                    word_fp,
                    phase_loss ^ dx,
                ));
                branches.push((
                    parent_idx,
                    BranchMutation::Pauli { op: 2, addr0 },
                    p_sum.clone() * p[1].clone(),
                    word_fp,
                    phase_loss ^ dy,
                ));
                branches.push((
                    parent_idx,
                    BranchMutation::Pauli { op: 3, addr0 },
                    p_sum.clone() * p[2].clone(),
                    word_fp,
                    phase_loss ^ dz,
                ));

                *p_sum *= T::Coeff::one() - p_total.clone();
            });

        let needs_normalize = self
            .entries
            .insert_or_merge_mutated_branches(branches, &self.sum_cutoff);
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
