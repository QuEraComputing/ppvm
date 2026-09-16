// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;

use bitvec::view::{BitView, BitViewSized};
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_pauli_word::pattern::NotIdentity;
use ppvm_tableau::{
    data::{GeneralizedTableau, QubitStatus},
    noise::is_admissible_correlated_loss,
    sparsevec::SparseVector,
    tableau_index::TableauIndex,
};
use ppvm_traits::config::Config;
use ppvm_traits::traits::{
    Clifford, CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel, PauliError,
    ResetLossChannel, TwoQubitPauliError,
};
use rand::{RngExt, rngs::SmallRng};

use crate::{
    data::GeneralizedTableauSum,
    storage::{
        Branch, BranchMutation, EntryStore, RowMasks, bit_at, loss_mask, pauli_branch_phase_loss,
    },
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
    if tab.is_lost(addr0) {
        // Don't branch if it's already lost
        return;
    }

    let tab_seed = rng.random::<u64>();
    let mut tab_branch = tab.fork(Some(tab_seed));
    tab_branch.qubit_status[addr0] = QubitStatus::Lost;
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
                if tab.is_lost(addr0) {
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
        // The store-word index / bit position of column `addr0` are the same for
        // every entry and row, so resolve them once (Lsb0 convention).
        let bits_per_word = std::mem::size_of::<<T::Storage as BitView>::Store>() * 8;
        let word_idx = addr0 / bits_per_word;
        let bit = addr0 % bits_per_word;
        let mut idx = 0usize;
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                let parent_idx = idx;
                idx += 1;
                if tab.is_lost(addr0) {
                    return;
                }

                let (mut dx, mut dy, mut dz) = (0u64, 0u64, 0u64);
                for (row, pw) in tab.tableau.data.iter().enumerate() {
                    let xw = pw.word.xbits.data.as_raw_slice();
                    let zw = pw.word.zbits.data.as_raw_slice();
                    let x: bool = bit_at(xw, word_idx, bit);
                    let z: bool = bit_at(zw, word_idx, bit);
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
                    BranchMutation::Pauli {
                        op: NotIdentity::X,
                        addr0,
                    },
                    p_sum.clone() * p[0].clone(),
                    word_fp,
                    phase_loss ^ dx,
                ));
                branches.push((
                    parent_idx,
                    BranchMutation::Pauli {
                        op: NotIdentity::Y,
                        addr0,
                    },
                    p_sum.clone() * p[1].clone(),
                    word_fp,
                    phase_loss ^ dy,
                ));
                branches.push((
                    parent_idx,
                    BranchMutation::Pauli {
                        op: NotIdentity::Z,
                        addr0,
                    },
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
                if tab.is_lost(addr0) || tab.is_lost(addr1) {
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
        debug_assert!(
            is_admissible_correlated_loss(&p),
            "correlated loss needs p0, p1 >= 0, p0 + 2*p1 <= 1, p2 in [0, 1]; got {p:?}"
        );
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            3 * self.entries.len(),
        );
        self.entries
            .for_each_mut_with_keys(|tab, p_sum, word_fp, phase_loss| {
                // if either is lost already, we just lose the other with probability p[2]
                if tab.is_lost(addr0) {
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
                } else if tab.is_lost(addr1) {
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
                // losing both (p[0]), one, or the other qubit. `p[1]` is the
                // probability that a *named* one of the pair is lost, so each
                // single-loss branch carries `p[1]` and the survivor keeps
                // `1 − p[0] − 2·p[1]`. See
                // `ppvm_traits::traits::CorrelatedLossChannel`.

                let tab_seed_both = self.rng.random::<u64>();
                let mut tab_lose_both = tab.fork(Some(tab_seed_both));
                tab_lose_both.qubit_status[addr0] = QubitStatus::Lost;
                tab_lose_both.qubit_status[addr1] = QubitStatus::Lost;

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
                tab_lose_0.qubit_status[addr0] = QubitStatus::Lost;
                branches.push((
                    tab_lose_0,
                    p_sum.clone() * p[1].clone(),
                    word_fp,
                    phase_loss ^ loss_mask(addr0),
                ));

                let tab_seed_1 = self.rng.random::<u64>();
                let mut tab_lose_1 = tab.fork(Some(tab_seed_1));
                tab_lose_1.qubit_status[addr1] = QubitStatus::Lost;
                branches.push((
                    tab_lose_1,
                    p_sum.clone() * p[1].clone(),
                    word_fp,
                    phase_loss ^ loss_mask(addr1),
                ));

                let p_total = p[0].clone() + p[1].clone() + p[1].clone();
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
        let mut branches = self.entries.drain_where(|tab| tab.is_lost(addr0));
        for (tab, _, _, phase_loss) in branches.iter_mut() {
            tab.qubit_status[addr0] = QubitStatus::Live;
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

#[cfg(test)]
mod tests {
    // === G-040 — the correlated-loss `p[1]` convention ===
    //
    // The paper (`ppvm-paper/main.tex:462`, `:523`, `:845`) is the definition of
    // record: `p[1]` is `p_LQ`, the probability that a **named** one of the pair
    // is lost, so each single-loss branch carries `p[1]`, the total weight on
    // "exactly one lost" is `2·p[1]`, and the survivor keeps `1 − p[0] − 2·p[1]`.
    // The same number is observable three ways — as a Heisenberg coefficient
    // (`ppvm-pauli-sum`), as a branch weight (this mixture) and as a sampling
    // frequency (`ppvm-tableau`'s trajectory) — and the three must agree, since
    // the cross-backend disagreement is what a Python user actually hits.

    use ppvm_pauli_sum::config::fxhash::ByteF64;
    use ppvm_pauli_sum::prelude::*;
    use ppvm_tableau::prelude::*;
    use ppvm_traits::traits::{CorrelatedLossChannel, NoStrategy};

    use crate::data::GeneralizedTableauSum;
    use crate::storage::EntryStore;

    type Cfg = ByteF64<1>;
    type TabSum = GeneralizedTableauSum<Cfg, u128>;
    type Tab = GeneralizedTableau<Cfg, u128>;
    type LossyTestSum = PauliSum<
        ppvm_pauli_sum::config::fxhash::Byte<
            1,
            f64,
            NoStrategy,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    >;

    /// Total mixture weight on branches with exactly one lost qubit.
    /// `sum_cutoff = 0.0`, so no branch is truncated away and no
    /// renormalization can hide a mis-weighted survivor.
    fn mixture_single_loss_weight(p: [f64; 3]) -> f64 {
        let mut sum: TabSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.0, 7);
        sum.correlated_loss_channel(0, 1, p);
        sum.entries
            .iter()
            .filter(|(tab, _)| tab.is_lost(0) ^ tab.is_lost(1))
            .map(|(_, probability)| *probability)
            .sum()
    }

    /// Total mixture weight on branches where both qubits are still present.
    fn mixture_survivor_weight(p: [f64; 3]) -> f64 {
        let mut sum: TabSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.0, 7);
        sum.correlated_loss_channel(0, 1, p);
        sum.entries
            .iter()
            .filter(|(tab, _)| !tab.is_lost(0) && !tab.is_lost(1))
            .map(|(_, probability)| *probability)
            .sum()
    }

    /// The Heisenberg scale factor `ppvm-pauli-sum` applies to a fully
    /// in-subspace observable, i.e. `1 − p[0] − P(exactly one lost)`.
    fn pauli_sum_survivor(p: [f64; 3]) -> f64 {
        let mut sum = LossyTestSum::builder().n_qubits(2).build();
        sum += ("ZZ", 1.0);
        sum.correlated_loss_channel(0, 1, p);
        let zz: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "ZZ".into();
        *sum.data()
            .get(&zz)
            .expect("the all-present term survives a pure rescale")
    }

    /// The trajectory's sampled fraction of runs that lose exactly one qubit.
    fn trajectory_single_loss_fraction(p: [f64; 3], trials: u64) -> f64 {
        let mut hits = 0u64;
        for seed in 0..trials {
            let mut tab: Tab = GeneralizedTableau::new_with_seed(2, 1e-12, seed);
            tab.correlated_loss_channel(0, 1, p);
            if tab.is_lost(0) ^ tab.is_lost(1) {
                hits += 1;
            }
        }
        hits as f64 / trials as f64
    }

    #[test]
    fn correlated_loss_exactly_one_lost_is_two_p1_on_every_backend() {
        let p1 = 0.3_f64;
        let p = [0.0, p1, 0.0];
        let expected = 2.0 * p1;

        let mixture = mixture_single_loss_weight(p);
        let pauli_sum = 1.0 - pauli_sum_survivor(p);
        let trajectory = trajectory_single_loss_fraction(p, 20_000);

        // Report every dissenting backend, so a failure names the split rather
        // than only its first symptom.
        let mut wrong = Vec::new();
        for (backend, value, tolerance) in [
            ("mixture", mixture, 1e-12),
            ("ppvm-pauli-sum", pauli_sum, 1e-12),
            ("trajectory (20k seeds)", trajectory, 0.02),
        ] {
            if (value - expected).abs() >= tolerance {
                wrong.push(format!("  {backend}: {value}"));
            }
        }
        assert!(
            wrong.is_empty(),
            "p = {p:?}: P(exactly one lost) must be 2*p[1] = {expected} on every \
             backend, got\n{}",
            wrong.join("\n")
        );
    }

    #[test]
    fn correlated_loss_survivor_weight_is_one_minus_p0_minus_two_p1() {
        let p = [0.2_f64, 0.3, 0.4];
        let expected = 1.0 - p[0] - 2.0 * p[1];
        let mixture = mixture_survivor_weight(p);
        assert!(
            (mixture - expected).abs() < 1e-12,
            "mixture survivor weight {mixture}, want 1 - p0 - 2*p1 = {expected}"
        );
        let pauli_sum = pauli_sum_survivor(p);
        assert!(
            (pauli_sum - expected).abs() < 1e-12,
            "ppvm-pauli-sum survivor {pauli_sum}, want {expected}"
        );
    }

    // G-043 — the admissible region `p0, p1 >= 0`, `p0 + 2·p1 <= 1`,
    // `p2 ∈ [0, 1]` is what makes the channel completely positive. Outside it
    // the mixture truncates a negative survivor weight and renormalizes,
    // silently.

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "p0 + 2*p1 <= 1")]
    fn correlated_loss_rejects_inadmissible_probabilities() {
        let mut sum: TabSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.0, 7);
        sum.correlated_loss_channel(0, 1, [0.6, 0.6, 0.0]);
    }

    /// The saturated boundary `p0 + 2·p1 == 1` is admissible and must not trip
    /// the guard.
    #[test]
    fn correlated_loss_saturated_boundary_is_admissible() {
        let mut sum: TabSum = GeneralizedTableauSum::new_with_seed(2, 1e-12, 0.0, 7);
        sum.correlated_loss_channel(0, 1, [0.2, 0.4, 1.0]);
        let survivor = mixture_survivor_weight([0.2, 0.4, 1.0]);
        assert!(survivor.abs() < 1e-15, "survivor {survivor} should vanish");
    }
}
