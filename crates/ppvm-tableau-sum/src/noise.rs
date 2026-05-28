use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::{
    config::Config,
    traits::{
        Clifford, Depolarizing, Depolarizing2, LossChannel, PauliError, ResetLossChannel,
        TwoQubitPauliError,
    },
};
use ppvm_tableau::{
    data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex,
};
use rand::RngExt;

use crate::{
    data::GeneralizedTableauSum,
    storage::{EntryStore, loss_mask, pauli_branch_phase_loss},
};

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
                if tab.is_lost[addr0] {
                    // Don't branch if it's already lost
                    return;
                }

                let tab_seed = self.rng.random::<u64>();
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

        // Non-identity two-qubit Pauli pairs on (addr0, addr1), in the same
        // order as the probability array: IX, IY, IZ, XI, XX, XY, XZ, YI,
        // YX, YY, YZ, ZI, ZX, ZY, ZZ. Encoding: 0 = I, 1 = X, 2 = Y, 3 = Z.
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
        if self.entries.reset_loss_and_merge(addr0) {
            self.normalize_probabilities();
        }
    }
}
