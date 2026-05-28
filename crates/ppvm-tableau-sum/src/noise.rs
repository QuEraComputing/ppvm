use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::{
    config::Config,
    traits::{Clifford, Depolarizing, LossChannel},
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
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff, u64, u64)>::with_capacity(
            3 * self.entries.len(),
        );
        let p_3 = p.clone() / 3.0.into();

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
                branches.push((tab_branch_x, p_sum.clone() * p_3.clone(), word_fp, hx));
                branches.push((tab_branch_y, p_sum.clone() * p_3.clone(), word_fp, hy));
                branches.push((tab_branch_z, p_sum.clone() * p_3.clone(), word_fp, hz));

                *p_sum *= T::Coeff::one() - p.clone();
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
