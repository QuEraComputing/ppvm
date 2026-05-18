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

use crate::data::GeneralizedTableauSum;

impl<T: Config, I: TableauIndex + Send + Sync, C: SparseVector<Complex<T::Coeff>, I>> LossChannel<T>
    for GeneralizedTableauSum<T, I, C>
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
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff)>::new();
        for (tab, p_sum) in self.entries.iter_mut() {
            if tab.is_lost[addr0] {
                // Don't branch if it's already lost
                continue;
            }

            let tab_seed = self.rng.random::<u64>();
            let mut tab_branch = tab.fork(Some(tab_seed));
            tab_branch.is_lost[addr0] = true;
            branches.push((tab_branch, p_sum.clone() * p.clone()));
            *p_sum *= T::Coeff::one() - p.clone();
        }

        self.insert_or_update_batch(&branches);
        self.truncate();
    }
}

impl<T: Config, I: TableauIndex + Send + Sync, C: SparseVector<Complex<T::Coeff>, I>>
    Depolarizing<T> for GeneralizedTableauSum<T, I, C>
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
        let mut branches = Vec::<(GeneralizedTableau<T, I, C>, T::Coeff)>::new();
        let p_3 = p.clone() / 3.0.into();

        for (tab, p_sum) in self.entries.iter_mut() {
            if tab.is_lost[addr0] {
                continue;
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

            branches.push((tab_branch_x, p_sum.clone() * p_3.clone()));
            branches.push((tab_branch_y, p_sum.clone() * p_3.clone()));
            branches.push((tab_branch_z, p_sum.clone() * p_3.clone()));

            *p_sum *= T::Coeff::one() - p.clone();
        }

        self.insert_or_update_batch(&branches);
        self.truncate();
    }
}
