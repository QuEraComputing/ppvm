use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::measure::MeasureScratch;
use ppvm_tableau::{sparsevec::SparseVector, tableau_index::TableauIndex};
use rand::RngExt;

use crate::data::GeneralizedTableauSum;

pub struct Sampler<T: Config, I, C: SparseVector<Complex<T::Coeff>, I> = Vec<(Complex64, I)>> {
    pub(crate) p_cumulative: Vec<T::Coeff>,
    pub generalized_tableau_sum: GeneralizedTableauSum<T, I, C>,
    /// Per-thread scratch buffers reused across all shots taken on this
    /// sampler. Keeps the case-a HashMap and b_entries Vec allocations off
    /// the per-shot critical path.
    pub(crate) scratch: MeasureScratch<I, T::Coeff>,
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> Sampler<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug + Send + Sync,
{
    pub fn sample(&mut self) -> Vec<Option<bool>> {
        let p = self.generalized_tableau_sum.rng.random::<f64>();
        let idx = self.p_cumulative.iter().position(|p_| *p_ > p);
        match idx {
            Some(i) => {
                let tab_seed = self.generalized_tableau_sum.rng.random::<u64>();
                let mut tab = self.generalized_tableau_sum.entries[i]
                    .0
                    .fork(Some(tab_seed));
                tab.measure_all_with_scratch(&mut self.scratch)
            }
            None => unreachable!("GeneralizedTableauSum normalization error!"),
        }
    }

    pub fn sample_shots(&mut self, n_shots: usize) -> Vec<Vec<Option<bool>>> {
        (0..n_shots).map(|_| self.sample()).collect()
    }
}
