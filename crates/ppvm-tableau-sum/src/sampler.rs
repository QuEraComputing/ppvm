use std::fmt::Debug;

#[cfg(feature = "rayon")]
use rayon::prelude::*;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;
use ppvm_tableau::{data::GeneralizedTableau, measure::MeasureScratch};
use ppvm_tableau::{sparsevec::SparseVector, tableau_index::TableauIndex};
use rand::{RngExt, rngs::SmallRng};

pub struct Sampler<T: Config, I, C: SparseVector<Complex<T::Coeff>, I> = Vec<(Complex64, I)>> {
    pub(crate) p_cumulative: Vec<T::Coeff>,
    pub entries: Vec<(GeneralizedTableau<T, I, C>, T::Coeff)>,
    pub(crate) rng: SmallRng,
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
        + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    pub fn sample(&mut self) -> Vec<Option<bool>> {
        let p = self.rng.random::<f64>();
        let idx = self.p_cumulative.iter().position(|p_| *p_ > p);
        match idx {
            Some(i) => {
                let tab_seed = self.rng.random::<u64>();
                let mut tab = self.entries[i].0.fork(Some(tab_seed));
                tab.measure_all_with_scratch(&mut self.scratch)
            }
            None => unreachable!("GeneralizedTableauSum normalization error!"),
        }
    }

    #[cfg(not(feature = "rayon"))]
    pub fn sample_shots(&mut self, n_shots: usize) -> Vec<Vec<Option<bool>>> {
        (0..n_shots).map(|_| self.sample()).collect()
    }

    #[cfg(feature = "rayon")]
    pub fn sample_shots(&mut self, n_shots: usize) -> Vec<Vec<Option<bool>>>
    where
        <T as Config>::BuildHasher: Sync,
        C: Send + Sync,
    {
        let sample_inds_and_seeds: Vec<(usize, u64)> = (0..n_shots)
            .map(|_| {
                let p = self.rng.random::<f64>();
                let idx = self
                    .p_cumulative
                    .iter()
                    .position(|p_| *p_ > p)
                    .expect("GeneralizedTableauSum normalization error!");
                (idx, self.rng.random::<u64>())
            })
            .collect();

        sample_inds_and_seeds
            .par_iter()
            .map_init(
                MeasureScratch::<I, T::Coeff>::new,
                |mut scratch, &(i, seed)| {
                    let mut tab = self.entries[i].0.fork(Some(seed));
                    tab.measure_all_with_scratch(&mut scratch)
                },
            )
            .collect()
    }
}
