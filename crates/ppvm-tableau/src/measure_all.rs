use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_runtime::config::Config;

use crate::measure::MeasureScratch;
use crate::{data::GeneralizedTableau, sparsevec::SparseVector, tableau_index::TableauIndex};

pub trait LossyMeasureAll {
    fn measure_all(&mut self) -> Vec<Option<bool>>;
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> LossyMeasureAll
    for GeneralizedTableau<T, I, C>
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
    fn measure_all(&mut self) -> Vec<Option<bool>> {
        // One scratch reused across all n measurements: keeps the HashMap and
        // b_entries Vec allocations across qubits, and caches the
        // odd-phase-destabilizer mask between case-a measurements (only
        // recomputed when `update_tableau_according_to_outcome` modifies
        // destabilizer phases).
        let mut scratch: MeasureScratch<I, T::Coeff> = MeasureScratch::new();
        self.measure_all_with_scratch(&mut scratch)
    }
}

impl<T: Config, I, C: SparseVector<Complex<T::Coeff>, I>> GeneralizedTableau<T, I, C>
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
    /// Same as [`LossyMeasureAll::measure_all`], but the caller supplies a
    /// `MeasureScratch` that's reused across the n per-qubit measurements
    /// (and, if the caller chooses, across many invocations / shots).
    ///
    /// This is the entry point samplers should use when running many shots:
    /// initialize one scratch alongside the sampler and thread it through
    /// every shot to amortize the case-a HashMap and b_entries allocations.
    pub fn measure_all_with_scratch(
        &mut self,
        scratch: &mut MeasureScratch<I, T::Coeff>,
    ) -> Vec<Option<bool>> {
        (0..self.n_qubits())
            .map(|idx| self.measure_with_scratch(idx, scratch))
            .collect()
    }
}
