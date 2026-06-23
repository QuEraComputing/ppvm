// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;

use bitvec::view::BitView;
use num::{
    Complex, One, PrimInt, ToPrimitive, Zero,
    complex::{Complex64, ComplexFloat},
};
use ppvm_traits::{char::Pauli, config::Config};

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
            .map(|idx| self.measure_one_with_scratch(idx, scratch))
            .collect()
    }

    /// Measure the given qubit `indices` in order, reusing a caller-supplied
    /// `MeasureScratch` across the per-index measurements (and, if the caller
    /// chooses, across many invocations / shots) — the explicit-index analogue
    /// of [`measure_all_with_scratch`](Self::measure_all_with_scratch). This is
    /// the scratch-reusing engine behind the tableau's `LossyMeasure::measure_many`
    /// override (in `measure.rs`), which the Stim executor and the Python
    /// `measure_many` binding both go through.
    pub fn measure_many_with_scratch(
        &mut self,
        indices: &[usize],
        scratch: &mut MeasureScratch<I, T::Coeff>,
    ) -> Vec<Option<bool>> {
        indices
            .iter()
            .map(|&idx| self.measure_one_with_scratch(idx, scratch))
            .collect()
    }

    /// Measure a single qubit `idx` in the Z basis, reusing `scratch`. A lost
    /// qubit pushes `None` onto the measurement record and returns `None`,
    /// exactly as the standalone `measure` path does. Shared by
    /// `measure_all_with_scratch` and `measure_many_with_scratch` so every
    /// batched path measures — and records — identically to per-qubit `measure`.
    fn measure_one_with_scratch(
        &mut self,
        idx: usize,
        scratch: &mut MeasureScratch<I, T::Coeff>,
    ) -> Option<bool> {
        if self.is_lost[idx] {
            self.measurement_record.push(None);
            return None;
        }
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(idx, Pauli::Z);
        self.measure_with_scratch(
            idx,
            scratch,
            phase_decomp,
            stab_anticomm_bits,
            destab_anticomm_bits,
        )
    }
}
