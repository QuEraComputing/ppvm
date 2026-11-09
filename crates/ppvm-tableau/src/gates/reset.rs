// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;

use crate::prelude::*;
use bitvec::view::BitView;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, PrimInt, ToPrimitive, Zero};

impl<T: Config> ResetState for Tableau<T>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    fn reset_state(&mut self, addr0: usize) {
        let m = self.measure(addr0);
        if m {
            self.x(addr0);
        }
    }
}

/// A plain tableau has no loss or leakage, so a full reset is a state reset.
impl<T: Config> Reset for Tableau<T>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
{
    fn reset(&mut self, addr0: usize) {
        self.reset_state(addr0);
    }
}

impl<T, I, C> ResetState for GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Debug,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
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
{
    fn reset_state(&mut self, addr0: usize) {
        // Skip qubits outside the computational subspace: leaked qubits must
        // not be re-zeroed, and lost qubits have no state to reset.
        if self.is_inactive(addr0) {
            return;
        }

        let m = self.measure(addr0);

        // A reset is not a measurement in stim's model: drop the record
        // entry that the internal `measure` just pushed so the reset is
        // measurement-record-neutral.
        self.measurement_record.pop();

        if let Some(true) = m {
            self.x(addr0);
        }
    }
}

impl<T, I, C> Reset for GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Debug,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
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
{
    fn reset(&mut self, addr0: usize) {
        // Un-lose / un-leak first, otherwise `reset_state` skips the qubit.
        // `reset_leakage_channel` already re-zeroes the qubit.
        if self.qubit_status[addr0] == QubitStatus::Leaked {
            self.reset_leakage_channel(addr0);
        } else {
            self.reset_loss_channel(addr0);
            self.reset_state(addr0);
        }
    }
}
