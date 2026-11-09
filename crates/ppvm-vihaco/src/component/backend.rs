// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::device_info::PPVMDeviceInfo;
use crate::measurements::MeasurementOutcome;
use bitvec::view::BitView;
use eyre::eyre;
use num::PrimInt;
use num::complex::Complex64;
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_pauli_sum::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
use ppvm_tableau::prelude::*;

/// Largest qubit count any backend can simulate. The widest size bucket is
/// backed by 2048-bit integers (`U2048` / `[u8; 256]`), so every constructor
/// rejects `n_qubits > MAX_QUBITS` rather than panicking.
pub const MAX_QUBITS: usize = 2048;

/// Truncation strategy used by every `PauliSum` / `LossyPauliSum` size bucket.
/// Coefficient-threshold pruning is always on; the Pauli-weight cap is set per
/// run from the header (defaults to `usize::MAX` = no cap).
pub(crate) type PauliSumStrategy = CombinedStrategy<CoefficientThreshold, MaxPauliWeight>;

/// `PauliSum<T>`'s `T` for the lossless backend: `[u8; N]` storage, fx hash,
/// f64 coefficients, the strategy above.
pub(crate) type PauliSumConfig<const N: usize> = ByteFxHashF64<N, PauliSumStrategy>;

/// Same as `PauliSumConfig` but with `LossyPauliWord` as the word type, so the
/// loss-channel methods are dispatchable on the resulting `PauliSum<T>`.
/// `LossyPauliWord`'s second type parameter (hasher) defaults to
/// `fxhash::FxBuildHasher`, matching `ByteFxHashF64`'s internal hasher.
pub(crate) type LossyPauliSumConfig<const N: usize> =
    ByteFxHashF64<N, PauliSumStrategy, LossyPauliWord<[u8; N]>>;

/// Build a `PauliSumStrategy` value from a `PPVMDeviceInfo`. Pulled out so the
/// six size-bucket constructors don't each repeat the strategy spelling.
pub(crate) fn paulisum_strategy(info: &PPVMDeviceInfo) -> PauliSumStrategy {
    CombinedStrategy(
        CoefficientThreshold(info.coefficient_threshold),
        MaxPauliWeight(info.max_pauli_weight.unwrap_or(usize::MAX)),
    )
}

/// Trait defining all methods required to implement the backend of the VM component
pub trait CircuitBackend {
    type Config: Config<Coeff = f64>;
    type State: Clifford
        + CliffordExtensions
        + CliffordBatch
        + CliffordExtensionsBatch
        + Reset
        + TGate<Self::Config>
        + U3Gate<Self::Config>
        + RotationOne<Self::Config>
        + RotXY<Self::Config>
        + RotationTwo<Self::Config>
        + Depolarizing<Self::Config>
        + Depolarizing2<Self::Config>
        + PauliError<Self::Config>
        + TwoQubitPauliError<Self::Config>
        + LossChannel<Self::Config>
        + CorrelatedLossChannel<Self::Config>;

    fn state(&mut self) -> &mut Self::State;

    // NOTE: methods that differ between backends get dedicated methods

    fn measure(&mut self, _addr0: usize) -> eyre::Result<MeasurementOutcome> {
        Err(eyre!("measure is not supported on this backend"))
    }

    fn leakage(&mut self, _addr0: usize, p0: f64, p1: f64) -> eyre::Result<()> {
        Err(eyre!("leakage is not supported on this backend"))
    }

    fn truncate(&mut self);

    /// `Σ_{P matches pattern} ⟨P⟩` on the current state.
    fn trace(&self, pattern: &PauliPattern) -> eyre::Result<f64>;

    fn reset(&mut self);
}

impl<T, I, C> CircuitBackend for GeneralizedTableau<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Config = T;
    type State = GeneralizedTableau<T, I, C>;

    fn state(&mut self) -> &mut Self::State {
        self
    }

    fn measure(&mut self, q: usize) -> eyre::Result<MeasurementOutcome> {
        let m = LossyMeasure::measure(self, q);
        let e = match m {
            Some(b) => {
                if b {
                    MeasurementOutcome::One
                } else {
                    MeasurementOutcome::Zero
                }
            }
            None => MeasurementOutcome::Lost,
        };
        Ok(e)
    }

    fn truncate(&mut self) {}

    fn trace(&self, pattern: &PauliPattern) -> eyre::Result<f64> {
        Ok(GeneralizedTableau::trace(self, pattern))
    }

    fn reset(&mut self) {
        self.reset_all();
    }

    fn leakage(&mut self, _addr0: usize, p0: f64, p1: f64) -> eyre::Result<()> {
        self.leakage_channel(_addr0, p0, p1);
        Ok(())
    }
}

/// PauliSum that also stores the initial observable so we can reset the component
pub struct PauliSumBackend<T: Config<Coeff = f64>> {
    pub state: PauliSum<T>,
    pub(crate) initial_state: PauliSum<T>,
}

impl<T: Config<Coeff = f64>> PauliSumBackend<T> {
    pub fn new(observable: PauliSum<T>) -> Self {
        Self {
            state: observable.clone(),
            initial_state: observable.clone(),
        }
    }
}

impl<T> CircuitBackend for PauliSumBackend<T>
where
    T: Config<Coeff = f64>,
    for<'a> PauliSum<T>: Trace<'a, PauliPattern, Output = f64>,
{
    type Config = T;
    type State = PauliSum<T>;

    fn state(&mut self) -> &mut Self::State {
        &mut self.state
    }

    fn reset(&mut self) {
        self.state = self.initial_state.clone();
    }

    fn truncate(&mut self) {
        self.state.truncate();
    }

    fn trace(&self, pattern: &PauliPattern) -> eyre::Result<f64> {
        Ok(Trace::trace(&self.state, pattern))
    }
}
