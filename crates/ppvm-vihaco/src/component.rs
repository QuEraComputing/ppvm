// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::device_info::PPVMDeviceInfo;
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementEffect, MeasurementOutcome, TraceEffect,
};
use bitvec::view::BitView;
use bnum::types::{U256, U512, U1024, U2048};
use eyre::{Result, eyre};
use num::PrimInt;
use num::complex::Complex64;
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;
use ppvm_pauli_sum::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
use ppvm_tableau::prelude::*;
use vihaco::{Effects, Execute, PopFromStack, component, dispatch, observe};
use vihaco_circuit_isa::circuit;

/// Largest qubit count any backend can simulate. The widest size bucket is
/// backed by 2048-bit integers (`U2048` / `[u8; 256]`), so every constructor
/// rejects `n_qubits > MAX_QUBITS` rather than panicking.
pub const MAX_QUBITS: usize = 2048;

/// Truncation strategy used by every `PauliSum` / `LossyPauliSum` size bucket.
/// Coefficient-threshold pruning is always on; the Pauli-weight cap is set per
/// run from the header (defaults to `usize::MAX` = no cap).
type PauliSumStrategy = CombinedStrategy<CoefficientThreshold, MaxPauliWeight>;

/// `PauliSum<T>`'s `T` for the lossless backend: `[u8; N]` storage, fx hash,
/// f64 coefficients, the strategy above.
type PauliSumConfig<const N: usize> = ByteFxHashF64<N, PauliSumStrategy>;

/// Same as `PauliSumConfig` but with `LossyPauliWord` as the word type, so the
/// loss-channel methods are dispatchable on the resulting `PauliSum<T>`.
/// `LossyPauliWord`'s second type parameter (hasher) defaults to
/// `fxhash::FxBuildHasher`, matching `ByteFxHashF64`'s internal hasher.
type LossyPauliSumConfig<const N: usize> =
    ByteFxHashF64<N, PauliSumStrategy, LossyPauliWord<[u8; N]>>;

/// Build a `PauliSumStrategy` value from a `PPVMDeviceInfo`. Pulled out so the
/// six size-bucket constructors don't each repeat the strategy spelling.
fn paulisum_strategy(info: &PPVMDeviceInfo) -> PauliSumStrategy {
    CombinedStrategy(
        CoefficientThreshold(info.coefficient_threshold),
        MaxPauliWeight(info.max_pauli_weight.unwrap_or(usize::MAX)),
    )
}

macro_rules! batch_for {
    ($tab:expr, $method:ident, $addrs:expr) => {
        for addr in $addrs { $tab.$method(*addr); }
    };
    ($tab:expr, $method:ident, $addrs:expr, $($arg:expr),+) => {
        for addr in $addrs { $tab.$method(*addr, $($arg),+); }
    };
}

/// Two-qubit sibling of [`batch_for!`]: drives a method that takes two qubit
/// addresses (plus optional extra args) over a slice of `(usize, usize)` pairs.
macro_rules! batch_pairs_for {
    ($state:expr, $method:ident, $pairs:expr) => {
        for &(a, b) in $pairs { $state.$method(a, b); }
    };
    ($state:expr, $method:ident, $pairs:expr, $($arg:expr),+) => {
        for &(a, b) in $pairs { $state.$method(a, b, $($arg),+); }
    };
}

/// Trait defining all methods required to implement the backend of the VM component
pub trait CircuitBackend {
    type Config: Config<Coeff = f64>;
    type State: Clifford
        + CliffordExtensions
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
    initial_state: PauliSum<T>,
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
}

#[component]
#[instructions { circuit }]
pub struct CircuitExecutor<B: CircuitBackend> {
    pub backend: B,
}

macro_rules! single_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl<B: CircuitBackend> Execute<circuit::$leaf> for CircuitExecutor<B> {
            type Effect = ();
            type Message = PopFromStack<(u64,)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.backend.state().$method(msg.0.0 as usize);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl<B: CircuitBackend> Execute<circuit::$leaf> for CircuitExecutor<B> {
            type Effect = ();
            type Message = PopFromStack<(u64, u64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.backend
                    .state()
                    .$method(msg.0.0 as usize, msg.0.1 as usize);
                Ok(())
            }
        }
    };
}

macro_rules! single_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl<B: CircuitBackend> Execute<circuit::$leaf> for CircuitExecutor<B> {
            type Effect = ();
            type Message = PopFromStack<(u64, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.backend.state().$method(msg.0.0 as usize, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl<B: CircuitBackend> Execute<circuit::$leaf> for CircuitExecutor<B> {
            type Effect = ();
            type Message = PopFromStack<(u64, u64, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.backend
                    .state()
                    .$method(msg.0.0 as usize, msg.0.1 as usize, msg.0.2);
                Ok(())
            }
        }
    };
}

single_qubit_execute!(X, x);
single_qubit_execute!(Y, y);
single_qubit_execute!(Z, z);
single_qubit_execute!(T, t);
single_qubit_execute!(S, s);
single_qubit_execute!(TAdj, t_dag);
single_qubit_execute!(SAdj, s_dag);
single_qubit_execute!(SqrtX, sqrt_x);
single_qubit_execute!(H, h);
single_qubit_execute!(SqrtXAdj, sqrt_x_dag);
single_qubit_execute!(SqrtY, sqrt_y);
single_qubit_execute!(SqrtYAdj, sqrt_y_dag);

two_qubit_execute!(CNOT, cnot);
two_qubit_execute!(CZ, cz);

single_qubit_rotation_execute!(RX, rx);
single_qubit_rotation_execute!(RY, ry);
single_qubit_rotation_execute!(RZ, rz);

two_qubit_rotation_execute!(RXX, rxx);
two_qubit_rotation_execute!(RYY, ryy);
two_qubit_rotation_execute!(RZZ, rzz);

single_qubit_execute!(Reset, reset);
single_qubit_rotation_execute!(Depolarize, depolarize1);
single_qubit_rotation_execute!(Loss, loss_channel);
two_qubit_rotation_execute!(Depolarize2, depolarize2);

// TODO: batched gates; need custom messages

impl<B: CircuitBackend> Execute<circuit::U3> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::U3, msg: Self::Message) -> eyre::Result<()> {
        self.backend
            .state()
            .u3(msg.0.0 as usize, msg.0.1, msg.0.2, msg.0.3);
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::R> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::R, msg: Self::Message) -> eyre::Result<()> {
        self.backend.state().r(msg.0.0 as usize, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::Leakage> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Leakage, msg: Self::Message) -> eyre::Result<()> {
        self.backend.leakage(msg.0.0 as usize, msg.0.1, msg.0.2)?;
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::PauliError> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::PauliError, msg: Self::Message) -> eyre::Result<()> {
        self.backend
            .state()
            .pauli_error(msg.0.0 as usize, [msg.0.1, msg.0.2, msg.0.3]);
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::TwoQubitPauliError> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, u64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::TwoQubitPauliError, msg: Self::Message) -> eyre::Result<()> {
        // TODO: probably needs a custom message since we can't pop float arrays
        unimplemented!("TODO: needs a custom message")
    }
}

impl<B: CircuitBackend> Execute<circuit::CorrelatedLoss> for CircuitExecutor<B> {
    type Effect = ();
    type Message = PopFromStack<(u64, u64, f64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::CorrelatedLoss, msg: Self::Message) -> eyre::Result<()> {
        self.backend.state().correlated_loss_channel(
            msg.0.0 as usize,
            msg.0.1 as usize,
            [msg.0.2, msg.0.3, msg.0.4],
        );
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::Measure> for CircuitExecutor<B> {
    type Effect = Effects<CircuitOutcomeEffect>;
    type Message = PopFromStack<(u64,)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Measure, msg: Self::Message) -> eyre::Result<Self::Effect> {
        let outcome = self.backend.measure(msg.0.0 as usize)?;
        Ok(Effects::one(CircuitOutcomeEffect::Measurement(
            MeasurementEffect {
                measurement_results: smallvec::smallvec![outcome],
            },
        )))
    }
}

impl<B: CircuitBackend> Execute<circuit::Truncate> for CircuitExecutor<B> {
    type Effect = ();
    type Message = vihaco::NoMessage;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Truncate, _: Self::Message) -> eyre::Result<()> {
        self.backend.truncate();
        Ok(())
    }
}

impl<B: CircuitBackend> Execute<circuit::Trace> for CircuitExecutor<B> {
    type Effect = ();
    type Message = vihaco::NoMessage;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Trace, _: Self::Message) -> eyre::Result<()> {
        unimplemented!("TODO: needs a custom message");
    }
}

impl<B: CircuitBackend> vihaco::Reset for CircuitExecutor<B> {
    fn reset(&mut self) {
        self.backend.reset();
    }
}
