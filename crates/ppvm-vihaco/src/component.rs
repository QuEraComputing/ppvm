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

#[component]
#[instructions { circuit }]
pub struct CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    pub tab: GeneralizedTableau<T, I, C>,
}

macro_rules! single_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl<T, I, C> Execute<circuit::$leaf> for CircuitExecutor<T, I, C>
        where
            T: Config<Coeff = f64>,
            <<T as Config>::Storage as BitView>::Store: PrimInt,
            I: TableauIndex + Send + Sync + std::fmt::Debug,
            C: SparseVector<Complex64, I> + std::fmt::Debug,
        {
            type Effect = ();
            type Message = PopFromStack<usize>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.tab.$method(msg.0);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl<T, I, C> Execute<circuit::$leaf> for CircuitExecutor<T, I, C>
        where
            T: Config<Coeff = f64>,
            <<T as Config>::Storage as BitView>::Store: PrimInt,
            I: TableauIndex + Send + Sync + std::fmt::Debug,
            C: SparseVector<Complex64, I> + std::fmt::Debug,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, usize)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.tab.$method(msg.0.0, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! single_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl<T, I, C> Execute<circuit::$leaf> for CircuitExecutor<T, I, C>
        where
            T: Config<Coeff = f64>,
            <<T as Config>::Storage as BitView>::Store: PrimInt,
            I: TableauIndex + Send + Sync + std::fmt::Debug,
            C: SparseVector<Complex64, I> + std::fmt::Debug,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.tab.$method(msg.0.0, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl<T, I, C> Execute<circuit::$leaf> for CircuitExecutor<T, I, C>
        where
            T: Config<Coeff = f64>,
            <<T as Config>::Storage as BitView>::Store: PrimInt,
            I: TableauIndex + Send + Sync + std::fmt::Debug,
            C: SparseVector<Complex64, I> + std::fmt::Debug,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, usize, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.tab.$method(msg.0.0, msg.0.1, msg.0.2);
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

impl<T, I, C> Execute<circuit::U3> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, f64, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::U3, msg: Self::Message) -> eyre::Result<()> {
        self.tab.u3(msg.0.0, msg.0.1, msg.0.2, msg.0.3);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::R> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::R, msg: Self::Message) -> eyre::Result<()> {
        self.tab.r(msg.0.0, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::Leakage> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, f64, f64)>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Leakage, msg: Self::Message) -> eyre::Result<()> {
        self.tab.leakage_channel(msg.0.0, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::PauliError> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, [f64; 3])>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::PauliError, msg: Self::Message) -> eyre::Result<()> {
        self.tab.pauli_error(msg.0.0, msg.0.1);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::TwoQubitPauliError> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, usize, [f64; 15])>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::TwoQubitPauliError, msg: Self::Message) -> eyre::Result<()> {
        self.tab.two_qubit_pauli_error(msg.0.0, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::CorrelatedLoss> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = PopFromStack<(usize, usize, [f64; 3])>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::CorrelatedLoss, msg: Self::Message) -> eyre::Result<()> {
        self.tab.correlated_loss_channel(msg.0.0, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::Measure> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = Effects<CircuitOutcomeEffect>;
    type Message = PopFromStack<usize>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Measure, msg: Self::Message) -> eyre::Result<Self::Effect> {
        let outcome: MeasurementOutcome = self.tab.measure(msg.0).into();
        Ok(Effects::one(CircuitOutcomeEffect::Measurement(
            MeasurementEffect {
                measurement_results: smallvec::smallvec![outcome],
            },
        )))
    }
}

// Truncate is a no-op on the tableau: gate methods already prune via the threshold.
impl<T, I, C> Execute<circuit::Truncate> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = vihaco::NoMessage;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Truncate, _: Self::Message) -> eyre::Result<()> {
        Ok(())
    }
}

impl<T, I, C> Execute<circuit::Trace> for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    type Effect = ();
    type Message = vihaco::NoMessage;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::Trace, _: Self::Message) -> eyre::Result<()> {
        unimplemented!("TODO: needs a custom message");
    }
}

impl<T, I, C> vihaco::Reset for CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    fn reset(&mut self) {
        self.tab.reset_all();
    }
}

/// PauliSum-backed executor (Heisenberg picture). Answers the same `circuit`
/// dialect as `CircuitExecutor`, but without measurement / reset / loss support.
#[component]
#[instructions { circuit }]
pub struct PauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
{
    pub state: PauliSum<T>,
    /// Snapshot of the seeded observable, restored by `reset`.
    initial: PauliSum<T>,
}

/// LossyPauliSum-backed executor. Same as `PauliSumExecutor` plus `Loss` /
/// `CorrelatedLoss`; `T` uses `LossyPauliWord` (see `LossyPauliSumConfig`).
#[component]
#[instructions { circuit }]
pub struct LossyPauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
{
    pub state: PauliSum<T>,
    /// Snapshot of the seeded observable, restored by `reset`.
    initial: PauliSum<T>,
}

macro_rules! paulisum_single_qubit_execute {
    ($exec:ident, $leaf:ident, $method:ident) => {
        impl<T> Execute<circuit::$leaf> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<usize>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.state.$method(msg.0);
                Ok(())
            }
        }
    };
}

macro_rules! paulisum_two_qubit_execute {
    ($exec:ident, $leaf:ident, $method:ident) => {
        impl<T> Execute<circuit::$leaf> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, usize)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.state.$method(msg.0.0, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! paulisum_single_qubit_rotation_execute {
    ($exec:ident, $leaf:ident, $method:ident) => {
        impl<T> Execute<circuit::$leaf> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.state.$method(msg.0.0, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! paulisum_two_qubit_rotation_execute {
    ($exec:ident, $leaf:ident, $method:ident) => {
        impl<T> Execute<circuit::$leaf> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, usize, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, msg: Self::Message) -> eyre::Result<()> {
                self.state.$method(msg.0.0, msg.0.1, msg.0.2);
                Ok(())
            }
        }
    };
}

/// Rejects an instruction the backend can't simulate. Message / effect types
/// match `CircuitExecutor` so programs stay backend-agnostic.
macro_rules! paulisum_unsupported_execute {
    ($exec:ident, $leaf:ident, $msg:ty, $effect:ty, $why:expr) => {
        impl<T> Execute<circuit::$leaf> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = $effect;
            type Message = $msg;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::$leaf, _: Self::Message) -> eyre::Result<$effect> {
                Err(eyre!($why))
            }
        }
    };
}

/// Every `Execute` impl shared by `PauliSumExecutor` and `LossyPauliSumExecutor`.
macro_rules! paulisum_common_execute {
    ($exec:ident, $backend:literal) => {
        paulisum_single_qubit_execute!($exec, X, x);
        paulisum_single_qubit_execute!($exec, Y, y);
        paulisum_single_qubit_execute!($exec, Z, z);
        paulisum_single_qubit_execute!($exec, H, h);
        paulisum_single_qubit_execute!($exec, S, s);
        paulisum_single_qubit_execute!($exec, SAdj, s_dag);
        paulisum_single_qubit_execute!($exec, SqrtX, sqrt_x);
        paulisum_single_qubit_execute!($exec, SqrtXAdj, sqrt_x_dag);
        paulisum_single_qubit_execute!($exec, SqrtY, sqrt_y);
        paulisum_single_qubit_execute!($exec, SqrtYAdj, sqrt_y_dag);

        paulisum_two_qubit_execute!($exec, CNOT, cnot);
        paulisum_two_qubit_execute!($exec, CZ, cz);

        paulisum_single_qubit_rotation_execute!($exec, RX, rx);
        paulisum_single_qubit_rotation_execute!($exec, RY, ry);
        paulisum_single_qubit_rotation_execute!($exec, RZ, rz);

        paulisum_two_qubit_rotation_execute!($exec, RXX, rxx);
        paulisum_two_qubit_rotation_execute!($exec, RYY, ryy);
        paulisum_two_qubit_rotation_execute!($exec, RZZ, rzz);

        paulisum_single_qubit_rotation_execute!($exec, Depolarize, depolarize1);
        paulisum_two_qubit_rotation_execute!($exec, Depolarize2, depolarize2);

        // Leakage is GeneralizedTableau-only (pinned |0⟩/|1⟩), not a qutrit |L⟩ channel.
        paulisum_unsupported_execute!($exec, Measure, PopFromStack<usize>, Effects<CircuitOutcomeEffect>,
            concat!("measure is not supported on the ", $backend, " backend"));
        paulisum_unsupported_execute!($exec, Reset, PopFromStack<usize>, (),
            concat!("reset is not supported on the ", $backend, " backend"));
        paulisum_unsupported_execute!($exec, Leakage, PopFromStack<(usize, f64, f64)>, (),
            concat!("leakage is not supported on the ", $backend, " backend"));

        // TODO: batched gates; need custom messages

        impl<T> Execute<circuit::T> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<usize>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::T, msg: Self::Message) -> eyre::Result<()> {
                self.state.rz(msg.0, std::f64::consts::PI / 8.0);
                Ok(())
            }
        }

        impl<T> Execute<circuit::TAdj> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<usize>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::TAdj, msg: Self::Message) -> eyre::Result<()> {
                self.state.rz(msg.0, -std::f64::consts::PI / 8.0);
                Ok(())
            }
        }

        // U3(θ, φ, λ) = RZ(φ)·RY(θ)·RZ(λ); Heisenberg picture, so the tableau's
        // rz(λ); ry(θ); rz(φ) runs reversed here.
        impl<T> Execute<circuit::U3> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, f64, f64, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::U3, msg: Self::Message) -> eyre::Result<()> {
                let (addr, theta, phi, lam) = msg.0;
                self.state.rz(addr, phi);
                self.state.ry(addr, theta);
                self.state.rz(addr, lam);
                Ok(())
            }
        }

        impl<T> Execute<circuit::R> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, f64, f64)>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::R, msg: Self::Message) -> eyre::Result<()> {
                self.state.r(msg.0.0, msg.0.1, msg.0.2);
                Ok(())
            }
        }

        impl<T> Execute<circuit::PauliError> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = PopFromStack<(usize, [f64; 3])>;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::PauliError, msg: Self::Message) -> eyre::Result<()> {
                self.state.pauli_error(msg.0.0, msg.0.1);
                Ok(())
            }
        }

        // CUSTOM MESSAGE: 2 addresses + 15 probabilities exceed vihaco's 16-element tuple limit.
        impl<T> Execute<circuit::TwoQubitPauliError> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = vihaco::NoMessage;
            type Error = eyre::Report;

            fn execute(
                &mut self,
                _: &circuit::TwoQubitPauliError,
                _: Self::Message,
            ) -> eyre::Result<()> {
                unimplemented!("TODO: needs a custom message");
            }
        }

        impl<T> Execute<circuit::Truncate> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = vihaco::NoMessage;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::Truncate, _: Self::Message) -> eyre::Result<()> {
                self.state.truncate();
                Ok(())
            }
        }

        // CUSTOM MESSAGE: the Pauli-pattern string can't come off the operand stack.
        impl<T> Execute<circuit::Trace> for $exec<T>
        where
            T: Config<Coeff = f64>,
        {
            type Effect = ();
            type Message = vihaco::NoMessage;
            type Error = eyre::Report;

            fn execute(&mut self, _: &circuit::Trace, _: Self::Message) -> eyre::Result<()> {
                unimplemented!("TODO: needs a custom message");
            }
        }
    };
}

paulisum_common_execute!(PauliSumExecutor, "PauliSum");
paulisum_common_execute!(LossyPauliSumExecutor, "LossyPauliSum");

paulisum_unsupported_execute!(PauliSumExecutor, Loss, PopFromStack<(usize, f64)>, (),
    "loss is not supported on the PauliSum backend; use the LossyPauliSum backend instead");
paulisum_unsupported_execute!(PauliSumExecutor, CorrelatedLoss, PopFromStack<(usize, usize, [f64; 3])>, (),
    "correlated_loss is not supported on the PauliSum backend; use the LossyPauliSum backend instead");

paulisum_single_qubit_rotation_execute!(LossyPauliSumExecutor, Loss, loss_channel);

impl<T> Execute<circuit::CorrelatedLoss> for LossyPauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
{
    type Effect = ();
    type Message = PopFromStack<(usize, usize, [f64; 3])>;
    type Error = eyre::Report;

    fn execute(&mut self, _: &circuit::CorrelatedLoss, msg: Self::Message) -> eyre::Result<()> {
        self.state.correlated_loss_channel(msg.0.0, msg.0.1, msg.0.2);
        Ok(())
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     fn info(n_qubits: usize) -> PPVMDeviceInfo {
//         PPVMDeviceInfo {
//             n_qubits,
//             ..Default::default()
//         }
//     }

//     /// Dispatch a single-qubit `Measure` and read the outcome out of the
//     /// returned effect.
//     fn measure(circuit: &mut Circuit, addr: usize) -> MeasurementOutcome {
//         let effects = circuit
//             .execute_instruction(&CircuitInstruction::Measure, &CircuitMessage::Qubit(addr))
//             .unwrap();
//         match effects.into_iter().next() {
//             Some(CircuitOutcomeEffect::Measurement(m)) => m.measurement_results[0],
//             other => panic!("expected a measurement effect, got {other:?}"),
//         }
//     }

//     /// Smoke test: a single-qubit gate dispatches to the tableau backend and
//     /// flips the qubit.
//     #[test]
//     fn tableau_backend_x_flips_qubit() {
//         let mut circuit = Circuit::tableau(&info(1)).unwrap();
//         circuit
//             .execute_instruction(&CircuitInstruction::X, &CircuitMessage::Qubit(0))
//             .unwrap();
//         assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
//     }

//     /// Smoke test: a two-qubit gate dispatches correctly — `X(0); CNOT(0, 1)`
//     /// leaves both qubits in |1⟩.
//     #[test]
//     fn tableau_backend_cnot_propagates_flip() {
//         let mut circuit = Circuit::tableau(&info(2)).unwrap();
//         circuit
//             .execute_instruction(&CircuitInstruction::X, &CircuitMessage::Qubit(0))
//             .unwrap();
//         circuit
//             .execute_instruction(&CircuitInstruction::CNOT, &CircuitMessage::TwoQubit(0, 1))
//             .unwrap();
//         assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
//         assert_eq!(measure(&mut circuit, 1), MeasurementOutcome::One);
//     }

//     // ─── Construction rejects more qubits than the backend ceiling ────────
//     //
//     // The widest executor bucket is 2048 qubits (U2048 / `[u8; 256]`); beyond
//     // that there is no backing width, so the constructors return an error
//     // rather than panicking — a panic would tear down the TUI that drives this.

//     #[test]
//     fn tableau_rejects_more_than_2048_qubits() {
//         assert!(Circuit::tableau(&info(MAX_QUBITS + 1)).is_err());
//     }

//     #[test]
//     fn tableau_with_seed_rejects_more_than_2048_qubits() {
//         assert!(Circuit::tableau_with_seed(&info(MAX_QUBITS + 1), 0).is_err());
//     }

//     #[test]
//     fn paulisum_rejects_more_than_2048_qubits() {
//         assert!(Circuit::paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
//     }

//     #[test]
//     fn lossy_paulisum_rejects_more_than_2048_qubits() {
//         assert!(Circuit::lossy_paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
//     }

//     #[test]
//     fn constructs_at_the_2048_qubit_boundary() {
//         // 2048 is the last valid bucket; it must still succeed.
//         assert!(Circuit::tableau(&info(MAX_QUBITS)).is_ok());
//     }
// }
