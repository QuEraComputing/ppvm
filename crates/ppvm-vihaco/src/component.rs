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
use vihaco::{Apply, Component, Effects, Execute, PopFromStack, component, dispatch, observe};
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

/// Generic backend component - use this if you need one concrete variant in your machine
#[component]
#[instructions { circuit }]
pub struct CircuitExecutor<B: CircuitBackend> {
    pub backend: B,
}

// CircuitOp trait: define how an instruction executes on a generic backend
// This is used to avoid having huge where constraints on the enums later
// - without this trait, we'd have to constrain every enum variant to all Execute impls
// TODO: re-evaluate whether this is needed once we refactor and try to get rid of enums
pub trait CircuitOp {
    type Message;
    type Effect;

    fn apply<B: CircuitBackend>(
        &self,
        backend: &mut B,
        msg: Self::Message,
    ) -> eyre::Result<Self::Effect>;
}

macro_rules! single_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl CircuitOp for circuit::$leaf {
            type Effect = ();
            type Message = PopFromStack<(u64,)>;

            fn apply<B: CircuitBackend>(
                &self,
                backend: &mut B,
                msg: Self::Message,
            ) -> eyre::Result<()> {
                backend.state().$method(msg.0.0 as usize);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_execute {
    ($leaf:ident, $method:ident) => {
        impl CircuitOp for circuit::$leaf {
            type Effect = ();
            type Message = PopFromStack<(u64, u64)>;

            fn apply<B: CircuitBackend>(
                &self,
                backend: &mut B,
                msg: Self::Message,
            ) -> eyre::Result<()> {
                backend.state().$method(msg.0.0 as usize, msg.0.1 as usize);
                Ok(())
            }
        }
    };
}

macro_rules! single_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl CircuitOp for circuit::$leaf {
            type Effect = ();
            type Message = PopFromStack<(u64, f64)>;

            fn apply<B: CircuitBackend>(
                &self,
                backend: &mut B,
                msg: Self::Message,
            ) -> eyre::Result<()> {
                backend.state().$method(msg.0.0 as usize, msg.0.1);
                Ok(())
            }
        }
    };
}

macro_rules! two_qubit_rotation_execute {
    ($leaf:ident, $method:ident) => {
        impl CircuitOp for circuit::$leaf {
            type Effect = ();
            type Message = PopFromStack<(u64, u64, f64)>;

            fn apply<B: CircuitBackend>(
                &self,
                backend: &mut B,
                msg: Self::Message,
            ) -> eyre::Result<()> {
                backend
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

impl CircuitOp for circuit::U3 {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        backend
            .state()
            .u3(msg.0.0 as usize, msg.0.1, msg.0.2, msg.0.3);
        Ok(())
    }
}

impl CircuitOp for circuit::R {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        backend.state().r(msg.0.0 as usize, msg.0.1, msg.0.2);
        Ok(())
    }
}

impl CircuitOp for circuit::Leakage {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        backend.leakage(msg.0.0 as usize, msg.0.1, msg.0.2)?;
        Ok(())
    }
}

impl CircuitOp for circuit::PauliError {
    type Effect = ();
    type Message = PopFromStack<(u64, f64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        backend
            .state()
            .pauli_error(msg.0.0 as usize, [msg.0.1, msg.0.2, msg.0.3]);
        Ok(())
    }
}

impl CircuitOp for circuit::TwoQubitPauliError {
    type Effect = ();
    type Message = PopFromStack<(u64, u64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        // TODO: probably needs a custom message since we can't pop float arrays
        unimplemented!("TODO: needs a custom message")
    }
}

impl CircuitOp for circuit::CorrelatedLoss {
    type Effect = ();
    type Message = PopFromStack<(u64, u64, f64, f64, f64)>;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, msg: Self::Message) -> eyre::Result<()> {
        backend.state().correlated_loss_channel(
            msg.0.0 as usize,
            msg.0.1 as usize,
            [msg.0.2, msg.0.3, msg.0.4],
        );
        Ok(())
    }
}

impl CircuitOp for circuit::Measure {
    type Effect = Effects<CircuitOutcomeEffect>;
    type Message = PopFromStack<(u64,)>;

    fn apply<B: CircuitBackend>(
        &self,
        backend: &mut B,
        msg: Self::Message,
    ) -> eyre::Result<Self::Effect> {
        let outcome = backend.measure(msg.0.0 as usize)?;
        Ok(Effects::one(CircuitOutcomeEffect::Measurement(
            MeasurementEffect {
                measurement_results: smallvec::smallvec![outcome],
            },
        )))
    }
}

impl CircuitOp for circuit::Truncate {
    type Effect = ();
    type Message = vihaco::NoMessage;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, _: Self::Message) -> eyre::Result<()> {
        backend.truncate();
        Ok(())
    }
}

impl CircuitOp for circuit::Trace {
    type Effect = ();
    type Message = vihaco::NoMessage;

    fn apply<B: CircuitBackend>(&self, backend: &mut B, _: Self::Message) -> eyre::Result<()> {
        unimplemented!("TODO: needs a custom message");
    }
}

impl<B: CircuitBackend> vihaco::Reset for CircuitExecutor<B> {
    fn reset(&mut self) {
        self.backend.reset();
    }
}

// The single execute impl for all CircuitOps
impl<B: CircuitBackend, I: CircuitOp> Execute<I> for CircuitExecutor<B> {
    type Effect = I::Effect;
    type Message = I::Message;
    type Error = eyre::Report;

    fn execute(
        &mut self,
        instruction: &I,
        message: Self::Message,
    ) -> Result<Self::Effect, Self::Error> {
        I::apply(instruction, &mut self.backend, message)
    }
}

// TODO: enum dispatch can go once we only set the index type in tableau #204
pub(crate) enum TableauBuckets {
    Bits64(CircuitExecutor<GeneralizedTableau<Byte8F64<1>, usize, Vec<(Complex64, usize)>>>),
    Bits128(CircuitExecutor<GeneralizedTableau<Byte8F64<2>, u128, Vec<(Complex64, u128)>>>),
    Bits256(CircuitExecutor<GeneralizedTableau<Byte8F64<4>, U256, Vec<(Complex64, U256)>>>),
    Bits512(CircuitExecutor<GeneralizedTableau<Byte8F64<8>, U512, Vec<(Complex64, U512)>>>),
    Bits1024(CircuitExecutor<GeneralizedTableau<Byte8F64<16>, U1024, Vec<(Complex64, U1024)>>>),
    Bits2048(CircuitExecutor<GeneralizedTableau<Byte8F64<32>, U2048, Vec<(Complex64, U2048)>>>),
}

impl TableauBuckets {
    fn apply<I: CircuitOp>(&mut self, op: &I, msg: I::Message) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Result<Self> {
        if n_qubits <= 64 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits64(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 128 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits128(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 256 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits256(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 512 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits512(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 1024 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits1024(CircuitExecutor { backend: tab }))
        } else if n_qubits <= MAX_QUBITS {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits2048(CircuitExecutor { backend: tab }))
        } else {
            Err(eyre!(
                "cannot simulate {n_qubits} qubits: maximum is {MAX_QUBITS}"
            ))
        }
    }

    /// Same as [`TableauCircuit::new`], but seed the RNG deterministically so a
    /// shot is reproducible.
    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: f64, seed: u64) -> Result<Self> {
        macro_rules! seeded {
            ($variant:ident) => {{
                let tab = GeneralizedTableau::new_with_seed(n_qubits, coefficient_threshold, seed);
                Ok(Self::$variant(CircuitExecutor { backend: tab }))
            }};
        }
        if n_qubits <= 64 {
            seeded!(Bits64)
        } else if n_qubits <= 128 {
            seeded!(Bits128)
        } else if n_qubits <= 256 {
            seeded!(Bits256)
        } else if n_qubits <= 512 {
            seeded!(Bits512)
        } else if n_qubits <= 1024 {
            seeded!(Bits1024)
        } else if n_qubits <= MAX_QUBITS {
            seeded!(Bits2048)
        } else {
            Err(eyre!(
                "cannot simulate {n_qubits} qubits: maximum is {MAX_QUBITS}"
            ))
        }
    }
}

pub(crate) enum PauliSumBuckets {
    Bits64(CircuitExecutor<PauliSumBackend<PauliSumConfig<8>>>),
    Bits128(CircuitExecutor<PauliSumBackend<PauliSumConfig<16>>>),
    Bits256(CircuitExecutor<PauliSumBackend<PauliSumConfig<32>>>),
    Bits512(CircuitExecutor<PauliSumBackend<PauliSumConfig<64>>>),
    Bits1024(CircuitExecutor<PauliSumBackend<PauliSumConfig<128>>>),
    Bits2048(CircuitExecutor<PauliSumBackend<PauliSumConfig<256>>>),
}

impl PauliSumBuckets {
    fn apply<I: CircuitOp>(&mut self, op: &I, msg: I::Message) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term:
    /// `for (word, coef) in terms { state += (word, coef); }`. Words must
    /// already be validated against `info.n_qubits` by the caller.
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<PauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                let initial_state = state.clone();
                let psum = PauliSumBackend::<PauliSumConfig<$N>> {
                    state,
                    initial_state,
                };
                Ok(Self::$variant(CircuitExecutor { backend: psum }))
            }};
        }
        if info.n_qubits <= 64 {
            build!(Bits64, 8)
        } else if info.n_qubits <= 128 {
            build!(Bits128, 16)
        } else if info.n_qubits <= 256 {
            build!(Bits256, 32)
        } else if info.n_qubits <= 512 {
            build!(Bits512, 64)
        } else if info.n_qubits <= 1024 {
            build!(Bits1024, 128)
        } else if info.n_qubits <= MAX_QUBITS {
            build!(Bits2048, 256)
        } else {
            Err(eyre!(
                "cannot simulate {} qubits: maximum is {MAX_QUBITS}",
                info.n_qubits
            ))
        }
    }
}

pub(crate) enum LossyPauliSumBuckets {
    Bits64(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<8>>>),
    Bits128(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<16>>>),
    Bits256(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<32>>>),
    Bits512(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<64>>>),
    Bits1024(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<128>>>),
    Bits2048(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<256>>>),
}

impl LossyPauliSumBuckets {
    fn apply<I: CircuitOp>(&mut self, op: &I, msg: I::Message) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<LossyPauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                let initial_state = state.clone();
                let psum = PauliSumBackend::<LossyPauliSumConfig<$N>> {
                    state,
                    initial_state,
                };
                Ok(Self::$variant(CircuitExecutor { backend: psum }))
            }};
        }
        if info.n_qubits <= 64 {
            build!(Bits64, 8)
        } else if info.n_qubits <= 128 {
            build!(Bits128, 16)
        } else if info.n_qubits <= 256 {
            build!(Bits256, 32)
        } else if info.n_qubits <= 512 {
            build!(Bits512, 64)
        } else if info.n_qubits <= 1024 {
            build!(Bits1024, 128)
        } else if info.n_qubits <= MAX_QUBITS {
            build!(Bits2048, 256)
        } else {
            Err(eyre!(
                "cannot simulate {} qubits: maximum is {MAX_QUBITS}",
                info.n_qubits
            ))
        }
    }
}

pub(crate) enum CircuitKind {
    Tableau(TableauBuckets),
    PauliSum(PauliSumBuckets),
    LossyPauliSum(LossyPauliSumBuckets),
}

impl CircuitKind {
    fn reset(&mut self) {
        match self {
            Self::Tableau(b) => b.reset(),
            Self::PauliSum(b) => b.reset(),
            Self::LossyPauliSum(b) => b.reset(),
        }
    }
}

/// circuit component that selects qubit size from a runtime value and can do all backends
#[component]
#[instructions { circuit }]
pub struct Circuit {
    inner: CircuitKind,
}

impl<I: CircuitOp> Execute<I> for Circuit {
    type Message = I::Message;
    type Effect = I::Effect;
    type Error = eyre::Report;

    fn execute(&mut self, instruction: &I, msg: I::Message) -> Result<Self::Effect, eyre::Report> {
        match &mut self.inner {
            CircuitKind::Tableau(b) => b.apply(instruction, msg),
            CircuitKind::PauliSum(b) => b.apply(instruction, msg),
            CircuitKind::LossyPauliSum(b) => b.apply(instruction, msg),
        }
    }
}

/// Tableau-only component with runtime configured size
#[component]
#[instructions { circuit }]
pub struct TableauCircuit {
    inner: TableauBuckets,
}

impl<I: CircuitOp> Execute<I> for TableauCircuit {
    type Message = I::Message;
    type Effect = I::Effect;
    type Error = eyre::Report;

    fn execute(
        &mut self,
        instruction: &I,
        message: Self::Message,
    ) -> Result<Self::Effect, Self::Error> {
        self.inner.apply(instruction, message)
    }
}

impl vihaco::Reset for Circuit {
    fn reset(&mut self) {
        self.inner.reset()
    }
}

impl vihaco::Reset for TableauCircuit {
    fn reset(&mut self) {
        self.inner.reset()
    }
}

impl Circuit {
    /// Build a Tableau-backed circuit. Tableau init only needs `n_qubits` and
    /// `coefficient_threshold` from `info`; no observable required.
    pub fn tableau(info: &PPVMDeviceInfo) -> Result<Self> {
        let tab = TableauBuckets::new(info.n_qubits, info.coefficient_threshold)?;
        Ok(Self {
            inner: CircuitKind::Tableau(tab),
        })
    }

    /// Same as [`Circuit.tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn tableau_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        let tab = TableauBuckets::new_with_seed(info.n_qubits, info.coefficient_threshold, seed)?;
        Ok(Self {
            inner: CircuitKind::Tableau(tab),
        })
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term in
    /// `terms`. Each `(word, coef)` is added via `state += (word, coef)`; the
    /// caller is responsible for having parsed/validated the words against
    /// `info.n_qubits` (see `parse_observable_terms` in `composite.rs`).
    pub fn paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        let psum = PauliSumBuckets::new(info, terms)?;
        Ok(Self {
            inner: CircuitKind::PauliSum(psum),
        })
    }

    /// Build a LossyPauliSum-backed circuit. Same contract as
    /// [`Circuit::paulisum`].
    pub fn lossy_paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        let psum = LossyPauliSumBuckets::new(info, terms)?;
        Ok(Self {
            inner: CircuitKind::LossyPauliSum(psum),
        })
    }
}

impl Default for Circuit {
    fn default() -> Self {
        // Default backend is Tableau with 0 qubits, which always fits the
        // smallest bucket, so construction here is infallible.
        Self::tableau(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}

impl TableauCircuit {
    /// Build a Tableau-backed circuit. Tableau init only needs `n_qubits` and
    /// `coefficient_threshold` from `info`; no observable required.
    pub fn new(info: &PPVMDeviceInfo) -> Result<Self> {
        let tab = TableauBuckets::new(info.n_qubits, info.coefficient_threshold)?;
        Ok(Self { inner: tab })
    }

    /// Same as [`Circuit.tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn new_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        let tab = TableauBuckets::new_with_seed(info.n_qubits, info.coefficient_threshold, seed)?;
        Ok(Self { inner: tab })
    }
}

impl Default for TableauCircuit {
    fn default() -> Self {
        Self::new(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}
