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
use vihaco::{Effects, component, observe};
use vihaco_circuit_isa::{CircuitEffect, CircuitInstruction, CircuitMessage};

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

pub struct CircuitExecutor<T: Config<Coeff = f64>, I: TableauIndex, C: SparseVector<Complex64, I>> {
    pub tab: GeneralizedTableau<T, I, C>,
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = CircuitOutcomeEffect)]
impl<T, I, C> CircuitExecutor<T, I, C>
where
    T: Config<Coeff = f64>,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    I: TableauIndex + Send + Sync + std::fmt::Debug,
    C: SparseVector<Complex64, I> + std::fmt::Debug,
{
    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        use CircuitInstruction::*;
        use CircuitMessage::*;

        match (inst, msg) {
            // Single-qubit Clifford
            (X, &Qubit(addr)) => self.tab.x(addr),
            (Y, &Qubit(addr)) => self.tab.y(addr),
            (Z, &Qubit(addr)) => self.tab.z(addr),
            (H, &Qubit(addr)) => self.tab.h(addr),
            (S, &Qubit(addr)) => self.tab.s(addr),
            (SAdj, &Qubit(addr)) => self.tab.s_dag(addr),
            (SqrtX, &Qubit(addr)) => self.tab.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => self.tab.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => self.tab.sqrt_x_dag(addr),
            (SqrtYAdj, &Qubit(addr)) => self.tab.sqrt_y_dag(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => self.tab.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => self.tab.cz(addr0, addr1),

            // T gate
            (T, &Qubit(addr)) => self.tab.t(addr),
            (TAdj, &Qubit(addr)) => self.tab.t_dag(addr),

            // Single-qubit rotations
            (RX, &QubitAndFloat(addr, angle)) => self.tab.rx(addr, angle),
            (RY, &QubitAndFloat(addr, angle)) => self.tab.ry(addr, angle),
            (RZ, &QubitAndFloat(addr, angle)) => self.tab.rz(addr, angle),

            // Two-qubit rotations
            (RXX, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.rxx(addr0, addr1, angle),
            (RYY, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.ryy(addr0, addr1, angle),
            (RZZ, &TwoQubitAndFloat(addr0, addr1, angle)) => self.tab.rzz(addr0, addr1, angle),

            // U3
            (U3, &QubitU3(addr, theta, phi, lam)) => self.tab.u3(addr, theta, phi, lam),

            // RXY: rotation about an axis in the x/y plane
            (R, &QubitAndTwoFloats(addr, axis_angle, theta)) => self.tab.r(addr, axis_angle, theta),

            // Measure & Reset
            (Measure, &Qubit(addr)) => {
                let outcome: MeasurementOutcome = self.tab.measure(addr).into();
                return Ok(Effects::one(CircuitOutcomeEffect::Measurement(
                    MeasurementEffect {
                        measurement_results: smallvec::smallvec![outcome],
                    },
                )));
            }
            (Reset, &Qubit(addr)) => self.tab.reset(addr),

            // Noise
            (Depolarize, &QubitAndFloat(addr, p)) => self.tab.depolarize1(addr, p),
            (Depolarize2, &TwoQubitAndFloat(addr0, addr1, p)) => {
                self.tab.depolarize2(addr0, addr1, p)
            }
            (PauliError, QubitAndFloatArr3(addr0, ps)) => self.tab.pauli_error(*addr0, *ps),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(addr0, addr1, ps)) => {
                self.tab.two_qubit_pauli_error(*addr0, *addr1, *ps)
            }

            // Loss
            (Loss, &QubitAndFloat(addr, p)) => self.tab.loss_channel(addr, p),
            (CorrelatedLoss, TwoQubitAndFloatArr3(addr0, addr1, ps)) => {
                self.tab.correlated_loss_channel(*addr0, *addr1, *ps)
            }

            /* BATCH OPERATIONS START HERE */
            // Batch: dedicated batch methods
            (SqrtX, QubitBatch(addrs)) => self.tab.sqrt_x_many(addrs),
            (SqrtY, QubitBatch(addrs)) => self.tab.sqrt_y_many(addrs),
            (SqrtXAdj, QubitBatch(addrs)) => self.tab.sqrt_x_dag_many(addrs),
            (SqrtYAdj, QubitBatch(addrs)) => self.tab.sqrt_y_dag_many(addrs),
            (H, QubitBatch(addrs)) => self.tab.h_many(addrs),
            (CZ, TwoQubitBatch(pairs)) => self.tab.cz_many(pairs),

            // TODO: replace things below by actual batched methods once they are available
            // Batch: single-qubit for loops
            (X, QubitBatch(addrs)) => self.tab.x_many(addrs),
            (Y, QubitBatch(addrs)) => self.tab.y_many(addrs),
            (Z, QubitBatch(addrs)) => self.tab.z_many(addrs),
            (S, QubitBatch(addrs)) => self.tab.s_many(addrs),
            (SAdj, QubitBatch(addrs)) => self.tab.s_dag_many(addrs),
            (T, QubitBatch(addrs)) => self.tab.t_many(addrs),
            (TAdj, QubitBatch(addrs)) => self.tab.t_dag_many(addrs),
            (Reset, QubitBatch(addrs)) => self.tab.reset_many(addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => self.tab.rx_many(addrs, *angle),
            (RY, QubitBatchAndFloat(addrs, angle)) => self.tab.ry_many(addrs, *angle),
            (RZ, QubitBatchAndFloat(addrs, angle)) => self.tab.rz_many(addrs, *angle),
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.tab, depolarize1, addrs, *p)
            }
            (Loss, QubitBatchAndFloat(addrs, p)) => batch_for!(self.tab, loss_channel, addrs, *p),
            (PauliError, QubitBatchAndFloatArr3(addrs, ps)) => {
                batch_for!(self.tab, pauli_error, addrs, *ps)
            }
            (U3, QubitBatchU3(addrs, theta, phi, lam)) => {
                batch_for!(self.tab, u3, addrs, *theta, *phi, *lam)
            }

            // Batch: two-qubit for loops
            (CNOT, TwoQubitBatch(pairs)) => {
                self.tab.cnot_many(pairs);
            }
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                self.tab.rxx_many(pairs, *angle);
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                self.tab.ryy_many(pairs, *angle);
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                self.tab.rzz_many(pairs, *angle);
            }
            (Depolarize2, TwoQubitBatchAndFloat(pairs, p)) => {
                for &(a, b) in pairs {
                    self.tab.depolarize2(a, b, *p);
                }
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(pairs, ps)) => {
                for &(a, b) in pairs {
                    self.tab.two_qubit_pauli_error(a, b, *ps);
                }
            }
            (CorrelatedLoss, TwoQubitBatchAndFloatArr3(pairs, ps)) => {
                for &(a, b) in pairs {
                    self.tab.correlated_loss_channel(a, b, *ps);
                }
            }

            // Batch: measure (emits per qubit)
            (Measure, QubitBatch(addrs)) => {
                let outcomes = self.tab.measure_many(addrs);
                return Ok(Effects::one(CircuitOutcomeEffect::Measurement(
                    MeasurementEffect {
                        measurement_results: outcomes.iter().map(|&o| o.into()).collect(),
                    },
                )));
            }

            // Truncate is a silent no-op on the Tableau backend — the tableau's
            // gate methods already prune via the configured coefficient
            // threshold, so there's nothing extra to do here.
            (Truncate, None) => {}

            // Trace: parse the resolved pattern and compute Σ_{P matches} ⟨ψ|P|ψ⟩
            // on the tableau state. Asymmetric with the PauliSum semantics by
            // design (Decision 9): on the tableau this is a sum of expectations,
            // not a coefficient filter.
            (Trace, PauliPatternStr(s)) => {
                let pat = PauliPattern::parse(s)
                    .map_err(|e| eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                let value = self.tab.trace(&pat);
                return Ok(Effects::one(CircuitOutcomeEffect::Trace(TraceEffect {
                    value,
                })));
            }

            // Fallback
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid circuit instruction arguments {:?} for instruction {:?}",
                    msg,
                    inst
                ));
            }
        };

        Ok(Effects::None)
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

/// Shared dispatch body for `PauliSumExecutor` and `LossyPauliSumExecutor`.
/// Every non-loss `CircuitInstruction` lands here. `LossyPauliSumExecutor`
/// matches `Loss` / `CorrelatedLoss` (single + batched) before invoking this
/// macro and never reaches the loss-rejection arm below.
///
/// `$self` is passed as an `ident` (typically `self`) so the macro's
/// expansion shares hygiene with the surrounding method's `self` parameter.
/// `$inst` / `$msg` are passed the same way; `$backend` is the human-readable
/// backend name baked into error messages.
macro_rules! dispatch_common_paulisum {
    ($self:ident, $inst:ident, $msg:ident, $backend:literal) => {{
        use CircuitInstruction::*;
        use CircuitMessage::*;
        match ($inst, $msg) {
            // Single-qubit Clifford
            (X, &Qubit(addr)) => $self.state.x(addr),
            (Y, &Qubit(addr)) => $self.state.y(addr),
            (Z, &Qubit(addr)) => $self.state.z(addr),
            (H, &Qubit(addr)) => $self.state.h(addr),
            (S, &Qubit(addr)) => $self.state.s(addr),
            (SAdj, &Qubit(addr)) => $self.state.s_dag(addr),
            (SqrtX, &Qubit(addr)) => $self.state.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => $self.state.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => $self.state.sqrt_x_dag(addr),
            (SqrtYAdj, &Qubit(addr)) => $self.state.sqrt_y_dag(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => $self.state.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => $self.state.cz(addr0, addr1),

            // Single-qubit rotations
            (RX, &QubitAndFloat(addr, angle)) => $self.state.rx(addr, angle),
            (RY, &QubitAndFloat(addr, angle)) => $self.state.ry(addr, angle),
            (RZ, &QubitAndFloat(addr, angle)) => $self.state.rz(addr, angle),

            // Two-qubit rotations
            (RXX, &TwoQubitAndFloat(addr0, addr1, angle)) => {
                $self.state.rxx(addr0, addr1, angle)
            }
            (RYY, &TwoQubitAndFloat(addr0, addr1, angle)) => {
                $self.state.ryy(addr0, addr1, angle)
            }
            (RZZ, &TwoQubitAndFloat(addr0, addr1, angle)) => {
                $self.state.rzz(addr0, addr1, angle)
            }

            // RXY: rotation about an axis in the x/y plane
            (R, &QubitAndTwoFloats(addr, axis_angle, theta)) => {
                $self.state.r(addr, axis_angle, theta)
            }

            // Noise
            (Depolarize, &QubitAndFloat(addr, p)) => $self.state.depolarize1(addr, p),
            (Depolarize2, &TwoQubitAndFloat(addr0, addr1, p)) => {
                $self.state.depolarize2(addr0, addr1, p)
            }
            (PauliError, QubitAndFloatArr3(addr0, ps)) => {
                $self.state.pauli_error(*addr0, *ps)
            }
            (TwoQubitPauliError, TwoQubitAndFloatArr15(addr0, addr1, ps)) => {
                $self.state.two_qubit_pauli_error(*addr0, *addr1, *ps)
            }

            // Truncate: pruning per the configured strategy.
            (Truncate, None) => $self.state.truncate(),

            // Batched arms: simple for-loop dispatch (no dedicated batch
            // methods on PauliSum<T>, unlike GeneralizedTableau).
            (X, QubitBatch(addrs)) => $self.state.x_many(addrs),
            (Y, QubitBatch(addrs)) => $self.state.y_many(addrs),
            (Z, QubitBatch(addrs)) => $self.state.z_many(addrs),
            (H, QubitBatch(addrs)) => $self.state.h_many(addrs),
            (S, QubitBatch(addrs)) => $self.state.s_many(addrs),
            (SAdj, QubitBatch(addrs)) => $self.state.s_dag_many(addrs),
            (SqrtX, QubitBatch(addrs)) => $self.state.sqrt_x_many(addrs),
            (SqrtY, QubitBatch(addrs)) => $self.state.sqrt_y_many(addrs),
            (SqrtXAdj, QubitBatch(addrs)) => $self.state.sqrt_x_dag_many(addrs),
            (SqrtYAdj, QubitBatch(addrs)) => $self.state.sqrt_y_dag_many(addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => {
                $self.state.rx_many(addrs, *angle)
            }
            (RY, QubitBatchAndFloat(addrs, angle)) => {
                $self.state.ry_many(addrs, *angle)
            }
            (RZ, QubitBatchAndFloat(addrs, angle)) => {
                $self.state.rz_many(addrs, *angle)
            }
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!($self.state, depolarize1, addrs, *p)
            }
            (PauliError, QubitBatchAndFloatArr3(addrs, ps)) => {
                batch_for!($self.state, pauli_error, addrs, *ps)
            }
            (CNOT, TwoQubitBatch(pairs)) => $self.state.cnot_many(pairs),
            (CZ, TwoQubitBatch(pairs)) => $self.state.cz_many(pairs),
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                $self.state.rxx_many(pairs, *angle)
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                $self.state.ryy_many(pairs, *angle)
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                $self.state.rzz_many(pairs, *angle)
            }
            (Depolarize2, TwoQubitBatchAndFloat(pairs, p)) => {
                batch_pairs_for!($self.state, depolarize2, pairs, *p)
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(pairs, ps)) => {
                batch_pairs_for!($self.state, two_qubit_pauli_error, pairs, *ps)
            }

            // Not supported on either backend (Decision 11 + Gate Support
            // Matrix). Loss / CorrelatedLoss handling differs by backend
            // and lives in the caller's impl block, not this macro.
            (Measure | Reset, _) => {
                return Err(eyre!("{} is not supported on the {} backend", $inst, $backend));
            }

            // T / T_dag / U3 are listed as supported on PauliSum in the
            // plan's Gate Support Matrix, but ppvm-runtime does not yet
            // implement TGate or U3Gate for PauliSum<T> (only for
            // GeneralizedTableau).
            (T | TAdj | U3, _) => {
                return Err(eyre!(
                    "{} on {} requires upstream ppvm-runtime support that is not yet implemented",
                    $inst,
                    $backend
                ));
            }

            // Trace: parse the resolved pattern string and compute the
            // trace. Per plan Decision 9, parsing happens on every
            // execution; no module-load caching.
            (Trace, PauliPatternStr(s)) => {
                let pat = PauliPattern::parse(s)
                    .map_err(|e| eyre!("invalid Pauli pattern `{}`: {:?}", s, e))?;
                let value = $self.state.trace(&pat);
                return Ok(Effects::one(CircuitOutcomeEffect::Trace(TraceEffect {
                    value,
                })));
            }

            // Fallback (mismatched shapes, etc.)
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid circuit instruction arguments {:?} for instruction {:?} on the {} backend",
                    msg,
                    inst,
                    $backend
                ));
            }
        };
        Ok(Effects::None)
    }};
}

/// PauliSum-backed executor (Heisenberg picture). Holds a `PauliSum<T>` and
/// answers the same `CircuitInstruction` vocabulary as `CircuitExecutor`,
/// but without measurement / reset / loss support.
pub struct PauliSumExecutor<T: Config<Coeff = f64>> {
    pub state: PauliSum<T>,
    /// Snapshot of the seeded observable, restored by `reset`.
    initial: PauliSum<T>,
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = CircuitOutcomeEffect)]
impl<T> PauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
    for<'a> PauliSum<T>: Trace<'a, PauliPattern, Output = f64>,
{
    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        use CircuitInstruction::*;

        if matches!(inst, Loss | CorrelatedLoss) {
            return Err(eyre!(
                "{inst} is not supported on the PauliSum backend; use the LossyPauliSum backend instead"
            ));
        }

        dispatch_common_paulisum!(self, inst, msg, "PauliSum")
    }
}

impl<T> vihaco::Reset for PauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
    PauliSum<T>: Clone,
{
    fn reset(&mut self) {
        self.state = self.initial.clone();
    }
}

/// LossyPauliSum-backed executor. Same dispatch as `PauliSumExecutor` plus
/// `Loss` / `CorrelatedLoss` channels. The concrete `T` used by the
/// enclosing `Circuit::LossyPauliSum` variant is a `Config` whose
/// `PauliWordType` is `LossyPauliWord` (see `LossyPauliSumConfig`).
pub struct LossyPauliSumExecutor<T: Config<Coeff = f64>> {
    pub state: PauliSum<T>,
    /// Snapshot of the seeded observable, restored by `reset`.
    initial: PauliSum<T>,
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = CircuitOutcomeEffect)]
impl<T> LossyPauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
    for<'a> PauliSum<T>: Trace<'a, PauliPattern, Output = f64>,
{
    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        use CircuitInstruction::*;
        use CircuitMessage::*;

        // Loss / CorrelatedLoss are the only instructions that differ from
        // PauliSum; handle them here then delegate everything else to the
        // shared dispatch.
        match (inst, msg) {
            (Loss, &QubitAndFloat(addr, p)) => {
                self.state.loss_channel(addr, p);
                return Ok(Effects::None);
            }
            (CorrelatedLoss, TwoQubitAndFloatArr3(addr0, addr1, ps)) => {
                self.state.correlated_loss_channel(*addr0, *addr1, *ps);
                return Ok(Effects::None);
            }
            (Loss, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.state, loss_channel, addrs, *p);
                return Ok(Effects::None);
            }
            (CorrelatedLoss, TwoQubitBatchAndFloatArr3(pairs, ps)) => {
                batch_pairs_for!(self.state, correlated_loss_channel, pairs, *ps);
                return Ok(Effects::None);
            }
            _ => {}
        }

        dispatch_common_paulisum!(self, inst, msg, "LossyPauliSum")
    }
}

impl<T> vihaco::Reset for LossyPauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
    PauliSum<T>: Clone,
{
    fn reset(&mut self) {
        self.state = self.initial.clone();
    }
}

/// Tableau-backed inner enum (Schrödinger picture). Carries the six
/// size-bucketed `CircuitExecutor` variants; bucket is picked from `n_qubits`.
pub enum TableauCircuit {
    Bits64(CircuitExecutor<Byte8F64<1>, usize, Vec<(Complex64, usize)>>),
    Bits128(CircuitExecutor<Byte8F64<2>, u128, Vec<(Complex64, u128)>>),
    Bits256(CircuitExecutor<Byte8F64<4>, U256, Vec<(Complex64, U256)>>),
    Bits512(CircuitExecutor<Byte8F64<8>, U512, Vec<(Complex64, U512)>>),
    Bits1024(CircuitExecutor<Byte8F64<16>, U1024, Vec<(Complex64, U1024)>>),
    Bits2048(CircuitExecutor<Byte8F64<32>, U2048, Vec<(Complex64, U2048)>>),
}

impl TableauCircuit {
    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Result<Self> {
        if n_qubits <= 64 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits64(CircuitExecutor { tab }))
        } else if n_qubits <= 128 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits128(CircuitExecutor { tab }))
        } else if n_qubits <= 256 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits256(CircuitExecutor { tab }))
        } else if n_qubits <= 512 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits512(CircuitExecutor { tab }))
        } else if n_qubits <= 1024 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits1024(CircuitExecutor { tab }))
        } else if n_qubits <= MAX_QUBITS {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits2048(CircuitExecutor { tab }))
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
                Ok(Self::$variant(CircuitExecutor { tab }))
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

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Bits64(ex) => ex.execute_instruction(inst, msg),
            Self::Bits128(ex) => ex.execute_instruction(inst, msg),
            Self::Bits256(ex) => ex.execute_instruction(inst, msg),
            Self::Bits512(ex) => ex.execute_instruction(inst, msg),
            Self::Bits1024(ex) => ex.execute_instruction(inst, msg),
            Self::Bits2048(ex) => ex.execute_instruction(inst, msg),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => ex.tab.to_string(),
            Self::Bits128(ex) => ex.tab.to_string(),
            Self::Bits256(ex) => ex.tab.to_string(),
            Self::Bits512(ex) => ex.tab.to_string(),
            Self::Bits1024(ex) => ex.tab.to_string(),
            Self::Bits2048(ex) => ex.tab.to_string(),
        }
    }
}

impl vihaco::Reset for TableauCircuit {
    fn reset(&mut self) {
        match self {
            Self::Bits64(ex) => ex.reset(),
            Self::Bits128(ex) => ex.reset(),
            Self::Bits256(ex) => ex.reset(),
            Self::Bits512(ex) => ex.reset(),
            Self::Bits1024(ex) => ex.reset(),
            Self::Bits2048(ex) => ex.reset(),
        };
    }
}

/// PauliSum-backed inner enum (Heisenberg picture). Per Decision 7 of the plan,
/// the size buckets carry `[u8; N]`-storage `ByteFxHashF64` configs (N = 8, 16,
/// …, 256) rather than the tableau's `[u64; N]` configs; bucket labels match
/// the semantic qubit count (`Bits64` = 64 qubits) so the outer enum's dispatch
/// is uniform across backends.
pub enum PauliSumCircuit {
    Bits64(PauliSumExecutor<PauliSumConfig<8>>),
    Bits128(PauliSumExecutor<PauliSumConfig<16>>),
    Bits256(PauliSumExecutor<PauliSumConfig<32>>),
    Bits512(PauliSumExecutor<PauliSumConfig<64>>),
    Bits1024(PauliSumExecutor<PauliSumConfig<128>>),
    Bits2048(PauliSumExecutor<PauliSumConfig<256>>),
}

impl PauliSumCircuit {
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
                let initial = state.clone();
                Ok(Self::$variant(PauliSumExecutor { state, initial }))
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

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Bits64(ex) => ex.execute_instruction(inst, msg),
            Self::Bits128(ex) => ex.execute_instruction(inst, msg),
            Self::Bits256(ex) => ex.execute_instruction(inst, msg),
            Self::Bits512(ex) => ex.execute_instruction(inst, msg),
            Self::Bits1024(ex) => ex.execute_instruction(inst, msg),
            Self::Bits2048(ex) => ex.execute_instruction(inst, msg),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => ex.state.to_string(),
            Self::Bits128(ex) => ex.state.to_string(),
            Self::Bits256(ex) => ex.state.to_string(),
            Self::Bits512(ex) => ex.state.to_string(),
            Self::Bits1024(ex) => ex.state.to_string(),
            Self::Bits2048(ex) => ex.state.to_string(),
        }
    }
}

impl vihaco::Reset for PauliSumCircuit {
    fn reset(&mut self) {
        match self {
            Self::Bits64(ex) => ex.reset(),
            Self::Bits128(ex) => ex.reset(),
            Self::Bits256(ex) => ex.reset(),
            Self::Bits512(ex) => ex.reset(),
            Self::Bits1024(ex) => ex.reset(),
            Self::Bits2048(ex) => ex.reset(),
        };
    }
}

/// LossyPauliSum-backed inner enum. Identical shape to [`PauliSumCircuit`]
/// but with `LossyPauliWord`-keyed configs so loss-channel methods dispatch.
pub enum LossyPauliSumCircuit {
    Bits64(LossyPauliSumExecutor<LossyPauliSumConfig<8>>),
    Bits128(LossyPauliSumExecutor<LossyPauliSumConfig<16>>),
    Bits256(LossyPauliSumExecutor<LossyPauliSumConfig<32>>),
    Bits512(LossyPauliSumExecutor<LossyPauliSumConfig<64>>),
    Bits1024(LossyPauliSumExecutor<LossyPauliSumConfig<128>>),
    Bits2048(LossyPauliSumExecutor<LossyPauliSumConfig<256>>),
}

impl LossyPauliSumCircuit {
    /// Build a LossyPauliSum-backed circuit, seeding every term via
    /// `state += (word, coef)`. Words must already be validated against
    /// `info.n_qubits` by the caller.
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<LossyPauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                let initial = state.clone();
                Ok(Self::$variant(LossyPauliSumExecutor { state, initial }))
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

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Bits64(ex) => ex.execute_instruction(inst, msg),
            Self::Bits128(ex) => ex.execute_instruction(inst, msg),
            Self::Bits256(ex) => ex.execute_instruction(inst, msg),
            Self::Bits512(ex) => ex.execute_instruction(inst, msg),
            Self::Bits1024(ex) => ex.execute_instruction(inst, msg),
            Self::Bits2048(ex) => ex.execute_instruction(inst, msg),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => ex.state.to_string(),
            Self::Bits128(ex) => ex.state.to_string(),
            Self::Bits256(ex) => ex.state.to_string(),
            Self::Bits512(ex) => ex.state.to_string(),
            Self::Bits1024(ex) => ex.state.to_string(),
            Self::Bits2048(ex) => ex.state.to_string(),
        }
    }
}

impl vihaco::Reset for LossyPauliSumCircuit {
    fn reset(&mut self) {
        match self {
            Self::Bits64(ex) => ex.reset(),
            Self::Bits128(ex) => ex.reset(),
            Self::Bits256(ex) => ex.reset(),
            Self::Bits512(ex) => ex.reset(),
            Self::Bits1024(ex) => ex.reset(),
            Self::Bits2048(ex) => ex.reset(),
        };
    }
}

/// Outer `Circuit` enum: backend selector. Picks one of the three inner enums
/// based on `info.backend` at construction time; from there, every per-step
/// call routes outer → inner → executor.
pub enum Circuit {
    Tableau(TableauCircuit),
    PauliSum(PauliSumCircuit),
    LossyPauliSum(LossyPauliSumCircuit),
}

#[component(instruction = CircuitInstruction, message = CircuitMessage, effect = CircuitOutcomeEffect)]
impl Circuit {
    /// Build a Tableau-backed circuit. Tableau init only needs `n_qubits` and
    /// `coefficient_threshold` from `info`; no observable required.
    pub fn tableau(info: &PPVMDeviceInfo) -> Result<Self> {
        Ok(Self::Tableau(TableauCircuit::new(
            info.n_qubits,
            info.coefficient_threshold,
        )?))
    }

    /// Same as [`Circuit::tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn tableau_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Result<Self> {
        Ok(Self::Tableau(TableauCircuit::new_with_seed(
            info.n_qubits,
            info.coefficient_threshold,
            seed,
        )?))
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term in
    /// `terms`. Each `(word, coef)` is added via `state += (word, coef)`; the
    /// caller is responsible for having parsed/validated the words against
    /// `info.n_qubits` (see `parse_observable_terms` in `composite.rs`).
    pub fn paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        Ok(Self::PauliSum(PauliSumCircuit::new(info, terms)?))
    }

    /// Build a LossyPauliSum-backed circuit. Same contract as
    /// [`Circuit::paulisum`].
    pub fn lossy_paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        Ok(Self::LossyPauliSum(LossyPauliSumCircuit::new(info, terms)?))
    }

    fn execute(
        &mut self,
        inst: CircuitInstruction,
        msg: CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&inst, &msg)
    }

    fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Tableau(c) => c.execute_instruction(inst, msg),
            Self::PauliSum(c) => c.execute_instruction(inst, msg),
            Self::LossyPauliSum(c) => c.execute_instruction(inst, msg),
        }
    }

    /// Render the current state. Used by the REPL's `show` command.
    pub fn state_string(&self) -> String {
        match self {
            Self::Tableau(c) => c.state_string(),
            Self::PauliSum(c) => c.state_string(),
            Self::LossyPauliSum(c) => c.state_string(),
        }
    }
}

#[observe(CircuitEffect, effect=CircuitOutcomeEffect)]
impl Circuit {
    fn observe_circuit_effect(
        &mut self,
        effect: &CircuitEffect,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        self.execute_instruction(&effect.inst, &effect.msg)
    }
}

impl vihaco::Reset for Circuit {
    fn reset(&mut self) {
        match self {
            Self::Tableau(c) => c.reset(),
            Self::PauliSum(c) => c.reset(),
            Self::LossyPauliSum(c) => c.reset(),
        };
    }
}

impl Default for Circuit {
    fn default() -> Self {
        // Default backend is Tableau with 0 qubits, which always fits the
        // smallest bucket, so construction here is infallible.
        Self::tableau(&PPVMDeviceInfo::default()).expect("0-qubit tableau is always constructible")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn info(n_qubits: usize) -> PPVMDeviceInfo {
        PPVMDeviceInfo {
            n_qubits,
            ..Default::default()
        }
    }

    /// Dispatch a single-qubit `Measure` and read the outcome out of the
    /// returned effect.
    fn measure(circuit: &mut Circuit, addr: usize) -> MeasurementOutcome {
        let effects = circuit
            .execute_instruction(&CircuitInstruction::Measure, &CircuitMessage::Qubit(addr))
            .unwrap();
        match effects.into_iter().next() {
            Some(CircuitOutcomeEffect::Measurement(m)) => m.measurement_results[0],
            other => panic!("expected a measurement effect, got {other:?}"),
        }
    }

    /// Smoke test: a single-qubit gate dispatches to the tableau backend and
    /// flips the qubit.
    #[test]
    fn tableau_backend_x_flips_qubit() {
        let mut circuit = Circuit::tableau(&info(1)).unwrap();
        circuit
            .execute_instruction(&CircuitInstruction::X, &CircuitMessage::Qubit(0))
            .unwrap();
        assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
    }

    /// Smoke test: a two-qubit gate dispatches correctly — `X(0); CNOT(0, 1)`
    /// leaves both qubits in |1⟩.
    #[test]
    fn tableau_backend_cnot_propagates_flip() {
        let mut circuit = Circuit::tableau(&info(2)).unwrap();
        circuit
            .execute_instruction(&CircuitInstruction::X, &CircuitMessage::Qubit(0))
            .unwrap();
        circuit
            .execute_instruction(&CircuitInstruction::CNOT, &CircuitMessage::TwoQubit(0, 1))
            .unwrap();
        assert_eq!(measure(&mut circuit, 0), MeasurementOutcome::One);
        assert_eq!(measure(&mut circuit, 1), MeasurementOutcome::One);
    }

    // ─── Construction rejects more qubits than the backend ceiling ────────
    //
    // The widest executor bucket is 2048 qubits (U2048 / `[u8; 256]`); beyond
    // that there is no backing width, so the constructors return an error
    // rather than panicking — a panic would tear down the TUI that drives this.

    #[test]
    fn tableau_rejects_more_than_2048_qubits() {
        assert!(Circuit::tableau(&info(MAX_QUBITS + 1)).is_err());
    }

    #[test]
    fn tableau_with_seed_rejects_more_than_2048_qubits() {
        assert!(Circuit::tableau_with_seed(&info(MAX_QUBITS + 1), 0).is_err());
    }

    #[test]
    fn paulisum_rejects_more_than_2048_qubits() {
        assert!(Circuit::paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
    }

    #[test]
    fn lossy_paulisum_rejects_more_than_2048_qubits() {
        assert!(Circuit::lossy_paulisum(&info(MAX_QUBITS + 1), &[]).is_err());
    }

    #[test]
    fn constructs_at_the_2048_qubit_boundary() {
        // 2048 is the last valid bucket; it must still succeed.
        assert!(Circuit::tableau(&info(MAX_QUBITS)).is_ok());
    }
}
