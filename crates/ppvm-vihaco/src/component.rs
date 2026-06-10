// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::composite::PPVMDeviceInfo;
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementEffect, MeasurementOutcome, TraceEffect,
};
use bitvec::view::BitView;
use bnum::types::{U256, U512, U1024, U2048};
use eyre::{Result, eyre};
use num::PrimInt;
use num::complex::Complex64;
use ppvm_runtime::config::fx64hash::Byte8F64;
use ppvm_runtime::config::indexmap::ByteFxHashF64;
use ppvm_runtime::strategy::{CoefficientThreshold, CombinedStrategy, MaxPauliWeight};
use ppvm_tableau::prelude::*;
use vihaco::{Effects, component, observe};
use vihaco_circuit_isa::{CircuitEffect, CircuitInstruction, CircuitMessage};

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
            (SAdj, &Qubit(addr)) => self.tab.s_adj(addr),
            (SqrtX, &Qubit(addr)) => self.tab.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => self.tab.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => self.tab.sqrt_x_adj(addr),
            (SqrtYAdj, &Qubit(addr)) => self.tab.sqrt_y_adj(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => self.tab.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => self.tab.cz(addr0, addr1),

            // T gate
            (T, &Qubit(addr)) => self.tab.t(addr),
            (TAdj, &Qubit(addr)) => self.tab.t_adj(addr),

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
            (Depolarize, &QubitAndFloat(addr, p)) => self.tab.depolarize(addr, p),
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
            (SqrtX, QubitBatch(addrs)) => self.tab.sqrt_x_batch(addrs),
            (SqrtY, QubitBatch(addrs)) => self.tab.sqrt_y_batch(addrs),
            (SqrtXAdj, QubitBatch(addrs)) => self.tab.sqrt_x_adj_batch(addrs),
            (SqrtYAdj, QubitBatch(addrs)) => self.tab.sqrt_y_adj_batch(addrs),
            (H, QubitBatch(addrs)) => self.tab.h_batch(addrs),
            (CZ, TwoQubitBatch(pairs)) => self.tab.cz_batch(pairs),

            // TODO: replace things below by actual batched methods once they are available
            // Batch: single-qubit for loops
            (X, QubitBatch(addrs)) => batch_for!(self.tab, x, addrs),
            (Y, QubitBatch(addrs)) => batch_for!(self.tab, y, addrs),
            (Z, QubitBatch(addrs)) => batch_for!(self.tab, z, addrs),
            (S, QubitBatch(addrs)) => batch_for!(self.tab, s, addrs),
            (SAdj, QubitBatch(addrs)) => batch_for!(self.tab, s_adj, addrs),
            (T, QubitBatch(addrs)) => batch_for!(self.tab, t, addrs),
            (TAdj, QubitBatch(addrs)) => batch_for!(self.tab, t_adj, addrs),
            (Reset, QubitBatch(addrs)) => batch_for!(self.tab, reset, addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, rx, addrs, *angle),
            (RY, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, ry, addrs, *angle),
            (RZ, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.tab, rz, addrs, *angle),
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.tab, depolarize, addrs, *p)
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
                for &(a, b) in pairs {
                    self.tab.cnot(a, b);
                }
            }
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.rxx(a, b, *angle);
                }
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.ryy(a, b, *angle);
                }
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                for &(a, b) in pairs {
                    self.tab.rzz(a, b, *angle);
                }
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
                let outcomes = addrs.iter().map(|&addr| self.tab.measure(addr).into());
                return Ok(Effects::one(CircuitOutcomeEffect::Measurement(
                    MeasurementEffect {
                        measurement_results: outcomes.collect(),
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
                    "Invalid gate arguments {:?} for gate {:?}",
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

/// PauliSum-backed executor (Heisenberg picture). Holds a `PauliSum<T>` and
/// answers the same `CircuitInstruction` vocabulary as `CircuitExecutor`, but
/// without measurement / reset support.
///
/// Skeleton only — `execute_instruction` is a no-op until Task 5 fills in the
/// gate-dispatch table per the plan's Gate Support Matrix. Not yet wired into
/// the `Circuit` enum (that happens in Task 4).
pub struct PauliSumExecutor<T: Config<Coeff = f64>> {
    pub state: PauliSum<T>,
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
        use CircuitMessage::*;

        match (inst, msg) {
            // Single-qubit Clifford
            (X, &Qubit(addr)) => self.state.x(addr),
            (Y, &Qubit(addr)) => self.state.y(addr),
            (Z, &Qubit(addr)) => self.state.z(addr),
            (H, &Qubit(addr)) => self.state.h(addr),
            (S, &Qubit(addr)) => self.state.s(addr),
            (SAdj, &Qubit(addr)) => self.state.s_adj(addr),
            (SqrtX, &Qubit(addr)) => self.state.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => self.state.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => self.state.sqrt_x_adj(addr),
            (SqrtYAdj, &Qubit(addr)) => self.state.sqrt_y_adj(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => self.state.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => self.state.cz(addr0, addr1),

            // Single-qubit rotations
            (RX, &QubitAndFloat(addr, angle)) => self.state.rx(addr, angle),
            (RY, &QubitAndFloat(addr, angle)) => self.state.ry(addr, angle),
            (RZ, &QubitAndFloat(addr, angle)) => self.state.rz(addr, angle),

            // Two-qubit rotations
            (RXX, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.rxx(addr0, addr1, angle),
            (RYY, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.ryy(addr0, addr1, angle),
            (RZZ, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.rzz(addr0, addr1, angle),

            // RXY: rotation about an axis in the x/y plane
            (R, &QubitAndTwoFloats(addr, axis_angle, theta)) => {
                self.state.r(addr, axis_angle, theta)
            }

            // Noise
            (Depolarize, &QubitAndFloat(addr, p)) => self.state.depolarize(addr, p),
            (Depolarize2, &TwoQubitAndFloat(addr0, addr1, p)) => {
                self.state.depolarize2(addr0, addr1, p)
            }
            (PauliError, QubitAndFloatArr3(addr0, ps)) => self.state.pauli_error(*addr0, *ps),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(addr0, addr1, ps)) => {
                self.state.two_qubit_pauli_error(*addr0, *addr1, *ps)
            }

            // Truncate: pruning per the configured strategy.
            (Truncate, None) => self.state.truncate(),

            // Batched arms: simple for-loop dispatch (no dedicated batch
            // methods on PauliSum<T>, unlike GeneralizedTableau).
            (X, QubitBatch(addrs)) => batch_for!(self.state, x, addrs),
            (Y, QubitBatch(addrs)) => batch_for!(self.state, y, addrs),
            (Z, QubitBatch(addrs)) => batch_for!(self.state, z, addrs),
            (H, QubitBatch(addrs)) => batch_for!(self.state, h, addrs),
            (S, QubitBatch(addrs)) => batch_for!(self.state, s, addrs),
            (SAdj, QubitBatch(addrs)) => batch_for!(self.state, s_adj, addrs),
            (SqrtX, QubitBatch(addrs)) => batch_for!(self.state, sqrt_x, addrs),
            (SqrtY, QubitBatch(addrs)) => batch_for!(self.state, sqrt_y, addrs),
            (SqrtXAdj, QubitBatch(addrs)) => batch_for!(self.state, sqrt_x_adj, addrs),
            (SqrtYAdj, QubitBatch(addrs)) => batch_for!(self.state, sqrt_y_adj, addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, rx, addrs, *angle),
            (RY, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, ry, addrs, *angle),
            (RZ, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, rz, addrs, *angle),
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.state, depolarize, addrs, *p)
            }
            (PauliError, QubitBatchAndFloatArr3(addrs, ps)) => {
                batch_for!(self.state, pauli_error, addrs, *ps)
            }
            (CNOT, TwoQubitBatch(pairs)) => batch_pairs_for!(self.state, cnot, pairs),
            (CZ, TwoQubitBatch(pairs)) => batch_pairs_for!(self.state, cz, pairs),
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, rxx, pairs, *angle)
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, ryy, pairs, *angle)
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, rzz, pairs, *angle)
            }
            (Depolarize2, TwoQubitBatchAndFloat(pairs, p)) => {
                batch_pairs_for!(self.state, depolarize2, pairs, *p)
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(pairs, ps)) => {
                batch_pairs_for!(self.state, two_qubit_pauli_error, pairs, *ps)
            }

            // Not supported on PauliSum (Decision 11 + Gate Support Matrix).
            (Measure | Reset, _) => {
                return Err(eyre!("{inst} is not supported on the PauliSum backend"));
            }
            (Loss | CorrelatedLoss, _) => {
                return Err(eyre!(
                    "{inst} is not supported on the PauliSum backend; use the LossyPauliSum backend instead"
                ));
            }

            // T / T_adj / U3 are listed as supported on PauliSum in the plan's
            // Gate Support Matrix, but ppvm-runtime does not yet implement
            // TGate or U3Gate for PauliSum<T> (only for GeneralizedTableau).
            // Flag this finding here; lifting the upstream impls is out of
            // scope for Task 5.
            (T | TAdj | U3, _) => {
                return Err(eyre!(
                    "{inst} on PauliSum requires upstream ppvm-runtime support that is not yet implemented"
                ));
            }

            // Trace: parse the resolved pattern string and compute the trace.
            // Per plan Decision 9, parsing happens on every execution; no
            // module-load caching.
            (Trace, PauliPatternStr(s)) => {
                let pat = PauliPattern::parse(s)
                    .map_err(|e| eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                let value = self.state.trace(&pat);
                return Ok(Effects::one(CircuitOutcomeEffect::Trace(TraceEffect {
                    value,
                })));
            }

            // Fallback (batched messages, mismatched shapes, etc.)
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid gate arguments {:?} for gate {:?} on the PauliSum backend",
                    msg,
                    inst
                ));
            }
        };

        Ok(Effects::None)
    }
}

impl<T> vihaco::Reset for PauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
{
    fn reset(&mut self) {
        // TODO(Task 5/6): rebuild self.state from the seeded observable.
    }
}

/// LossyPauliSum-backed executor. Same shape as `PauliSumExecutor`; the
/// distinction lives at the dispatch level (this executor accepts `Loss` /
/// `CorrelatedLoss`) and at the concrete `T` used by the enclosing
/// `Circuit::LossyPauliSum` variant (a Config whose `PauliWordType` is
/// `LossyPauliWord`, picked in Task 4).
///
/// Skeleton only — Task 5 fills in the dispatch.
pub struct LossyPauliSumExecutor<T: Config<Coeff = f64>> {
    pub state: PauliSum<T>,
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

        match (inst, msg) {
            // Single-qubit Clifford
            (X, &Qubit(addr)) => self.state.x(addr),
            (Y, &Qubit(addr)) => self.state.y(addr),
            (Z, &Qubit(addr)) => self.state.z(addr),
            (H, &Qubit(addr)) => self.state.h(addr),
            (S, &Qubit(addr)) => self.state.s(addr),
            (SAdj, &Qubit(addr)) => self.state.s_adj(addr),
            (SqrtX, &Qubit(addr)) => self.state.sqrt_x(addr),
            (SqrtY, &Qubit(addr)) => self.state.sqrt_y(addr),
            (SqrtXAdj, &Qubit(addr)) => self.state.sqrt_x_adj(addr),
            (SqrtYAdj, &Qubit(addr)) => self.state.sqrt_y_adj(addr),

            // Controlled gates
            (CNOT, &TwoQubit(addr0, addr1)) => self.state.cnot(addr0, addr1),
            (CZ, &TwoQubit(addr0, addr1)) => self.state.cz(addr0, addr1),

            // Single-qubit rotations
            (RX, &QubitAndFloat(addr, angle)) => self.state.rx(addr, angle),
            (RY, &QubitAndFloat(addr, angle)) => self.state.ry(addr, angle),
            (RZ, &QubitAndFloat(addr, angle)) => self.state.rz(addr, angle),

            // Two-qubit rotations
            (RXX, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.rxx(addr0, addr1, angle),
            (RYY, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.ryy(addr0, addr1, angle),
            (RZZ, &TwoQubitAndFloat(addr0, addr1, angle)) => self.state.rzz(addr0, addr1, angle),

            // RXY: rotation about an axis in the x/y plane
            (R, &QubitAndTwoFloats(addr, axis_angle, theta)) => {
                self.state.r(addr, axis_angle, theta)
            }

            // Noise
            (Depolarize, &QubitAndFloat(addr, p)) => self.state.depolarize(addr, p),
            (Depolarize2, &TwoQubitAndFloat(addr0, addr1, p)) => {
                self.state.depolarize2(addr0, addr1, p)
            }
            (PauliError, QubitAndFloatArr3(addr0, ps)) => self.state.pauli_error(*addr0, *ps),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(addr0, addr1, ps)) => {
                self.state.two_qubit_pauli_error(*addr0, *addr1, *ps)
            }

            // Loss (accepted on LossyPauliSum; rejected on plain PauliSum)
            (Loss, &QubitAndFloat(addr, p)) => self.state.loss_channel(addr, p),
            (CorrelatedLoss, TwoQubitAndFloatArr3(addr0, addr1, ps)) => {
                self.state.correlated_loss_channel(*addr0, *addr1, *ps)
            }

            // Truncate: pruning per the configured strategy.
            (Truncate, None) => self.state.truncate(),

            // Batched arms: simple for-loop dispatch (no dedicated batch
            // methods on PauliSum<T>).
            (X, QubitBatch(addrs)) => batch_for!(self.state, x, addrs),
            (Y, QubitBatch(addrs)) => batch_for!(self.state, y, addrs),
            (Z, QubitBatch(addrs)) => batch_for!(self.state, z, addrs),
            (H, QubitBatch(addrs)) => batch_for!(self.state, h, addrs),
            (S, QubitBatch(addrs)) => batch_for!(self.state, s, addrs),
            (SAdj, QubitBatch(addrs)) => batch_for!(self.state, s_adj, addrs),
            (SqrtX, QubitBatch(addrs)) => batch_for!(self.state, sqrt_x, addrs),
            (SqrtY, QubitBatch(addrs)) => batch_for!(self.state, sqrt_y, addrs),
            (SqrtXAdj, QubitBatch(addrs)) => batch_for!(self.state, sqrt_x_adj, addrs),
            (SqrtYAdj, QubitBatch(addrs)) => batch_for!(self.state, sqrt_y_adj, addrs),
            (RX, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, rx, addrs, *angle),
            (RY, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, ry, addrs, *angle),
            (RZ, QubitBatchAndFloat(addrs, angle)) => batch_for!(self.state, rz, addrs, *angle),
            (Depolarize, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.state, depolarize, addrs, *p)
            }
            (Loss, QubitBatchAndFloat(addrs, p)) => {
                batch_for!(self.state, loss_channel, addrs, *p)
            }
            (PauliError, QubitBatchAndFloatArr3(addrs, ps)) => {
                batch_for!(self.state, pauli_error, addrs, *ps)
            }
            (CNOT, TwoQubitBatch(pairs)) => batch_pairs_for!(self.state, cnot, pairs),
            (CZ, TwoQubitBatch(pairs)) => batch_pairs_for!(self.state, cz, pairs),
            (RXX, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, rxx, pairs, *angle)
            }
            (RYY, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, ryy, pairs, *angle)
            }
            (RZZ, TwoQubitBatchAndFloat(pairs, angle)) => {
                batch_pairs_for!(self.state, rzz, pairs, *angle)
            }
            (Depolarize2, TwoQubitBatchAndFloat(pairs, p)) => {
                batch_pairs_for!(self.state, depolarize2, pairs, *p)
            }
            (CorrelatedLoss, TwoQubitBatchAndFloatArr3(pairs, ps)) => {
                batch_pairs_for!(self.state, correlated_loss_channel, pairs, *ps)
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(pairs, ps)) => {
                batch_pairs_for!(self.state, two_qubit_pauli_error, pairs, *ps)
            }

            // Not supported on LossyPauliSum (Decision 11 + Gate Support Matrix).
            (Measure | Reset, _) => {
                return Err(eyre!(
                    "{inst} is not supported on the LossyPauliSum backend"
                ));
            }

            // See PauliSumExecutor: T/T_adj/U3 require upstream ppvm-runtime
            // impls that don't exist yet.
            (T | TAdj | U3, _) => {
                return Err(eyre!(
                    "{inst} on LossyPauliSum requires upstream ppvm-runtime support that is not yet implemented"
                ));
            }

            // Trace: parse the resolved pattern string and compute the trace.
            // Per plan Decision 9, parsing happens on every execution; no
            // module-load caching.
            (Trace, PauliPatternStr(s)) => {
                let pat = PauliPattern::parse(s)
                    .map_err(|e| eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                let value = self.state.trace(&pat);
                return Ok(Effects::one(CircuitOutcomeEffect::Trace(TraceEffect {
                    value,
                })));
            }

            // Fallback (batched messages, mismatched shapes, etc.)
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid gate arguments {:?} for gate {:?} on the LossyPauliSum backend",
                    msg,
                    inst
                ));
            }
        };

        Ok(Effects::None)
    }
}

impl<T> vihaco::Reset for LossyPauliSumExecutor<T>
where
    T: Config<Coeff = f64>,
{
    fn reset(&mut self) {
        // TODO(Task 5/6): rebuild self.state from the seeded observable.
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
    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Self {
        if n_qubits <= 64 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits64(CircuitExecutor { tab })
        } else if n_qubits <= 128 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits128(CircuitExecutor { tab })
        } else if n_qubits <= 256 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits256(CircuitExecutor { tab })
        } else if n_qubits <= 512 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits512(CircuitExecutor { tab })
        } else if n_qubits <= 1024 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits1024(CircuitExecutor { tab })
        } else if n_qubits <= 2048 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Self::Bits2048(CircuitExecutor { tab })
        } else {
            panic!("No matching executor for {} qubits", n_qubits);
        }
    }

    /// Same as [`TableauCircuit::new`], but seed the RNG deterministically so a
    /// shot is reproducible.
    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: f64, seed: u64) -> Self {
        macro_rules! seeded {
            ($variant:ident) => {{
                let tab = GeneralizedTableau::new_with_seed(n_qubits, coefficient_threshold, seed);
                Self::$variant(CircuitExecutor { tab })
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
        } else if n_qubits <= 2048 {
            seeded!(Bits2048)
        } else {
            panic!("No matching executor for {} qubits", n_qubits);
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
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Self {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<PauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                Self::$variant(PauliSumExecutor { state })
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
        } else if info.n_qubits <= 2048 {
            build!(Bits2048, 256)
        } else {
            panic!("No matching PauliSum executor for {} qubits", info.n_qubits);
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
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Self {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<LossyPauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                Self::$variant(LossyPauliSumExecutor { state })
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
        } else if info.n_qubits <= 2048 {
            build!(Bits2048, 256)
        } else {
            panic!(
                "No matching LossyPauliSum executor for {} qubits",
                info.n_qubits
            );
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
    pub fn tableau(info: &PPVMDeviceInfo) -> Self {
        Self::Tableau(TableauCircuit::new(
            info.n_qubits,
            info.coefficient_threshold,
        ))
    }

    /// Same as [`Circuit::tableau`], but seed the tableau RNG deterministically
    /// so a shot is reproducible.
    pub fn tableau_with_seed(info: &PPVMDeviceInfo, seed: u64) -> Self {
        Self::Tableau(TableauCircuit::new_with_seed(
            info.n_qubits,
            info.coefficient_threshold,
            seed,
        ))
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term in
    /// `terms`. Each `(word, coef)` is added via `state += (word, coef)`; the
    /// caller is responsible for having parsed/validated the words against
    /// `info.n_qubits` (see `parse_observable_terms` in `composite.rs`).
    pub fn paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Self {
        Self::PauliSum(PauliSumCircuit::new(info, terms))
    }

    /// Build a LossyPauliSum-backed circuit. Same contract as
    /// [`Circuit::paulisum`].
    pub fn lossy_paulisum(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Self {
        Self::LossyPauliSum(LossyPauliSumCircuit::new(info, terms))
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
        // Default backend is Tableau, which doesn't require an observable.
        Self::tableau(&PPVMDeviceInfo::default())
    }
}
