// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use eyre::{Result, eyre};
use vihaco::Effects;
use vihaco_circuit_isa::{CircuitInstruction, CircuitMessage};

use super::MAX_QUBITS;
#[cfg(feature = "legacy")]
use super::backend::LossyMeasure as _;
#[cfg(feature = "traits-2")]
use super::backend::Measure as _;
use super::backend::{
    Backend, Clifford as _, CliffordBatch as _, CliffordExtensions as _,
    CliffordExtensionsBatch as _, CorrelatedLossChannel as _, Depolarizing as _,
    Depolarizing2 as _, LossChannel as _, PauliError as _, PauliPattern, Reset as _, RotXY as _,
    RotationOne as _, RotationTwo as _, TGate as _, Tab64, Tab128, Tab256, Tab512, Tab1024,
    Tab2048, TwoQubitPauliError as _, U3Gate as _,
};
use super::dispatch::batch_for;
use crate::measurements::{
    CircuitOutcomeEffect, MeasurementEffect, MeasurementOutcome, TraceEffect,
};

pub struct CircuitExecutor<T> {
    pub tab: T,
}

pub enum TableauCircuit {
    Bits64(CircuitExecutor<Tab64>),
    Bits128(CircuitExecutor<Tab128>),
    Bits256(CircuitExecutor<Tab256>),
    Bits512(CircuitExecutor<Tab512>),
    Bits1024(CircuitExecutor<Tab1024>),
    Bits2048(CircuitExecutor<Tab2048>),
}

macro_rules! dispatch {
    ($tab:expr, $inst:expr, $msg:expr) => {{
        use vihaco_circuit_isa::{CircuitInstruction::*, CircuitMessage::*};
        match ($inst, $msg) {
            (X, &Qubit(q)) => $tab.x(q),
            (Y, &Qubit(q)) => $tab.y(q),
            (Z, &Qubit(q)) => $tab.z(q),
            (H, &Qubit(q)) => $tab.h(q),
            (S, &Qubit(q)) => $tab.s(q),
            (SAdj, &Qubit(q)) => $tab.s_dag(q),
            (SqrtX, &Qubit(q)) => $tab.sqrt_x(q),
            (SqrtY, &Qubit(q)) => $tab.sqrt_y(q),
            (SqrtXAdj, &Qubit(q)) => $tab.sqrt_x_dag(q),
            (SqrtYAdj, &Qubit(q)) => $tab.sqrt_y_dag(q),
            (CNOT, &TwoQubit(a, b)) => $tab.cnot(a, b),
            (CZ, &TwoQubit(a, b)) => $tab.cz(a, b),
            (T, &Qubit(q)) => $tab.t(q),
            (TAdj, &Qubit(q)) => $tab.t_dag(q),
            (RX, &QubitAndFloat(q, a)) => $tab.rx(q, a),
            (RY, &QubitAndFloat(q, a)) => $tab.ry(q, a),
            (RZ, &QubitAndFloat(q, a)) => $tab.rz(q, a),
            (RXX, &TwoQubitAndFloat(a, b, t)) => $tab.rxx(a, b, t),
            (RYY, &TwoQubitAndFloat(a, b, t)) => $tab.ryy(a, b, t),
            (RZZ, &TwoQubitAndFloat(a, b, t)) => $tab.rzz(a, b, t),
            (U3, &QubitU3(q, t, p, l)) => $tab.u3(q, t, p, l),
            (R, &QubitAndTwoFloats(q, a, t)) => $tab.r(q, a, t),
            (Measure, &Qubit(q)) => {
                let outcome: MeasurementOutcome = $tab.measure(q).into();
                return Ok(Effects::one(CircuitOutcomeEffect::Measurement(
                    MeasurementEffect {
                        measurement_results: smallvec::smallvec![outcome],
                    },
                )));
            }
            (Reset, &Qubit(q)) => $tab.reset(q),
            (Depolarize, &QubitAndFloat(q, p)) => $tab.depolarize1(q, p),
            (Depolarize2, &TwoQubitAndFloat(a, b, p)) => $tab.depolarize2(a, b, p),
            (PauliError, QubitAndFloatArr3(q, p)) => $tab.pauli_error(*q, *p),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(a, b, p)) => {
                $tab.two_qubit_pauli_error(*a, *b, *p)
            }
            (Loss, &QubitAndFloat(q, p)) => $tab.loss_channel(q, p),
            (CorrelatedLoss, TwoQubitAndFloatArr3(a, b, p)) => {
                $tab.correlated_loss_channel(*a, *b, *p)
            }
            (X, QubitBatch(qs)) => $tab.x_many(qs),
            (Y, QubitBatch(qs)) => $tab.y_many(qs),
            (Z, QubitBatch(qs)) => $tab.z_many(qs),
            (H, QubitBatch(qs)) => $tab.h_many(qs),
            (S, QubitBatch(qs)) => $tab.s_many(qs),
            (SAdj, QubitBatch(qs)) => $tab.s_dag_many(qs),
            (SqrtX, QubitBatch(qs)) => $tab.sqrt_x_many(qs),
            (SqrtY, QubitBatch(qs)) => $tab.sqrt_y_many(qs),
            (SqrtXAdj, QubitBatch(qs)) => $tab.sqrt_x_dag_many(qs),
            (SqrtYAdj, QubitBatch(qs)) => $tab.sqrt_y_dag_many(qs),
            (T, QubitBatch(qs)) => $tab.t_many(qs),
            (TAdj, QubitBatch(qs)) => $tab.t_dag_many(qs),
            (Reset, QubitBatch(qs)) => $tab.reset_many(qs),
            (RX, QubitBatchAndFloat(qs, a)) => $tab.rx_many(qs, *a),
            (RY, QubitBatchAndFloat(qs, a)) => $tab.ry_many(qs, *a),
            (RZ, QubitBatchAndFloat(qs, a)) => $tab.rz_many(qs, *a),
            (U3, QubitBatchU3(qs, t, p, l)) => {
                for &q in qs {
                    $tab.u3(q, *t, *p, *l);
                }
            }
            (CNOT, TwoQubitBatch(ps)) => $tab.cnot_many(ps),
            (CZ, TwoQubitBatch(ps)) => $tab.cz_many(ps),
            (RXX, TwoQubitBatchAndFloat(ps, a)) => $tab.rxx_many(ps, *a),
            (RYY, TwoQubitBatchAndFloat(ps, a)) => $tab.ryy_many(ps, *a),
            (RZZ, TwoQubitBatchAndFloat(ps, a)) => $tab.rzz_many(ps, *a),
            (Depolarize, QubitBatchAndFloat(qs, p)) => batch_for!($tab, depolarize1, qs, *p),
            (Loss, QubitBatchAndFloat(qs, p)) => batch_for!($tab, loss_channel, qs, *p),
            (PauliError, QubitBatchAndFloatArr3(qs, p)) => {
                batch_for!($tab, pauli_error, qs, *p)
            }
            (Depolarize2, TwoQubitBatchAndFloat(ps, p)) => {
                for &(a, b) in ps {
                    $tab.depolarize2(a, b, *p);
                }
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(ps, p)) => {
                for &(a, b) in ps {
                    $tab.two_qubit_pauli_error(a, b, *p);
                }
            }
            (CorrelatedLoss, TwoQubitBatchAndFloatArr3(ps, p)) => {
                for &(a, b) in ps {
                    $tab.correlated_loss_channel(a, b, *p);
                }
            }
            (Measure, QubitBatch(qs)) => {
                let outcomes = $tab.measure_many(qs);
                return Ok(Effects::one(CircuitOutcomeEffect::Measurement(
                    MeasurementEffect {
                        measurement_results: outcomes.into_iter().map(Into::into).collect(),
                    },
                )));
            }
            (Truncate, None) => {}
            (Trace, PauliPatternStr(s)) => {
                let pattern = PauliPattern::parse(s)
                    .map_err(|e| eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                return Ok(Effects::one(CircuitOutcomeEffect::Trace(TraceEffect {
                    value: $tab.trace(&pattern),
                })));
            }
            (inst, msg) => {
                return Err(eyre!(
                    "Invalid circuit instruction arguments {:?} for instruction {:?}",
                    msg,
                    inst
                ));
            }
        }
        Ok(Effects::None)
    }};
}

impl TableauCircuit {
    pub fn new(n: usize, threshold: f64) -> Result<Self> {
        Self::build(n, threshold, None)
    }

    pub fn new_with_seed(n: usize, threshold: f64, seed: u64) -> Result<Self> {
        Self::build(n, threshold, Some(seed))
    }

    fn build(n: usize, threshold: f64, seed: Option<u64>) -> Result<Self> {
        macro_rules! make {
            ($variant:ident, $new:ident, $seeded:ident) => {
                Ok(Self::$variant(CircuitExecutor {
                    tab: match seed {
                        Some(seed) => Backend::$seeded(n, threshold, seed),
                        None => Backend::$new(n, threshold),
                    },
                }))
            };
        }
        match n {
            0..=64 => make!(Bits64, tab64, tab64_seeded),
            65..=128 => make!(Bits128, tab128, tab128_seeded),
            129..=256 => make!(Bits256, tab256, tab256_seeded),
            257..=512 => make!(Bits512, tab512, tab512_seeded),
            513..=1024 => make!(Bits1024, tab1024, tab1024_seeded),
            1025..=MAX_QUBITS => make!(Bits2048, tab2048, tab2048_seeded),
            _ => Err(eyre!("cannot simulate {n} qubits: maximum is {MAX_QUBITS}")),
        }
    }

    pub(super) fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        match self {
            Self::Bits64(ex) => dispatch!(ex.tab, inst, msg),
            Self::Bits128(ex) => dispatch!(ex.tab, inst, msg),
            Self::Bits256(ex) => dispatch!(ex.tab, inst, msg),
            Self::Bits512(ex) => dispatch!(ex.tab, inst, msg),
            Self::Bits1024(ex) => dispatch!(ex.tab, inst, msg),
            Self::Bits2048(ex) => dispatch!(ex.tab, inst, msg),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => Backend::render(&ex.tab),
            Self::Bits128(ex) => Backend::render(&ex.tab),
            Self::Bits256(ex) => Backend::render(&ex.tab),
            Self::Bits512(ex) => Backend::render(&ex.tab),
            Self::Bits1024(ex) => Backend::render(&ex.tab),
            Self::Bits2048(ex) => Backend::render(&ex.tab),
        }
    }
}

impl vihaco::Reset for TableauCircuit {
    fn reset(&mut self) {
        match self {
            Self::Bits64(ex) => ex.tab.reset_all(),
            Self::Bits128(ex) => ex.tab.reset_all(),
            Self::Bits256(ex) => ex.tab.reset_all(),
            Self::Bits512(ex) => ex.tab.reset_all(),
            Self::Bits1024(ex) => ex.tab.reset_all(),
            Self::Bits2048(ex) => ex.tab.reset_all(),
        }
    }
}
