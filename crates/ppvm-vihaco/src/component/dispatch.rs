// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// [`draw!`](super::backend::draw) over a batch of single-qubit targets.
macro_rules! batch_draw {
    ($ex:expr, $field:ident, $method:ident, $items:expr $(, $arg:expr)*) => {
        for item in $items {
            draw!($ex, $field, $method(*item $(, $arg)*));
        }
    };
}

/// [`draw!`](super::backend::draw) over a batch of qubit pairs.
macro_rules! batch_pairs_draw {
    ($ex:expr, $field:ident, $method:ident, $items:expr $(, $arg:expr)*) => {
        for &(a, b) in $items {
            draw!($ex, $field, $method(a, b $(, $arg)*));
        }
    };
}

macro_rules! dispatch_pauli_sum {
    ($ex:expr, $state:ident, $inst:expr, $msg:expr, $name:literal) => {{
        use crate::component::backend::{
            Clifford as _, CliffordBatch as _, CliffordExtensions as _,
            CliffordExtensionsBatch as _, Depolarizing as _, Depolarizing2 as _, PauliError as _,
            PauliPattern, RotXY as _, RotationOne as _, RotationTwo as _, Trace as _,
            TwoQubitPauliError as _,
        };
        use vihaco_circuit_isa::{CircuitInstruction::*, CircuitMessage::*};

        match ($inst, $msg) {
            (X, &Qubit(q)) => $ex.$state.x(q),
            (Y, &Qubit(q)) => $ex.$state.y(q),
            (Z, &Qubit(q)) => $ex.$state.z(q),
            (H, &Qubit(q)) => $ex.$state.h(q),
            (S, &Qubit(q)) => $ex.$state.s(q),
            (SAdj, &Qubit(q)) => $ex.$state.s_dag(q),
            (SqrtX, &Qubit(q)) => $ex.$state.sqrt_x(q),
            (SqrtY, &Qubit(q)) => $ex.$state.sqrt_y(q),
            (SqrtXAdj, &Qubit(q)) => $ex.$state.sqrt_x_dag(q),
            (SqrtYAdj, &Qubit(q)) => $ex.$state.sqrt_y_dag(q),
            (CNOT, &TwoQubit(a, b)) => $ex.$state.cnot(a, b),
            (CZ, &TwoQubit(a, b)) => $ex.$state.cz(a, b),
            (RX, &QubitAndFloat(q, a)) => $ex.$state.rx(q, a),
            (RY, &QubitAndFloat(q, a)) => $ex.$state.ry(q, a),
            (RZ, &QubitAndFloat(q, a)) => $ex.$state.rz(q, a),
            (RXX, &TwoQubitAndFloat(a, b, t)) => $ex.$state.rxx(a, b, t),
            (RYY, &TwoQubitAndFloat(a, b, t)) => $ex.$state.ryy(a, b, t),
            (RZZ, &TwoQubitAndFloat(a, b, t)) => $ex.$state.rzz(a, b, t),
            (R, &QubitAndTwoFloats(q, a, t)) => $ex.$state.r(q, a, t),
            (Depolarize, &QubitAndFloat(q, p)) => draw!($ex, $state, depolarize1(q, p)),
            (Depolarize2, &TwoQubitAndFloat(a, b, p)) => draw!($ex, $state, depolarize2(a, b, p)),
            (PauliError, QubitAndFloatArr3(q, p)) => draw!($ex, $state, pauli_error(*q, *p)),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(a, b, p)) => {
                draw!($ex, $state, two_qubit_pauli_error(*a, *b, *p))
            }
            (Truncate, None) => $ex.$state.truncate(),
            (X, QubitBatch(qs)) => $ex.$state.x_many(qs),
            (Y, QubitBatch(qs)) => $ex.$state.y_many(qs),
            (Z, QubitBatch(qs)) => $ex.$state.z_many(qs),
            (H, QubitBatch(qs)) => $ex.$state.h_many(qs),
            (S, QubitBatch(qs)) => $ex.$state.s_many(qs),
            (SAdj, QubitBatch(qs)) => $ex.$state.s_dag_many(qs),
            (SqrtX, QubitBatch(qs)) => $ex.$state.sqrt_x_many(qs),
            (SqrtY, QubitBatch(qs)) => $ex.$state.sqrt_y_many(qs),
            (SqrtXAdj, QubitBatch(qs)) => $ex.$state.sqrt_x_dag_many(qs),
            (SqrtYAdj, QubitBatch(qs)) => $ex.$state.sqrt_y_dag_many(qs),
            (RX, QubitBatchAndFloat(qs, a)) => $ex.$state.rx_many(qs, *a),
            (RY, QubitBatchAndFloat(qs, a)) => $ex.$state.ry_many(qs, *a),
            (RZ, QubitBatchAndFloat(qs, a)) => $ex.$state.rz_many(qs, *a),
            (Depolarize, QubitBatchAndFloat(qs, p)) => {
                batch_draw!($ex, $state, depolarize1, qs, *p)
            }
            (PauliError, QubitBatchAndFloatArr3(qs, p)) => {
                batch_draw!($ex, $state, pauli_error, qs, *p)
            }
            (CNOT, TwoQubitBatch(ps)) => $ex.$state.cnot_many(ps),
            (CZ, TwoQubitBatch(ps)) => $ex.$state.cz_many(ps),
            (RXX, TwoQubitBatchAndFloat(ps, a)) => $ex.$state.rxx_many(ps, *a),
            (RYY, TwoQubitBatchAndFloat(ps, a)) => $ex.$state.ryy_many(ps, *a),
            (RZZ, TwoQubitBatchAndFloat(ps, a)) => $ex.$state.rzz_many(ps, *a),
            (Depolarize2, TwoQubitBatchAndFloat(ps, p)) => {
                batch_pairs_draw!($ex, $state, depolarize2, ps, *p)
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(ps, p)) => {
                batch_pairs_draw!($ex, $state, two_qubit_pauli_error, ps, *p)
            }
            (Measure | Reset, _) => {
                return Err(eyre::eyre!(
                    "{} is not supported on the {} backend",
                    $inst,
                    $name
                ));
            }
            (T, Qubit(q)) => $ex.$state.rz(*q, std::f64::consts::PI / 8.0),
            (TAdj, Qubit(q)) => $ex.$state.rz(*q, -std::f64::consts::PI / 8.0),
            (T, QubitBatch(qs)) => $ex.$state.rz_many(qs, std::f64::consts::PI / 8.0),
            (TAdj, QubitBatch(qs)) => $ex.$state.rz_many(qs, -std::f64::consts::PI / 8.0),
            (U3, &QubitU3(q, theta, phi, lam)) => {
                $ex.$state.rz(q, phi);
                $ex.$state.ry(q, theta);
                $ex.$state.rz(q, lam);
            }
            (U3, QubitBatchU3(qs, theta, phi, lam)) => {
                $ex.$state.rz_many(qs, *phi);
                $ex.$state.ry_many(qs, *theta);
                $ex.$state.rz_many(qs, *lam);
            }
            (Trace, PauliPatternStr(s)) => {
                let pattern = PauliPattern::parse(s)
                    .map_err(|e| eyre::eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                let value = $ex.$state.trace(&pattern);
                return Ok(vihaco::Effects::one(
                    crate::measurements::CircuitOutcomeEffect::Trace(
                        crate::measurements::TraceEffect { value },
                    ),
                ));
            }
            (inst, msg) => {
                return Err(eyre::eyre!(
                    "Invalid circuit instruction arguments {:?} for instruction {:?} on the {} backend",
                    msg,
                    inst,
                    $name
                ));
            }
        }
        Ok(vihaco::Effects::None)
    }};
}

pub(crate) use {batch_draw, batch_pairs_draw, dispatch_pauli_sum};
