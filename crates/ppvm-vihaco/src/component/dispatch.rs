// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

macro_rules! batch_for {
    ($state:expr, $method:ident, $items:expr $(, $arg:expr)*) => {
        for item in $items {
            $state.$method(*item $(, $arg)*);
        }
    };
}

macro_rules! batch_pairs_for {
    ($state:expr, $method:ident, $items:expr $(, $arg:expr)*) => {
        for &(a, b) in $items {
            $state.$method(a, b $(, $arg)*);
        }
    };
}

macro_rules! dispatch_pauli_sum {
    ($state:expr, $inst:expr, $msg:expr, $name:literal) => {{
        use crate::component::backend::{
            Clifford as _, CliffordBatch as _, CliffordExtensions as _,
            CliffordExtensionsBatch as _, Depolarizing as _, Depolarizing2 as _, PauliError as _,
            PauliPattern, RotXY as _, RotationOne as _, RotationTwo as _, Trace as _,
            TwoQubitPauliError as _,
        };
        use vihaco_circuit_isa::{CircuitInstruction::*, CircuitMessage::*};

        match ($inst, $msg) {
            (X, &Qubit(q)) => $state.x(q),
            (Y, &Qubit(q)) => $state.y(q),
            (Z, &Qubit(q)) => $state.z(q),
            (H, &Qubit(q)) => $state.h(q),
            (S, &Qubit(q)) => $state.s(q),
            (SAdj, &Qubit(q)) => $state.s_dag(q),
            (SqrtX, &Qubit(q)) => $state.sqrt_x(q),
            (SqrtY, &Qubit(q)) => $state.sqrt_y(q),
            (SqrtXAdj, &Qubit(q)) => $state.sqrt_x_dag(q),
            (SqrtYAdj, &Qubit(q)) => $state.sqrt_y_dag(q),
            (CNOT, &TwoQubit(a, b)) => $state.cnot(a, b),
            (CZ, &TwoQubit(a, b)) => $state.cz(a, b),
            (RX, &QubitAndFloat(q, a)) => $state.rx(q, a),
            (RY, &QubitAndFloat(q, a)) => $state.ry(q, a),
            (RZ, &QubitAndFloat(q, a)) => $state.rz(q, a),
            (RXX, &TwoQubitAndFloat(a, b, t)) => $state.rxx(a, b, t),
            (RYY, &TwoQubitAndFloat(a, b, t)) => $state.ryy(a, b, t),
            (RZZ, &TwoQubitAndFloat(a, b, t)) => $state.rzz(a, b, t),
            (R, &QubitAndTwoFloats(q, a, t)) => $state.r(q, a, t),
            (Depolarize, &QubitAndFloat(q, p)) => $state.depolarize1(q, p),
            (Depolarize2, &TwoQubitAndFloat(a, b, p)) => $state.depolarize2(a, b, p),
            (PauliError, QubitAndFloatArr3(q, p)) => $state.pauli_error(*q, *p),
            (TwoQubitPauliError, TwoQubitAndFloatArr15(a, b, p)) => {
                $state.two_qubit_pauli_error(*a, *b, *p)
            }
            (Truncate, None) => $state.truncate(),
            (X, QubitBatch(qs)) => $state.x_many(qs),
            (Y, QubitBatch(qs)) => $state.y_many(qs),
            (Z, QubitBatch(qs)) => $state.z_many(qs),
            (H, QubitBatch(qs)) => $state.h_many(qs),
            (S, QubitBatch(qs)) => $state.s_many(qs),
            (SAdj, QubitBatch(qs)) => $state.s_dag_many(qs),
            (SqrtX, QubitBatch(qs)) => $state.sqrt_x_many(qs),
            (SqrtY, QubitBatch(qs)) => $state.sqrt_y_many(qs),
            (SqrtXAdj, QubitBatch(qs)) => $state.sqrt_x_dag_many(qs),
            (SqrtYAdj, QubitBatch(qs)) => $state.sqrt_y_dag_many(qs),
            (RX, QubitBatchAndFloat(qs, a)) => $state.rx_many(qs, *a),
            (RY, QubitBatchAndFloat(qs, a)) => $state.ry_many(qs, *a),
            (RZ, QubitBatchAndFloat(qs, a)) => $state.rz_many(qs, *a),
            (Depolarize, QubitBatchAndFloat(qs, p)) => {
                batch_for!($state, depolarize1, qs, *p)
            }
            (PauliError, QubitBatchAndFloatArr3(qs, p)) => {
                batch_for!($state, pauli_error, qs, *p)
            }
            (CNOT, TwoQubitBatch(ps)) => $state.cnot_many(ps),
            (CZ, TwoQubitBatch(ps)) => $state.cz_many(ps),
            (RXX, TwoQubitBatchAndFloat(ps, a)) => $state.rxx_many(ps, *a),
            (RYY, TwoQubitBatchAndFloat(ps, a)) => $state.ryy_many(ps, *a),
            (RZZ, TwoQubitBatchAndFloat(ps, a)) => $state.rzz_many(ps, *a),
            (Depolarize2, TwoQubitBatchAndFloat(ps, p)) => {
                batch_pairs_for!($state, depolarize2, ps, *p)
            }
            (TwoQubitPauliError, TwoQubitBatchAndFloatArr15(ps, p)) => {
                batch_pairs_for!($state, two_qubit_pauli_error, ps, *p)
            }
            (Measure | Reset, _) => {
                return Err(eyre::eyre!(
                    "{} is not supported on the {} backend",
                    $inst,
                    $name
                ));
            }
            (T, Qubit(q)) => $state.rz(*q, std::f64::consts::PI / 8.0),
            (TAdj, Qubit(q)) => $state.rz(*q, -std::f64::consts::PI / 8.0),
            (T, QubitBatch(qs)) => $state.rz_many(qs, std::f64::consts::PI / 8.0),
            (TAdj, QubitBatch(qs)) => $state.rz_many(qs, -std::f64::consts::PI / 8.0),
            (U3, &QubitU3(q, theta, phi, lam)) => {
                $state.rz(q, phi);
                $state.ry(q, theta);
                $state.rz(q, lam);
            }
            (U3, QubitBatchU3(qs, theta, phi, lam)) => {
                $state.rz_many(qs, *phi);
                $state.ry_many(qs, *theta);
                $state.rz_many(qs, *lam);
            }
            (Trace, PauliPatternStr(s)) => {
                let pattern = PauliPattern::parse(s)
                    .map_err(|e| eyre::eyre!("invalid Pauli pattern `{s}`: {e:?}"))?;
                let value = $state.trace(&pattern);
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

pub(crate) use {batch_for, batch_pairs_for, dispatch_pauli_sum};
