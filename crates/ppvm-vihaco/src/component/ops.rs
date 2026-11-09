// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::backend::CircuitBackend;
use crate::measurements::{CircuitOutcomeEffect, MeasurementEffect};
use ppvm_tableau::prelude::*;
use vihaco::{Effects, PopFromStack};
use vihaco_circuit_isa::circuit;

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
