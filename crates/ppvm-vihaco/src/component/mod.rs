// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod backend;
mod buckets;
mod executors;
mod messages;
mod ops;

pub use backend::{CircuitBackend, MAX_QUBITS, PauliSumBackend};
pub use executors::{Circuit, CircuitExecutor, TableauCircuit};
pub use messages::{PopParams, QubitBatch, TraceMessage, TwoQubitPauliErrorMsg, pop_qubit_batch};
pub use ops::CircuitOp;

// vihaco resolves `#[device]` children through these generated companions.
pub use executors::*;
