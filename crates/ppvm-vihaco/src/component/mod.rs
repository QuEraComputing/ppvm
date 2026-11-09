// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

mod backend;
mod buckets;
mod executors;
mod ops;

pub use backend::{CircuitBackend, MAX_QUBITS, PauliSumBackend};
pub use executors::{Circuit, CircuitExecutor, TableauCircuit};
pub use ops::CircuitOp;
pub use vihaco_circuit_isa::messages::{
    PopParams, QubitBatch, QubitPairBatch, TraceMessage, TwoQubitPauliErrorMsg, pop_qubit_batch,
    pop_qubit_pair_batch,
};

// vihaco resolves `#[device]` children through these generated companions.
pub use executors::*;
