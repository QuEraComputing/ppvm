// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The circuit component of the PPVM: a `Circuit` that dispatches
//! [`vihaco_circuit_isa::CircuitInstruction`]s to a tableau or PauliSum
//! backend. The composite machine that drives it is added in a later PR.

pub mod component;
pub mod device_info;
pub mod measurements;

/// Re-exported so consumers can name gates for the circuit component without
/// depending on the ISA crate directly.
pub use vihaco_circuit_isa::CircuitInstruction;

pub mod prelude {
    pub use crate::component::Circuit;
    pub use crate::device_info::{BackendKind, PPVMDeviceInfo};
}
