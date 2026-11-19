// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Gate, measurement, reset, and noise-channel interfaces.

mod apply_gate_rule;
mod channel;
mod clifford;
mod measure;
mod rot;

pub use apply_gate_rule::ApplyGateRule;
pub use channel::{
    AmplitudeDamping, AsymmetricLossChannel, CorrelatedLossChannel, Depolarizing, Depolarizing2,
    LossChannel, PauliError, ResetLossChannel, TwoQubitPauliError,
};
pub use clifford::{Clifford, CliffordBatch};
pub use measure::{Measure, Projection, Reset, ResetState};
/// The value type of [`ApplyGateRule`], re-exported so containers need not depend on the kernels.
pub use ppvm_pauli_kernels::Bits;
pub use rot::{
    CRx, RotXY, RotationOne, RotationOneBatch, RotationTwo, RotationTwoBatch, TGate, U3Gate,
};
