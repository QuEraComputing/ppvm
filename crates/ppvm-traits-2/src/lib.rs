// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared arithmetic, word, gate, and container interfaces.
//! Includes small algebraic types, default gates, and `Vec`/`HashMap` implementations.

pub mod algebra;
pub mod arithmetic;
pub mod containers;
pub mod gates;
pub mod loss;
pub mod pauli;
pub mod word;

pub use algebra::{Conjugate, ImaginaryUnit, KeyProduct, Phase};
pub use arithmetic::{Angle, Coefficient, Halvable};
pub use containers::{
    Accumulate, Columnar, IdentityBuildHasher, IdentityHasher, Indexable, KeyBatch, KeyColumn,
    KeyColumnMut, LossColumn, Multiply, Pair, PauliColumn, Retain, Scale, Support, TermBatch,
    TermSink,
};
pub use gates::{
    AmplitudeDamping, AsymmetricLossChannel, CRx, Clifford, CliffordBatch, CorrelatedLossChannel,
    Depolarizing, Depolarizing2, LossChannel, Measure, PauliError, Projection, Reset,
    ResetLossChannel, RotXY, RotationOne, RotationOneBatch, RotationTwo, RotationTwoBatch, TGate,
    TwoQubitPauliError, U3Gate,
};
pub use loss::LossState;
pub use pauli::{Pauli, PhaseTrack, SymplecticColumns};
pub use word::{PauliBits, Word};

/// Common traits and types for implementing and using the interfaces.
/// Globbed per module so it cannot drift from the re-exports above.
pub mod prelude {
    pub use crate::{
        algebra::*, arithmetic::*, containers::*, gates::*, loss::*, pauli::*, word::*,
    };
}
