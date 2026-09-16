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
    Multiply, Pair, Retain, Scale, Support, TermBatch, TermSink,
};
pub use gates::{
    AmplitudeDamping, AsymmetricLossChannel, CRx, Clifford, CliffordBatch, CliffordExtensions,
    CliffordExtensionsBatch, CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel,
    Measure, PauliError, PauliErrorAll, PauliErrorFactors, Projection, Reset, ResetLossChannel,
    RotXY, RotationOne, RotationOneBatch, RotationTwo, RotationTwoBatch, TGate, TwoQubitPauliError,
    U3Gate,
};
pub use loss::LossState;
pub use pauli::{BlanketClifford, Pauli, PhaseTrack, SymplecticColumns};
pub use word::{PauliBits, Word};

/// Common traits and types for implementing and using the interfaces.
pub mod prelude {
    pub use crate::{
        Accumulate, AmplitudeDamping, Angle, AsymmetricLossChannel, BlanketClifford, CRx, Clifford,
        CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, Coefficient, Columnar,
        Conjugate, CorrelatedLossChannel, Depolarizing, Depolarizing2, Halvable,
        IdentityBuildHasher, IdentityHasher, ImaginaryUnit, Indexable, KeyBatch, KeyColumn,
        KeyProduct, LossChannel, LossState, Measure, Multiply, Pair, Pauli, PauliBits, PauliError,
        PauliErrorAll, PauliErrorFactors, Phase, PhaseTrack, Projection, Reset, ResetLossChannel,
        Retain, RotXY, RotationOne, RotationOneBatch, RotationTwo, RotationTwoBatch, Scale,
        Support, SymplecticColumns, TGate, TermBatch, TermSink, TwoQubitPauliError, U3Gate, Word,
    };
}
