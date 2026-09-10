// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared arithmetic, word, gate, and container interfaces.
//!
//! Includes small algebraic types, default gate implementations, and
//! container implementations for `Vec` and `HashMap`.

pub mod algebra;
pub mod arithmetic;
pub mod containers;
pub mod fermion_factor;
pub mod gates;
pub mod loss;
pub mod pauli;
pub mod word;

pub use algebra::{Conjugate, ImaginaryUnit, KeyProduct, Phase};
pub use arithmetic::{Angle, Coefficient, Halvable};
pub use containers::{
    Accumulate, Columnar, IdentityBuildHasher, IdentityHasher, Indexable, KeyBatch, KeyColumn,
    Multiply, Pair, RekeyStrategy, Retain, Scale, Support, TermBatch, TermProducer, TermSink,
    Trace,
};
pub use fermion_factor::{FermionAction, FermionSite};
pub use gates::{
    AmplitudeDamping, AsymmetricLossChannel, CRx, Clifford, CliffordBatch, CliffordExtensions,
    CliffordExtensionsBatch, CorrelatedLossChannel, Depolarizing, Depolarizing2, LossChannel,
    Measure, PauliError, PauliErrorAll, PauliErrorFactors, Projection, Reset, ResetLossChannel,
    RotXY, RotationOne, RotationTwo, TGate, TwoQubitPauliError, U3Gate,
};
pub use loss::LossState;
pub use pauli::{BlanketClifford, Pauli, PhaseTrack, StabilizerFrame, SymplecticColumns};
pub use word::{PauliBits, Word};

/// Common traits and types for implementing and using the interfaces.
pub mod prelude {
    pub use crate::{
        Accumulate, AmplitudeDamping, Angle, AsymmetricLossChannel, BlanketClifford, CRx, Clifford,
        CliffordBatch, CliffordExtensions, CliffordExtensionsBatch, Coefficient, Columnar,
        Conjugate, CorrelatedLossChannel, Depolarizing, Depolarizing2, FermionAction, FermionSite,
        Halvable, IdentityBuildHasher, IdentityHasher, ImaginaryUnit, Indexable, KeyBatch,
        KeyColumn, KeyProduct, LossChannel, LossState, Measure, Multiply, Pair, Pauli, PauliBits,
        PauliError, PauliErrorAll, PauliErrorFactors, Phase, PhaseTrack, Projection, RekeyStrategy,
        Reset, ResetLossChannel, Retain, RotXY, RotationOne, RotationTwo, Scale, StabilizerFrame,
        Support, SymplecticColumns, TGate, TermBatch, TermProducer, TermSink, Trace,
        TwoQubitPauliError, U3Gate, Word,
    };
}
