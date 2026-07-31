// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! `ppvm-traits-2`: the trait foundation for the second trait-system experiment.
//!
//! This crate holds only trait *definitions* and small leaf types — the split
//! `Coefficient`/`Angle`, `Word`/`Indexable`/`PauliBits`, the Pauli algebra
//! traits, the gate/noise traits, the graded map layers, the algebra
//! capabilities (`KeyProduct`/`ImaginaryUnit`/`Conjugate`), and the batch
//! contract. See [`docs/design/traits-2-configuration-and-hashing.md`] and
//! [`docs/design/traits-2-implementation-plan.md`].
//!
//! The concrete `PauliWord`/`Tableau` representations are deliberately *not*
//! here; they land in the `ppvm-pauli-word-2` / `ppvm-tableau-2` crates
//! (implementation-plan Phases 2 and 4). The only executable code in this crate
//! is a handful of leaf types and the blanket
//! `impl<T: SymplecticColumns + PhaseTrack> Clifford for T`.

pub mod algebra;
pub mod batch;
pub mod coefficient;
pub mod gates;
pub mod graded;
pub mod hash;
pub mod pauli;
pub mod word;

pub use algebra::{Conjugate, ImaginaryUnit, KeyProduct, Phase};
pub use batch::{Columnar, KeyBatch, KeyColumn, TermBatch, TermProducer, TermSink};
pub use coefficient::{Angle, Coefficient, Halvable};
pub use gates::{Clifford, Measure, PauliError, RotationOne};
pub use graded::{Accumulate, Multiply, Pair, Retain, Scale, Support};
pub use hash::{IdentityBuildHasher, IdentityHasher, Indexable};
pub use pauli::{BlanketClifford, PhaseTrack, StabilizerFrame, SymplecticColumns};
pub use word::{FermionAction, FermionSite, LossySite, Pauli, PauliBits, Word};

/// Convenient re-export of the whole trait surface.
pub mod prelude {
    pub use crate::algebra::{Conjugate, ImaginaryUnit, KeyProduct, Phase};
    pub use crate::batch::{Columnar, KeyBatch, KeyColumn, TermBatch, TermProducer, TermSink};
    pub use crate::coefficient::{Angle, Coefficient, Halvable};
    pub use crate::gates::{Clifford, Measure, PauliError, RotationOne};
    pub use crate::graded::{Accumulate, Multiply, Pair, Retain, Scale, Support};
    pub use crate::hash::{IdentityBuildHasher, IdentityHasher, Indexable};
    pub use crate::pauli::{BlanketClifford, PhaseTrack, StabilizerFrame, SymplecticColumns};
    pub use crate::word::{FermionAction, FermionSite, LossySite, Pauli, PauliBits, Word};
}
