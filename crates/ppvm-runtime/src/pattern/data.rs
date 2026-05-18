// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

// pattern intermediate representation

#[cfg(feature = "bincode")]
use bincode::{Decode, Encode};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
#[repr(u8)]
pub(crate) enum NotIdentity {
    X = 1,
    Z = 2,
    Y = 3,
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
pub(crate) enum OpPattern {
    Identity,                                   // I
    Single(NotIdentity),                        // X, Y, Z
    Double(NotIdentity, NotIdentity),           // [XY]
    AnyNonIdentity,                             // [XYZ]
    SingleOrIdentity(NotIdentity),              // X?
    DoubleOrIdentity(NotIdentity, NotIdentity), // [XY]?
    AnyPauliOrIdentity,                         // [XYZ]?
}

impl OpPattern {
    pub fn add_identity(self) -> Self {
        match self {
            OpPattern::Single(not_identity) => OpPattern::SingleOrIdentity(not_identity),
            OpPattern::Double(left, right) => OpPattern::DoubleOrIdentity(left, right),
            OpPattern::AnyNonIdentity => OpPattern::AnyPauliOrIdentity,
            _ => self,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
pub(crate) enum Decorated {
    Position(OpPattern, usize), // X1, Y2, Z3
    Star(OpPattern),            // [XY]*
    Repeat(OpPattern, usize),   // <OpPattern>{2}
}

/// A Pauli pattern representing a set of matching Pauli words.
///
/// # Example
/// `[XY]2Z3` represents `X2Z3` or `Y2Z3`.
/// `Z?` represents any Pauli word with `I` or `Z`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
pub struct PauliPattern(pub(super) Vec<Decorated>);
