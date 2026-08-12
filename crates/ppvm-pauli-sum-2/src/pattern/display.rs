// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt;

use ppvm_traits_2::Pauli;

use super::data::{Decorated, OpPattern, PauliPattern};

impl fmt::Display for OpPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Identity => f.write_str("I"),
            Self::Single(pauli) => write!(f, "{}", letter(*pauli)),
            Self::Double(a, b) => write!(f, "[{}{}]", letter(*a), letter(*b)),
            Self::AnyNonIdentity => f.write_str("[XYZ]"),
            Self::SingleOrIdentity(pauli) => write!(f, "{}?", letter(*pauli)),
            Self::DoubleOrIdentity(a, b) => {
                write!(f, "[{}{}]?", letter(*a), letter(*b))
            }
            Self::AnyPauliOrIdentity => f.write_str("_"),
        }
    }
}

impl fmt::Display for Decorated {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Position(op, position) => write!(f, "{op}{position}"),
            Self::Star(op) => write!(f, "{op}*"),
            Self::Repeat(op, count) => write!(f, "{op}{{{count}}}"),
        }
    }
}

impl fmt::Display for PauliPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for decorated in &self.0 {
            write!(f, "{decorated}")?;
        }
        Ok(())
    }
}

impl<S: AsRef<str>> From<S> for PauliPattern {
    fn from(value: S) -> Self {
        Self::parse(value).expect("Failed to parse Pauli pattern")
    }
}

fn letter(pauli: Pauli) -> char {
    match pauli {
        Pauli::I => 'I',
        Pauli::X => 'X',
        Pauli::Y => 'Y',
        Pauli::Z => 'Z',
    }
}
