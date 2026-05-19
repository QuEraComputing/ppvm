// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::data::{Decorated, NotIdentity, OpPattern, PauliPattern};

impl std::fmt::Display for NotIdentity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NotIdentity::X => write!(f, "X"),
            NotIdentity::Y => write!(f, "Y"),
            NotIdentity::Z => write!(f, "Z"),
        }
    }
}

impl std::fmt::Display for OpPattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OpPattern::Identity => write!(f, "I"),
            OpPattern::Single(not_identity) => write!(f, "{}", not_identity),
            OpPattern::Double(left, right) => write!(f, "[{}{}]", left, right),
            OpPattern::AnyNonIdentity => write!(f, "[XYZ]"),
            OpPattern::SingleOrIdentity(not_identity) => write!(f, "{}?", not_identity),
            OpPattern::DoubleOrIdentity(left, right) => write!(f, "[{}{}]?", left, right),
            OpPattern::AnyPauliOrIdentity => write!(f, "_"),
        }
    }
}

impl std::fmt::Display for Decorated {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Decorated::Position(op, pos) => write!(f, "{}{}", op, pos),
            Decorated::Star(op) => write!(f, "{}*", op),
            Decorated::Repeat(op, count) => write!(f, "{}{{{}}}", op, count),
        }
    }
}

impl std::fmt::Display for PauliPattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for pattern in &self.0 {
            write!(f, "{}", pattern)?;
        }
        Ok(())
    }
}
