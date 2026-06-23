// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::data::{NotIdentity, PauliPattern};
use ppvm_traits::char::Pauli;

impl<S: AsRef<str>> From<S> for PauliPattern {
    fn from(s: S) -> Self {
        PauliPattern::parse(s.as_ref()).expect("Failed to parse Pauli pattern")
    }
}

impl From<NotIdentity> for Pauli {
    fn from(op: NotIdentity) -> Self {
        unsafe { std::mem::transmute(op as u8) }
    }
}
