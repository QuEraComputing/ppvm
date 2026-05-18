// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

/// A single-qubit Pauli symbol.
///
/// The four standard Paulis are encoded so that the low two bits identify
/// the operator (`I = 0b00`, `X = 0b01`, `Z = 0b10`, `Y = 0b11`), which is
/// the same encoding stabilizer-formalism tools commonly use. The extra
/// variant [`Pauli::L`] marks a qubit as *lost*, used by the loss-aware
/// portions of the runtime.
///
/// # Examples
///
/// ```
/// use ppvm_runtime::char::Pauli;
///
/// assert_eq!(Pauli::I.to_string(), "I");
/// assert_eq!(Pauli::X.to_string(), "X");
/// assert_eq!(Pauli::Y.to_string(), "Y");
/// assert_eq!(Pauli::Z.to_string(), "Z");
/// assert_eq!(Pauli::L.to_string(), "L");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Pauli {
    /// Identity `I`.
    I = 0, // 0b00
    /// Pauli `X`.
    X = 1, // 0b01
    /// Pauli `Z`.
    Z = 2, // 0b10
    /// Pauli `Y`.
    Y = 3, // 0b11
    /// Loss marker — the qubit has been lost and no longer participates.
    L = 4, // 0b100
}

impl std::fmt::Display for Pauli {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pauli::I => write!(f, "I"),
            Pauli::X => write!(f, "X"),
            Pauli::Y => write!(f, "Y"),
            Pauli::Z => write!(f, "Z"),
            Pauli::L => write!(f, "L"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paul_word_display() {
        insta::assert_yaml_snapshot!(Pauli::I.to_string());
        insta::assert_yaml_snapshot!(Pauli::X.to_string());
        insta::assert_yaml_snapshot!(Pauli::Y.to_string());
        insta::assert_yaml_snapshot!(Pauli::Z.to_string());
    }
}
