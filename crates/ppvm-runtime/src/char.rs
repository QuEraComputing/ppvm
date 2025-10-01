/// Represents a Pauli operator (I, X, Y, Z)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Pauli {
    I = 0, // 0b00
    X = 1, // 0b01
    Z = 2, // 0b10
    Y = 3, // 0b11
}

impl std::fmt::Display for Pauli {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pauli::I => write!(f, "I"),
            Pauli::X => write!(f, "X"),
            Pauli::Y => write!(f, "Y"),
            Pauli::Z => write!(f, "Z"),
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
