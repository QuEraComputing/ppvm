// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parsing of Pauli-observable strings for
//! [`peek_observable_expectation`](crate::data::GeneralizedTableau::peek_observable_expectation).
//!
//! Two notations are accepted, each with an optional leading `+`/`-` sign:
//!
//! * **Sparse product** (Stim `MPP` syntax): factors joined by `*`, e.g.
//!   `"X0*X3*Z5*Y7"`, `"-Z0*Y1"`, `"Z0"`. Each factor is a Pauli letter
//!   followed by a 0-based qubit index.
//! * **Dense**: a string of `I`/`X`/`Y`/`Z` of length `n_qubits`, e.g.
//!   `"IXIZ"`, where position `i` is the Pauli on qubit `i`. `_` is accepted as
//!   an alias for `I` so Stim-style `PauliString` text (e.g. `"+__X_Y_Z"`)
//!   parses directly.
//!
//! The empty string (optionally just a sign) is the identity observable.

use ppvm_traits::char::Pauli;

/// Error returned when an observable string cannot be parsed against a tableau
/// of a given size.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ObservableParseError {
    /// A sparse factor was empty or malformed (e.g. `"X0**Y1"`, `"Q0"`, `"X"`).
    BadToken(String),
    /// A qubit index is `>= n_qubits`.
    QubitOutOfRange { qubit: usize, n_qubits: usize },
    /// The same qubit appears in more than one factor.
    RepeatedQubit(usize),
    /// A dense observable's length does not match `n_qubits`.
    DenseLengthMismatch { got: usize, n_qubits: usize },
}

impl std::fmt::Display for ObservableParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObservableParseError::BadToken(t) => {
                write!(f, "malformed Pauli factor {t:?}; expected e.g. \"X0\"")
            }
            ObservableParseError::QubitOutOfRange { qubit, n_qubits } => {
                write!(
                    f,
                    "qubit {qubit} out of range for a {n_qubits}-qubit system"
                )
            }
            ObservableParseError::RepeatedQubit(q) => {
                write!(f, "qubit {q} appears in more than one factor")
            }
            ObservableParseError::DenseLengthMismatch { got, n_qubits } => write!(
                f,
                "dense observable has length {got}, expected {n_qubits} (one Pauli per qubit)"
            ),
        }
    }
}

impl std::error::Error for ObservableParseError {}

fn pauli_from_char(c: char) -> Option<Pauli> {
    match c.to_ascii_uppercase() {
        // `_` is Stim's identity glyph in dense `PauliString` text (e.g.
        // "+__X_Y_Z"); accept it interchangeably with `I`.
        'I' | '_' => Some(Pauli::I),
        'X' => Some(Pauli::X),
        'Y' => Some(Pauli::Y),
        'Z' => Some(Pauli::Z),
        _ => None,
    }
}

/// Parse an observable string into `(negate, factors)`.
///
/// `negate` reflects a leading `-`. `factors` lists the non-identity Paulis as
/// `(qubit, pauli)` pairs; identity factors are dropped. See the
/// [module docs](self) for the accepted notations.
pub fn parse_observable(
    s: &str,
    n_qubits: usize,
) -> Result<(bool, Vec<(usize, Pauli)>), ObservableParseError> {
    let s = s.trim();
    let (negate, rest) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let rest = rest.trim();

    if rest.is_empty() {
        // Identity observable (`""`, `"+"`, or `"-"`).
        return Ok((negate, Vec::new()));
    }

    // Sparse iff it carries explicit indices (digits) or factor separators.
    let is_sparse = rest.contains('*') || rest.chars().any(|c| c.is_ascii_digit());
    if is_sparse {
        parse_sparse(rest, n_qubits, negate)
    } else {
        parse_dense(rest, n_qubits, negate)
    }
}

fn parse_sparse(
    rest: &str,
    n_qubits: usize,
    negate: bool,
) -> Result<(bool, Vec<(usize, Pauli)>), ObservableParseError> {
    let mut factors: Vec<(usize, Pauli)> = Vec::new();
    let mut seen: Vec<usize> = Vec::new();

    for raw in rest.split('*') {
        let token = raw.trim();
        let mut chars = token.chars();
        let pauli = match chars.next().and_then(pauli_from_char) {
            Some(p) => p,
            None => return Err(ObservableParseError::BadToken(token.to_string())),
        };
        let idx_str: String = chars.collect();
        let qubit: usize = idx_str
            .parse()
            .map_err(|_| ObservableParseError::BadToken(token.to_string()))?;

        if qubit >= n_qubits {
            return Err(ObservableParseError::QubitOutOfRange { qubit, n_qubits });
        }
        if seen.contains(&qubit) {
            return Err(ObservableParseError::RepeatedQubit(qubit));
        }
        seen.push(qubit);
        if pauli != Pauli::I {
            factors.push((qubit, pauli));
        }
    }
    Ok((negate, factors))
}

fn parse_dense(
    rest: &str,
    n_qubits: usize,
    negate: bool,
) -> Result<(bool, Vec<(usize, Pauli)>), ObservableParseError> {
    if rest.chars().count() != n_qubits {
        return Err(ObservableParseError::DenseLengthMismatch {
            got: rest.chars().count(),
            n_qubits,
        });
    }
    let mut factors: Vec<(usize, Pauli)> = Vec::new();
    for (qubit, c) in rest.chars().enumerate() {
        let pauli =
            pauli_from_char(c).ok_or_else(|| ObservableParseError::BadToken(c.to_string()))?;
        if pauli != Pauli::I {
            factors.push((qubit, pauli));
        }
    }
    Ok((negate, factors))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_with_sign_and_indices() {
        assert_eq!(
            parse_observable("-X0*Z2", 3),
            Ok((true, vec![(0, Pauli::X), (2, Pauli::Z)]))
        );
        assert_eq!(parse_observable("Y7", 8), Ok((false, vec![(7, Pauli::Y)])));
        assert_eq!(
            parse_observable("+Z0*Y1", 2),
            Ok((false, vec![(0, Pauli::Z), (1, Pauli::Y)]))
        );
    }

    #[test]
    fn dense_form() {
        assert_eq!(
            parse_observable("IXZ", 3),
            Ok((false, vec![(1, Pauli::X), (2, Pauli::Z)]))
        );
        assert_eq!(parse_observable("-III", 3), Ok((true, vec![])));
    }

    #[test]
    fn dense_accepts_stim_underscore_identity() {
        // Stim's `PauliString` text writes identity as `_` (e.g. "+__X_Y_Z"),
        // so `_` is accepted interchangeably with `I`.
        assert_eq!(
            parse_observable("+__X", 3),
            Ok((false, vec![(2, Pauli::X)]))
        );
        assert_eq!(parse_observable("-_Y_", 3), Ok((true, vec![(1, Pauli::Y)])));
        assert_eq!(parse_observable("___", 3), Ok((false, vec![])));
        // Mixed `_` and `I` identities in one dense string.
        assert_eq!(parse_observable("I_Z", 3), Ok((false, vec![(2, Pauli::Z)])));
    }

    #[test]
    fn identity_variants() {
        assert_eq!(parse_observable("", 3), Ok((false, vec![])));
        assert_eq!(parse_observable("+", 3), Ok((false, vec![])));
        assert_eq!(parse_observable("-", 3), Ok((true, vec![])));
    }

    #[test]
    fn sparse_drops_identity_factors_but_tracks_qubit() {
        assert_eq!(
            parse_observable("I0*Z1", 2),
            Ok((false, vec![(1, Pauli::Z)]))
        );
        // Repeated qubit even when one factor is identity.
        assert_eq!(
            parse_observable("I0*Z0", 2),
            Err(ObservableParseError::RepeatedQubit(0))
        );
    }

    #[test]
    fn errors() {
        assert_eq!(
            parse_observable("Z5", 3),
            Err(ObservableParseError::QubitOutOfRange {
                qubit: 5,
                n_qubits: 3
            })
        );
        assert_eq!(
            parse_observable("Z0*Z0", 3),
            Err(ObservableParseError::RepeatedQubit(0))
        );
        assert!(matches!(
            parse_observable("Q0", 3),
            Err(ObservableParseError::BadToken(_))
        ));
        // Sparse factor missing its index.
        assert!(matches!(
            parse_observable("X0*Y", 3),
            Err(ObservableParseError::BadToken(_))
        ));
        // Dangling separator yields an empty factor.
        assert!(matches!(
            parse_observable("X0*", 3),
            Err(ObservableParseError::BadToken(_))
        ));
        // A bare letter run is dense; here it just mismatches the qubit count.
        assert_eq!(
            parse_observable("ZZ", 3),
            Err(ObservableParseError::DenseLengthMismatch {
                got: 2,
                n_qubits: 3
            })
        );
        // ... and is valid dense when the length matches.
        assert_eq!(parse_observable("X", 1), Ok((false, vec![(0, Pauli::X)])));
    }
}
