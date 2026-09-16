// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Error type for [`crate::LindbladSpec`] construction and stepping.

use crate::MAX_QUBITS;
use std::fmt;

/// Errors raised when constructing a [`crate::LindbladSpec`].
#[derive(Debug, Clone)]
pub enum Error {
    TooManyQubits {
        got: usize,
    },
    LengthMismatch {
        what: &'static str,
        a: usize,
        b: usize,
    },
    InvalidPauliCode {
        code: u8,
    },
    InvalidPauliChar {
        c: char,
    },
    WrongLength {
        expected: usize,
        got: usize,
    },
    NegativeRate {
        index: usize,
        rate: f64,
    },
    EmptyLincomb {
        index: usize,
    },
    /// Row `row` of the Kossakowski matrix is not `n_ops` wide.
    KMatrixRowLength {
        row: usize,
        expected: usize,
        got: usize,
    },
    /// `K_nm ≠ conj(K_mn)`: not a valid GKSL pair matrix.
    KMatrixNotHermitian {
        n: usize,
        m: usize,
    },
    Internal(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::TooManyQubits { got } => {
                write!(
                    f,
                    "LindbladSpec supports n_qubits ≤ {MAX_QUBITS}; got {got}"
                )
            }
            Error::LengthMismatch { what, a, b } => {
                write!(f, "{what}: expected matching lengths, got {a} and {b}")
            }
            Error::InvalidPauliCode { code } => write!(
                f,
                "Pauli code must be 0 (I), 1 (X), 2 (Z), or 3 (Y); got {code}"
            ),
            Error::InvalidPauliChar { c } => {
                write!(f, "invalid Pauli character '{c}'; expected I, X, Y, or Z")
            }
            Error::WrongLength { expected, got } => {
                write!(f, "Pauli string has length {got} but n_qubits = {expected}")
            }
            Error::NegativeRate { index, rate } => {
                write!(f, "jump rate must be non-negative; got γ_{index} = {rate}")
            }
            Error::EmptyLincomb { index } => {
                write!(
                    f,
                    "jump {index}: lincomb must contain at least one Pauli term"
                )
            }
            Error::KMatrixRowLength { row, expected, got } => write!(
                f,
                "kossakowski K row {row} has length {got}; expected {expected} (one per operator)"
            ),
            Error::KMatrixNotHermitian { n, m } => write!(
                f,
                "kossakowski K must be Hermitian; K[{n}][{m}] ≠ conj(K[{m}][{n}])"
            ),
            Error::Internal(msg) => write!(f, "internal error: {msg}"),
        }
    }
}

impl std::error::Error for Error {}
