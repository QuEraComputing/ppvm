// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Exact sampling of CSS circuits with atom loss and leakage, many shots per machine word.
//!
//! A circuit is *CSS up to Hadamards* when conjugating some of its qubit segments (from one
//! reset to the next) by `H` turns every two-qubit gate into a `CX`. Preparations and
//! measurements may use either basis. Following Gidney & Jones, *CSS Circuits: useful quantum
//! circuits that can be simulated with O(1) cost per gate* ([arXiv:2611.14944]), the all-zero
//! record is then a valid noiseless outcome for every pattern of gates cancelled by loss, so a
//! Pauli frame per shot samples loss and leakage exactly. Gates touching a lost or leaked qubit
//! are skipped, with the same semantics as the generalized-tableau backend ([`sample`]).
//!
//! Supported: `R`, `RX`, `M`, `MX`, `MR`, `MRX`, `MPAD`, `X`, `Y`, `Z`, `H`, `SQRT_Y`,
//! `SQRT_Y_DAG`, `CX`, `CZ`, classical feedback (`CX` / `CY` / `CZ` controlled by `rec[-k]`),
//! Pauli noise, measurement flips, `REPEAT`, and the `I_ERROR` loss, leakage and reset tags.
//!
//! ```
//! use ppvm_stim::{css::CssProgram, parse_extended};
//!
//! let program = parse_extended("R 0 1\nH 0\nI_ERROR[loss](0.1) 0\nCX 0 1\nM 0 1")?;
//! let css = CssProgram::compile(&program)?;
//! let shots = css.sample(1000, 7); // 1000 x 2 bytes: 0 / 1 outcome, 2 = lost
//! assert_eq!(shots.len(), 2000);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! [arXiv:2611.14944]: https://arxiv.org/abs/2611.14944
//! [`sample`]: crate::sample

mod bits;
mod compile;
mod sample;

use crate::validate::ExecError;

/// A Stim program compiled for the CSS frame sampler; see the [module docs](self).
#[derive(Debug, Clone)]
pub struct CssProgram {
    ops: Vec<compile::Op>,
    num_qubits: usize,
    num_measurements: usize,
}

impl CssProgram {
    /// Number of measurement-record entries per shot.
    pub fn num_measurements(&self) -> usize {
        self.num_measurements
    }
}

/// Why a program cannot run on the CSS frame sampler.
#[derive(Debug, PartialEq, thiserror::Error)]
pub enum CssError {
    /// Invalid program, or an instruction outside the CSS gate set.
    #[error(transparent)]
    Exec(#[from] ExecError),
    /// The two-qubit gates admit no consistent choice of Hadamard-conjugated segments.
    #[error("the circuit is not CSS up to Hadamard gates (line {line})")]
    NotCss { line: usize },
}
