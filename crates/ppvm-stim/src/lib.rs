// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parse and execute Stim circuits against a [`GeneralizedTableau`].
//!
//! Two-stage pipeline:
//!
//! 1. [`parse_extended`] — `&str` → [`ExtendedProgram`] (re-exported from
//!    `stim_parser`).
//! 2. [`execute`] / [`sample`] — call [`validate`](fn@validate) to validate the
//!    [`ExtendedProgram`], then apply it to a [`GeneralizedTableau`].
//!
//! Multi-shot usage should call [`parse_extended`] once and pass the parsed
//! program to [`sample`]. The [`run_string`] / [`run_file`] convenience helpers
//! re-parse on every call and are intended for single-shot demos only.
//!
//! # Multi-shot pattern (recommended)
//!
//! ```ignore
//! use ppvm_stim::{parse_extended, sample};
//! use ppvm_tableau::prelude::*;
//!
//! let prog = parse_extended(circuit_src)?;
//! let shots = sample(&prog, 10_000, |_| {
//!     GeneralizedTableau::<_, usize, _>::new(n_qubits, 1e-10)
//! })?;
//! # Ok::<(), ppvm_stim::Error>(())
//! ```
//!
//! [`run_string`] / [`run_file`] re-parse on every call and exist only for
//! single-shot demos — never call them from a shot loop.
//!
//! [`ExtendedProgram`]: stim_parser::prelude::ExtendedProgram
//! [`GeneralizedTableau`]: ppvm_tableau::prelude::GeneralizedTableau

pub mod executor;
pub mod validate;

pub use stim_parser::prelude::*;

pub use executor::{
    execute, execute_validated, sample, sample_cached, sample_cached_validated, sample_serial,
    sample_serial_validated, sample_validated,
};
#[cfg(feature = "rayon")]
pub use executor::{sample_parallel, sample_parallel_validated};
pub use validate::{ExecError, validate};

use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A parse/validate/lower failure from `stim_parser`, reported as a
    /// [`Diagnostics`] aggregate.
    #[error(transparent)]
    Parse(#[from] Diagnostics),
    #[error(transparent)]
    Exec(#[from] ExecError),
    #[error("failed to read stim file {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
}

/// Parse → validate → execute in one shot. Re-parses each call; do **not**
/// use in shot loops — use [`parse_extended`] + [`sample`] instead.
pub fn run_string<T, I, C>(
    src: &str,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_pauli_sum::prelude::Config,
    <<T as ppvm_pauli_sum::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One
        + num::Zero
        + Clone
        + num::Num
        + num::ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + PartialOrd
        + Send
        + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + num::One
        + num::complex::ComplexFloat
        + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let prog = parse_extended(src)?;
    let results = execute(&prog, tab)?;
    Ok(results)
}

pub fn run_file<T, I, C>(
    path: &Path,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_pauli_sum::prelude::Config,
    <<T as ppvm_pauli_sum::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One
        + num::Zero
        + Clone
        + num::Num
        + num::ToPrimitive
        + std::fmt::Debug
        + std::ops::Mul<f64>
        + PartialOrd<f64>
        + PartialOrd
        + Send
        + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + num::One
        + num::complex::ComplexFloat
        + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let src = std::fs::read_to_string(path).map_err(|source| Error::Io {
        path: path.to_path_buf(),
        source,
    })?;
    run_string(&src, tab)
}
