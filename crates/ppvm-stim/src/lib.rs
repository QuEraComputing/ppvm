// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parse and execute Stim circuits against a supported generalized tableau.
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
//! use ppvm_stim::backend::prelude::*;
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
//! [`GeneralizedTableau`]: backend::prelude::GeneralizedTableau

#[cfg(all(feature = "legacy", feature = "traits-2"))]
compile_error!("features `legacy` and `traits-2` are mutually exclusive");
#[cfg(not(any(feature = "legacy", feature = "traits-2")))]
compile_error!("enable exactly one ppvm-stim backend: `legacy` or `traits-2`");

/// Types and traits from the selected tableau backend.
pub mod backend {
    #[cfg(feature = "legacy")]
    pub use ppvm_pauli_sum_legacy as pauli_sum;
    #[cfg(feature = "traits-2")]
    pub use ppvm_tableau_2 as tableau;
    #[cfg(feature = "legacy")]
    pub use ppvm_tableau_legacy as tableau;

    /// Configuration names used by the selected backend.
    #[cfg(feature = "legacy")]
    pub use ppvm_pauli_sum_legacy::config;
    /// Storage aliases matching legacy call sites. The `-2` tableau stores
    /// dynamically sized rows, so the legacy byte-width parameter is inert.
    #[cfg(feature = "traits-2")]
    pub mod config {
        pub mod indexmap {
            /// Byte-packed row storage matching legacy `ByteFxHashF64<N>`.
            pub type ByteFxHashF64<const N: usize> = [u8; N];
        }
        pub mod fx64hash {
            /// Native-word row storage matching legacy `Byte8F64<N>`.
            pub type Byte8F64<const N: usize> = [usize; N];
        }
    }

    /// Convenience imports from the selected tableau backend.
    pub mod prelude {
        #[cfg(feature = "traits-2")]
        pub use ppvm_tableau_2::prelude::*;
        #[cfg(feature = "legacy")]
        pub use ppvm_tableau_legacy::prelude::*;
    }
}

pub mod executor;
pub mod validate;

pub use stim_parser::prelude::*;

pub use executor::{
    BackendTableau, StimTableau, execute, execute_validated, execute_validated_with_rng,
    execute_with_rng, sample, sample_serial, sample_serial_validated,
    sample_serial_validated_with_rng, sample_serial_with_rng, sample_validated,
    sample_validated_with_rng, sample_with_rng,
};
#[cfg(feature = "rayon")]
pub use executor::{
    sample_parallel, sample_parallel_validated, sample_parallel_validated_with_rng,
    sample_parallel_with_rng,
};
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
pub fn run_string<C, I, S>(
    src: &str,
    tab: &mut BackendTableau<C, I, S>,
) -> Result<Vec<Option<bool>>, Error>
where
    executor::SelectedBackend: executor::TableauType<C, I, S>,
{
    let prog = parse_extended(src)?;
    let results = execute(&prog, tab)?;
    Ok(results)
}

/// Parse, validate, and execute with caller-supplied randomness.
pub fn run_string_with_rng<C, I, S, R: rand::Rng + ?Sized>(
    src: &str,
    tab: &mut BackendTableau<C, I, S>,
    rng: &mut R,
) -> Result<Vec<Option<bool>>, Error>
where
    executor::SelectedBackend: executor::TableauType<C, I, S>,
{
    let prog = parse_extended(src)?;
    Ok(execute_with_rng(&prog, tab, rng)?)
}

pub fn run_file<C, I, S>(
    path: &Path,
    tab: &mut BackendTableau<C, I, S>,
) -> Result<Vec<Option<bool>>, Error>
where
    executor::SelectedBackend: executor::TableauType<C, I, S>,
{
    let src = std::fs::read_to_string(path).map_err(|source| Error::Io {
        path: path.to_path_buf(),
        source,
    })?;
    run_string(&src, tab)
}

/// Read, parse, validate, and execute with caller-supplied randomness.
pub fn run_file_with_rng<C, I, S, R: rand::Rng + ?Sized>(
    path: &Path,
    tab: &mut BackendTableau<C, I, S>,
    rng: &mut R,
) -> Result<Vec<Option<bool>>, Error>
where
    executor::SelectedBackend: executor::TableauType<C, I, S>,
{
    let src = std::fs::read_to_string(path).map_err(|source| Error::Io {
        path: path.to_path_buf(),
        source,
    })?;
    run_string_with_rng(&src, tab, rng)
}
