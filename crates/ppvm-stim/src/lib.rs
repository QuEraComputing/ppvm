//! Stim circuit format parser and executor for the `ppvm` simulator.
//!
//! Three-stage pipeline (filled in across Tasks 2–14):
//!
//! 1. `parse` — `&str` → `Program` (pure-Stim AST, tags preserved).
//! 2. `normalize::to_tableau` — `Program` → `TableauProgram` (dialect resolved,
//!    phase-1-unsupported instructions rejected).
//! 3. `execute` — apply a `TableauProgram` to a `GeneralizedTableau`.
//!
//! Multi-shot usage should call `parse` and `normalize::to_tableau` once and
//! call `sample` for the shot loop. The `run_string` / `run_file`
//! convenience helpers re-parse on every call and are intended for
//! single-shot demos only.
//!
//! # Multi-shot pattern (recommended)
//!
//! ```ignore
//! use ppvm_stim::{parse, normalize, sample};
//! use ppvm_tableau::prelude::*;
//!
//! let prog = parse(circuit_src)?;
//! let tprog = normalize::to_tableau(&prog)?;
//! let shots = sample(&tprog, 10_000, || {
//!     GeneralizedTableau::<_, usize, _>::new(n_qubits, 1e-10)
//! })?;
//! # Ok::<(), ppvm_stim::Error>(())
//! ```
//!
//! [`run_string`] / [`run_file`] re-parse on every call and exist only for
//! single-shot demos — never call them from a shot loop.

pub mod parser;
pub mod tableau_program;
pub mod normalize;
pub mod executor;

pub use parser::ast::{
    AnnotationKind, GateName, MeasureName, NoiseName, ParseError, Program,
    RawInstruction, Tag, TagParam,
};
pub use parser::parse;

pub use tableau_program::{
    GateKind, Instruction, MeasureKind, NoiseKind, TableauProgram,
};

pub use normalize::NormalizeError;

pub use executor::{ExecError, execute, sample};

use std::path::Path;

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error(transparent)]
    Normalize(#[from] NormalizeError),
    #[error(transparent)]
    Exec(#[from] ExecError),
}

/// Parse → normalize → execute in one shot. Re-parses each call; do **not**
/// use in shot loops — use [`parse`] + [`normalize::to_tableau`] + [`sample`]
/// instead.
pub fn run_string<T, I, C>(
    src: &str,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_runtime::prelude::Config,
    <<T as ppvm_runtime::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One + num::Zero + Clone + num::Num + num::ToPrimitive
        + std::fmt::Debug + std::ops::Mul<f64> + PartialOrd<f64> + Send + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64> + std::ops::MulAssign + std::ops::AddAssign
        + num::One + num::complex::ComplexFloat + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let prog = parse(src)?;
    let tprog = normalize::to_tableau(&prog)?;
    let results = execute(&tprog, tab)?;
    Ok(results)
}

pub fn run_file<T, I, C>(
    path: &Path,
    tab: &mut ppvm_tableau::prelude::GeneralizedTableau<T, I, C>,
) -> Result<Vec<Option<bool>>, Error>
where
    T: ppvm_runtime::prelude::Config,
    <<T as ppvm_runtime::prelude::Config>::Storage as bitvec::view::BitView>::Store: num::PrimInt,
    C: ppvm_tableau::prelude::SparseVector<num::Complex<T::Coeff>, I> + std::fmt::Debug,
    T::Coeff: num::One + num::Zero + Clone + num::Num + num::ToPrimitive
        + std::fmt::Debug + std::ops::Mul<f64> + PartialOrd<f64> + Send + Sync,
    num::Complex<T::Coeff>: std::ops::Mul<Output = num::Complex<T::Coeff>>
        + From<num::complex::Complex64> + std::ops::MulAssign + std::ops::AddAssign
        + num::One + num::complex::ComplexFloat + Copy,
    I: ppvm_tableau::prelude::TableauIndex + std::fmt::Debug + Send + Sync,
{
    let src = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("failed to read stim file {}: {e}", path.display()));
    run_string(&src, tab)
}
