// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::ops::Deref;

use ppvm_stim::{ExtendedProgram, StimPrint, parse_extended, validate};

/// Python-facing wrapper around a validated extended Stim program.
#[pyclass(name = "StimProgram", module = "ppvm._core")]
pub struct PyStimProgram(pub ExtendedProgram);

#[pymethods]
impl PyStimProgram {
    /// Parse and validate a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr)?;
        validate(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self(program))
    }

    /// Read a `.stim` file and parse it.
    #[staticmethod]
    pub fn from_file(path: &str) -> PyResult<Self> {
        let src = std::fs::read_to_string(path)
            .map_err(|e| PyIOError::new_err(format!("failed to read {path}: {e}")))?;
        Self::parse(&src)
    }

    /// Render the program back to canonical Stim text.
    ///
    /// Parsing is round-trippable: `StimProgram.parse(s).to_stim()` reparses
    /// to the same string (a print/parse fixpoint). The output is normalized
    /// canonical Stim — comments and original whitespace/number spelling are
    /// not preserved.
    fn to_stim(&self) -> String {
        self.0.to_stim()
    }

    /// `str(program)` yields the canonical Stim text (same as `to_stim()`).
    fn __str__(&self) -> String {
        self.0.to_stim()
    }

    fn __repr__(&self) -> String {
        format!(
            "<StimProgram instructions={} measurements={}>",
            self.0.instructions.len(),
            self.0.measurement_count()
        )
    }
}

impl Deref for PyStimProgram {
    type Target = ExtendedProgram;
    fn deref(&self) -> &ExtendedProgram {
        &self.0
    }
}

fn stim_to_pyerr(e: ppvm_stim::Diagnostics) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_exec(e: ppvm_stim::ExecError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
