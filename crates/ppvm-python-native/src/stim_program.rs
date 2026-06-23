// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::ops::Deref;

use ppvm_stim::{
    ExtendedProgram, StimPrint, highlight_ansi, highlight_html, parse_extended, validate,
};

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

    /// `str(program)` / `print(program)` yield canonical, round-trippable
    /// Stim text: `StimProgram.parse(str(p))` reproduces `p`. Whitespace,
    /// comments, and number spelling are normalized to the canonical form.
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

    /// Jupyter rich display: syntax-highlighted Stim source. Only invoked in
    /// IPython/Jupyter; plain `str()`/`print()` stay uncoloured elsewhere.
    fn _repr_html_(&self) -> String {
        highlight_html(&self.0.to_stim())
    }

    /// IPython terminal pretty-printer: writes ANSI-coloured Stim source.
    fn _repr_pretty_(&self, printer: &Bound<'_, PyAny>, _cycle: bool) -> PyResult<()> {
        printer.call_method1("text", (highlight_ansi(&self.0.to_stim()),))?;
        Ok(())
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
