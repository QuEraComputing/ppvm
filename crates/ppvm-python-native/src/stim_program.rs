// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyByteArray;
use std::ops::Deref;
use std::sync::OnceLock;

use ppvm_stim::css::{CssError, CssProgram};
use ppvm_stim::{
    ExtendedProgram, StimPrint, highlight_ansi, highlight_html, parse_extended, validate,
};

/// Python-facing wrapper around a validated extended Stim program, with its CSS frame-sampler
/// compilation built on first use.
#[pyclass(name = "StimProgram", module = "ppvm._core")]
pub struct PyStimProgram {
    pub program: ExtendedProgram,
    css: OnceLock<Result<CssProgram, CssError>>,
}

#[pymethods]
impl PyStimProgram {
    /// Parse and validate a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr)?;
        validate(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self {
            program,
            css: OnceLock::new(),
        })
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
        self.program.to_stim()
    }

    fn __repr__(&self) -> String {
        format!(
            "<StimProgram instructions={} measurements={}>",
            self.program.instructions.len(),
            self.program.measurement_count()
        )
    }

    /// Number of qubits the program operates on: one past the highest qubit
    /// index any instruction references (`0` if it touches no qubits). Lets
    /// callers size a tableau when no explicit qubit count is given.
    #[getter]
    fn num_qubits(&self) -> usize {
        self.program.num_qubits()
    }

    /// Whether the CSS frame sampler can run this program: it is CSS up to Hadamard gates and
    /// uses only the gates, noise, loss and leakage that sampler supports. A CSS circuit that uses
    /// other instructions, such as `MPP`, is not covered.
    #[getter]
    fn is_css(&self) -> bool {
        self.css().is_ok()
    }

    /// Sample with the CSS frame sampler. Returns the outcome codes (0/1/2 = zero/one/lost) as
    /// one flat, shot-major `bytearray` plus its `(num_shots, n_measurements)` shape.
    #[pyo3(signature = (num_shots, seed = None))]
    fn _sample_css<'py>(
        &self,
        py: Python<'py>,
        num_shots: usize,
        seed: Option<u64>,
    ) -> PyResult<(Bound<'py, PyByteArray>, (usize, usize))> {
        let css = self
            .css()
            .as_ref()
            .map_err(|e| PyValueError::new_err(format!("not supported by the CSS sampler: {e}")))?;
        let seed = seed.unwrap_or_else(rand::random);
        let m = css.num_measurements();
        let len = num_shots
            .checked_mul(m)
            .ok_or_else(|| PyValueError::new_err("num_shots is too large"))?;
        let shots = PyByteArray::new_with(py, len, |out| {
            py.detach(|| css.sample_into(out, seed));
            Ok(())
        })?;
        Ok((shots, (num_shots, m)))
    }

    /// Jupyter rich display: syntax-highlighted Stim source. Only invoked in
    /// IPython/Jupyter; plain `str()`/`print()` stay uncoloured elsewhere.
    fn _repr_html_(&self) -> String {
        highlight_html(&self.program.to_stim())
    }

    /// IPython terminal pretty-printer: writes ANSI-coloured Stim source.
    fn _repr_pretty_(&self, printer: &Bound<'_, PyAny>, _cycle: bool) -> PyResult<()> {
        printer.call_method1("text", (highlight_ansi(&self.program.to_stim()),))?;
        Ok(())
    }
}

impl PyStimProgram {
    fn css(&self) -> &Result<CssProgram, CssError> {
        self.css.get_or_init(|| CssProgram::compile(&self.program))
    }
}

impl Deref for PyStimProgram {
    type Target = ExtendedProgram;
    fn deref(&self) -> &ExtendedProgram {
        &self.program
    }
}

fn stim_to_pyerr(e: ppvm_stim::Diagnostics) -> PyErr {
    // Rich, source-pointing report (offending line + caret), plain text.
    PyValueError::new_err(e.render())
}

fn stim_to_pyerr_exec(e: ppvm_stim::ExecError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
