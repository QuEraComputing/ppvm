use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::ops::Deref;

use ppvm_stim::{ExtendedProgram, parse_extended, prepare};

/// Python-facing wrapper around a prepared extended Stim program.
#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram(pub ExtendedProgram);

#[pymethods]
impl PyStimProgram {
    /// Parse and prepare a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let program = parse_extended(src).map_err(stim_to_pyerr)?;
        prepare(&program).map_err(stim_to_pyerr_exec)?;
        Ok(Self(program))
    }

    /// Read a `.stim` file and parse it.
    #[staticmethod]
    pub fn from_file(path: &str) -> PyResult<Self> {
        let src = std::fs::read_to_string(path)
            .map_err(|e| PyIOError::new_err(format!("failed to read {path}: {e}")))?;
        Self::parse(&src)
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

fn stim_to_pyerr(e: ppvm_stim::ExtendedParseError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_exec(e: ppvm_stim::ExecError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
