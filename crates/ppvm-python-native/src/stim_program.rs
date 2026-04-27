use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;

use ppvm_stim::{Error as StimError, TableauProgram, parse, normalize};

/// Python-facing wrapper around `ppvm_stim::TableauProgram`.
#[pyclass(name = "StimProgram", module = "ppvm_python_native")]
pub struct PyStimProgram {
    pub(crate) inner: TableauProgram,
}

#[pymethods]
impl PyStimProgram {
    /// Parse + normalize a Stim circuit string.
    #[staticmethod]
    pub fn parse(src: &str) -> PyResult<Self> {
        let prog = parse(src).map_err(stim_to_pyerr)?;
        let tprog = normalize::to_tableau(&prog).map_err(stim_to_pyerr_norm)?;
        Ok(Self { inner: tprog })
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
            self.inner.instructions.len(),
            self.inner.expected_measurement_count
        )
    }
}

fn stim_to_pyerr(e: ppvm_stim::ParseError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

fn stim_to_pyerr_norm(e: ppvm_stim::NormalizeError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}

#[allow(dead_code)]
fn full_to_pyerr(e: StimError) -> PyErr {
    PyValueError::new_err(format!("{e}"))
}
