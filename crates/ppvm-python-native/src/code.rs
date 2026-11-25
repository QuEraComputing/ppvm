// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use ppvm_code::{Code, Row};
use pyo3::{exceptions::PyValueError, prelude::*};

type IndexedParities = Vec<(Vec<usize>, bool)>;

#[pyclass(frozen, module = "ppvm._core")]
pub struct CodeSpec(Code);

fn parse(rows: Vec<String>) -> PyResult<Vec<Row>> {
    rows.iter()
        .map(|row| row.parse().map_err(PyValueError::new_err))
        .collect()
}

#[pymethods]
impl CodeSpec {
    #[new]
    #[pyo3(signature = (n, stabilizers, logical_xs=None, logical_zs=None))]
    fn new(
        n: usize,
        stabilizers: Vec<String>,
        logical_xs: Option<Vec<String>>,
        logical_zs: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let code = Code::new(n, parse(stabilizers)?).map_err(PyValueError::new_err)?;
        Self::finish(code, logical_xs, logical_zs)
    }

    #[staticmethod]
    #[pyo3(signature = (xs, zs, signs, logical_xs=None, logical_zs=None))]
    fn from_tableau(
        xs: PyReadonlyArray2<'_, bool>,
        zs: PyReadonlyArray2<'_, bool>,
        signs: PyReadonlyArray1<'_, bool>,
        logical_xs: Option<Vec<String>>,
        logical_zs: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let x = xs.as_array();
        let z = zs.as_array();
        let signs = signs.as_array();
        if x.dim() != z.dim() || signs.len() != x.nrows() {
            return Err(PyValueError::new_err(
                "expected matching X/Z matrices and one sign per row",
            ));
        }
        let mut rows = Vec::with_capacity(x.nrows());
        for ((xr, zr), &sign) in x.rows().into_iter().zip(z.rows()).zip(signs) {
            rows.push(
                Row::from_bits(&xr.to_vec(), &zr.to_vec(), sign).map_err(PyValueError::new_err)?,
            );
        }
        let code = Code::new(x.ncols(), rows).map_err(PyValueError::new_err)?;
        Self::finish(code, logical_xs, logical_zs)
    }

    fn contains_stabilizer(&self, op: &str) -> PyResult<bool> {
        let row = op.parse().map_err(PyValueError::new_err)?;
        self.0
            .contains_stabilizer(&row)
            .map_err(PyValueError::new_err)
    }

    fn terminal_readout(
        &self,
        axes: &str,
        logical_axes: &str,
    ) -> PyResult<(IndexedParities, IndexedParities)> {
        let (checks, obs) = self
            .0
            .terminal_readout(axes.as_bytes(), logical_axes.as_bytes())
            .map_err(PyValueError::new_err)?;
        let pack =
            |rows: Vec<ppvm_code::Parity>| rows.into_iter().map(|p| (p.qubits, p.offset)).collect();
        Ok((pack(checks), pack(obs)))
    }

    #[getter]
    fn n(&self) -> usize {
        self.0.n()
    }
    #[getter]
    fn k(&self) -> usize {
        self.0.k()
    }
    #[getter]
    fn stabilizers(&self) -> Vec<String> {
        self.0
            .stabilizers()
            .iter()
            .map(ToString::to_string)
            .collect()
    }
    #[getter]
    fn independent_stabilizers(&self) -> Vec<String> {
        self.0
            .independent_stabilizers()
            .map(ToString::to_string)
            .collect()
    }
    #[getter]
    fn logical_xs(&self) -> Vec<String> {
        self.0
            .logical_pairs()
            .iter()
            .map(|(x, _)| x.to_string())
            .collect()
    }
    #[getter]
    fn logical_zs(&self) -> Vec<String> {
        self.0
            .logical_pairs()
            .iter()
            .map(|(_, z)| z.to_string())
            .collect()
    }
}

impl CodeSpec {
    fn finish(
        code: Code,
        logical_xs: Option<Vec<String>>,
        logical_zs: Option<Vec<String>>,
    ) -> PyResult<Self> {
        match (logical_xs, logical_zs) {
            (None, None) => Ok(Self(code)),
            (Some(xs), Some(zs)) => code
                .with_logicals(parse(xs)?, parse(zs)?)
                .map(Self)
                .map_err(PyValueError::new_err),
            _ => Err(PyValueError::new_err(
                "supply both logical_xs and logical_zs",
            )),
        }
    }
}
