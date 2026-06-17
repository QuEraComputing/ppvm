// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! PyO3 wrapper around [`ppvm_lindblad::LindbladSpec`].
//!
//! All algorithmic work — Pauli arithmetic, active-site iteration, cache
//! management, and the dissipator branches (Hermitian Pauli fast path vs
//! general complex Pauli sum) — lives in the [`ppvm_lindblad`] crate. This
//! module is responsible only for the Python boundary: decoding the
//! `(N, n_qubits)` numpy uint8 arrays into [`ppvm_lindblad::Word`] vectors,
//! and re-encoding outputs back into numpy.

use std::collections::HashMap;

use num::Complex;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use ppvm_lindblad::{
    JumpInput, LindbladSpec as CoreSpec, PcStepConfig, Word, codes_from_word, word_from_codes,
};
use pyo3::{exceptions::PyValueError, prelude::*};

type PyPauliMap<'py> = (Bound<'py, PyArray2<u8>>, Bound<'py, PyArray1<f64>>);
type PyCoo<'py> = (
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<f64>>,
);

fn map_err(e: ppvm_lindblad::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Reject a basis that contains the same Pauli word at two distinct rows.
/// Duplicate rows would silently overwrite each other in the generator's
/// row-index map and produce an incorrect sparse matrix.
fn assert_basis_unique(basis: &[Word]) -> PyResult<()> {
    let mut seen: HashMap<&Word, usize> = HashMap::with_capacity(basis.len());
    for (i, w) in basis.iter().enumerate() {
        if let Some(prev) = seen.insert(w, i) {
            return Err(PyValueError::new_err(format!(
                "basis contains duplicate Pauli word at row {prev} and row {i}"
            )));
        }
    }
    Ok(())
}

/// Decode a `(N, n_qubits)` uint8 ndarray view into `N` packed [`Word`]s.
fn decode_basis(view: &numpy::ndarray::ArrayView2<u8>, n_qubits: usize) -> PyResult<Vec<Word>> {
    let n_basis = view.shape()[0];
    let n_cols = view.shape()[1];
    if n_cols != n_qubits {
        return Err(PyValueError::new_err(format!(
            "basis has {n_cols} columns but spec.n_qubits = {n_qubits}"
        )));
    }
    let mut out = Vec::with_capacity(n_basis);
    let mut row_buf = vec![0u8; n_qubits];
    for i in 0..n_basis {
        let row = view.row(i);
        for (q, slot) in row_buf.iter_mut().enumerate() {
            *slot = row[q];
        }
        out.push(word_from_codes(&row_buf).map_err(map_err)?);
    }
    Ok(out)
}

/// Pack `Vec<(Word, f64)>` into the standard PyO3 return shape.
fn pack_pauli_map<'py>(
    py: Python<'py>,
    pairs: Vec<(Word, f64)>,
    n_qubits: usize,
) -> PyResult<PyPauliMap<'py>> {
    let m = pairs.len();
    let mut basis = vec![0u8; m * n_qubits];
    let mut coeffs = vec![0f64; m];
    for (i, (w, c)) in pairs.into_iter().enumerate() {
        codes_from_word(&w, &mut basis[i * n_qubits..(i + 1) * n_qubits]);
        coeffs[i] = c;
    }
    let basis_arr = basis
        .into_pyarray(py)
        .reshape([m, n_qubits])
        .map_err(|e| PyValueError::new_err(format!("reshape failed: {e}")))?;
    Ok((basis_arr, coeffs.into_pyarray(py)))
}

/// PyO3 facade exposing the Lindbladian shim to Python.
#[pyclass]
pub struct LindbladSpec {
    inner: CoreSpec,
}

#[pymethods]
impl LindbladSpec {
    /// Construct a Lindbladian spec from Hamiltonian terms and jump operators.
    ///
    /// `jump_lincombs[k]` is a list of `(pauli_string, real, imag)` triples
    /// encoding `L_k = Σ_a (re + i·im) P_a`. A length-1 jump with `im == 0`
    /// is routed to the Hermitian-Pauli fast path (with rate scaled by `re²`).
    #[new]
    #[pyo3(signature = (n_qubits, h_terms, h_coeffs, jump_lincombs, jump_rates))]
    fn new(
        n_qubits: usize,
        h_terms: Vec<String>,
        h_coeffs: Vec<f64>,
        jump_lincombs: Vec<Vec<(String, f64, f64)>>,
        jump_rates: Vec<f64>,
    ) -> PyResult<Self> {
        if h_terms.len() != h_coeffs.len() {
            return Err(PyValueError::new_err(
                "h_terms and h_coeffs must have the same length",
            ));
        }
        if jump_lincombs.len() != jump_rates.len() {
            return Err(PyValueError::new_err(
                "jump_lincombs and jump_rates must have the same length",
            ));
        }
        let h: Vec<(String, f64)> = h_terms.into_iter().zip(h_coeffs).collect();
        let jumps: Vec<JumpInput> = jump_lincombs
            .into_iter()
            .zip(jump_rates)
            .map(|(lincomb, rate)| JumpInput {
                lincomb: lincomb
                    .into_iter()
                    .map(|(s, re, im)| (s, Complex::new(re, im)))
                    .collect(),
                rate,
            })
            .collect();
        let inner = CoreSpec::new(n_qubits, &h, &jumps).map_err(map_err)?;
        Ok(Self { inner })
    }

    #[getter]
    fn n_qubits(&self) -> usize {
        self.inner.n_qubits()
    }

    #[getter]
    fn num_h_terms(&self) -> usize {
        self.inner.num_h_terms()
    }

    #[getter]
    fn num_jump_terms(&self) -> usize {
        self.inner.num_jump_terms()
    }

    /// Apply `L*` to a single Pauli string `p`.
    fn action<'py>(
        &self,
        py: Python<'py>,
        p: PyReadonlyArray1<'py, u8>,
    ) -> PyResult<PyPauliMap<'py>> {
        let p_slice = p.as_slice()?;
        let p_word = word_from_codes(p_slice).map_err(map_err)?;
        let pairs = self.inner.action(&p_word);
        pack_pauli_map(py, pairs, self.inner.n_qubits())
    }

    /// Off-basis component of `L*( Σ_j coeffs[j] · basis[j] )`.
    #[pyo3(signature = (basis, coeffs, protected = None))]
    fn leakage<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
        coeffs: PyReadonlyArray1<'py, f64>,
        protected: Option<PyReadonlyArray2<'py, u8>>,
    ) -> PyResult<PyPauliMap<'py>> {
        let n_q = self.inner.n_qubits();
        let basis_view = basis.as_array();
        let basis_words = decode_basis(&basis_view, n_q)?;
        let coeffs_slice = coeffs.as_slice()?;
        if coeffs_slice.len() != basis_words.len() {
            return Err(PyValueError::new_err(format!(
                "coeffs has length {} but basis has {} rows",
                coeffs_slice.len(),
                basis_words.len()
            )));
        }
        let protected_words: Vec<Word> = if let Some(ref prot) = protected {
            let pv = prot.as_array();
            decode_basis(&pv, n_q)?
        } else {
            Vec::new()
        };
        let pairs = self
            .inner
            .leakage(&basis_words, coeffs_slice, &protected_words)
            .map_err(map_err)?;
        pack_pauli_map(py, pairs, n_q)
    }

    /// One predictor-corrector adaptive step.
    ///
    /// Internally: expand basis with first-hop leakage, predictor step
    /// (`exp(dt·M)`), expand again with second-hop leakage from the
    /// predicted state, then redo the step from the pre-step coefficients
    /// on the doubly-enlarged basis. The matrix exponential is computed in
    /// Rust via Al-Mohy & Higham scaling-and-squaring; no scipy required.
    ///
    /// Returns `(new_basis, new_coeffs)`.
    #[pyo3(signature = (
        basis, coeffs, dt, tau_add,
        drop_tol = 0.0,
        protected = None,
        num_threads = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn pc_step<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
        coeffs: PyReadonlyArray1<'py, f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: Option<PyReadonlyArray2<'py, u8>>,
        num_threads: Option<usize>,
    ) -> PyResult<PyPauliMap<'py>> {
        let n_q = self.inner.n_qubits();
        let basis_view = basis.as_array();
        let mut basis_words = decode_basis(&basis_view, n_q)?;
        assert_basis_unique(&basis_words)?;
        let mut coeffs_vec = coeffs.as_slice()?.to_vec();
        if coeffs_vec.len() != basis_words.len() {
            return Err(PyValueError::new_err(format!(
                "coeffs has length {} but basis has {} rows",
                coeffs_vec.len(),
                basis_words.len()
            )));
        }
        let protected_words: Vec<Word> = if let Some(ref p) = protected {
            decode_basis(&p.as_array(), n_q)?
        } else {
            Vec::new()
        };
        self.inner
            .pc_step(
                &mut basis_words,
                &mut coeffs_vec,
                dt,
                &protected_words,
                PcStepConfig {
                    tau_add,
                    drop_tol,
                    num_threads,
                },
            )
            .map_err(map_err)?;

        // Pack output. Basis may have grown; coeffs has the same new length.
        let pairs: Vec<(Word, f64)> = basis_words.into_iter().zip(coeffs_vec).collect();
        pack_pauli_map(py, pairs, n_q)
    }

    /// Same as [`Self::pc_step`] but also returns a dict mapping phase
    /// name → microseconds spent in that phase, for profiling.
    #[pyo3(signature = (
        basis, coeffs, dt, tau_add,
        drop_tol = 0.0,
        protected = None,
        num_threads = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn pc_step_timed<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
        coeffs: PyReadonlyArray1<'py, f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: Option<PyReadonlyArray2<'py, u8>>,
        num_threads: Option<usize>,
    ) -> PyResult<(PyPauliMap<'py>, Bound<'py, pyo3::types::PyDict>)> {
        let n_q = self.inner.n_qubits();
        let basis_view = basis.as_array();
        let mut basis_words = decode_basis(&basis_view, n_q)?;
        assert_basis_unique(&basis_words)?;
        let mut coeffs_vec = coeffs.as_slice()?.to_vec();
        if coeffs_vec.len() != basis_words.len() {
            return Err(PyValueError::new_err(format!(
                "coeffs has length {} but basis has {} rows",
                coeffs_vec.len(),
                basis_words.len()
            )));
        }
        let protected_words: Vec<Word> = if let Some(ref p) = protected {
            decode_basis(&p.as_array(), n_q)?
        } else {
            Vec::new()
        };
        let timings = self
            .inner
            .pc_step_timed(
                &mut basis_words,
                &mut coeffs_vec,
                dt,
                &protected_words,
                PcStepConfig {
                    tau_add,
                    drop_tol,
                    num_threads,
                },
            )
            .map_err(map_err)?;

        let pairs: Vec<(Word, f64)> = basis_words.into_iter().zip(coeffs_vec).collect();
        let map = pack_pauli_map(py, pairs, n_q)?;
        let d = pyo3::types::PyDict::new(py);
        d.set_item("leakage1_us", timings.leakage1_us)?;
        d.set_item("expand1_us", timings.expand1_us)?;
        d.set_item("gencsr1_us", timings.gencsr1_us)?;
        d.set_item("expm1_us", timings.expm1_us)?;
        d.set_item("leakage2_us", timings.leakage2_us)?;
        d.set_item("expand2_us", timings.expand2_us)?;
        d.set_item("gencsr2_us", timings.gencsr2_us)?;
        d.set_item("expm2_us", timings.expm2_us)?;
        Ok((map, d))
    }

    /// One classical RK4 step on the adjoint Lindbladian. Matrix-free: no
    /// CSR, no Krylov, no predictor-corrector enrichment. Four action
    /// evaluations per step, basis grows naturally, magnitude-prune at end.
    #[pyo3(signature = (
        basis, coeffs, dt,
        drop_tol = 0.0,
        protected = None,
        num_threads = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn rk4_step<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
        coeffs: PyReadonlyArray1<'py, f64>,
        dt: f64,
        drop_tol: f64,
        protected: Option<PyReadonlyArray2<'py, u8>>,
        num_threads: Option<usize>,
    ) -> PyResult<PyPauliMap<'py>> {
        let n_q = self.inner.n_qubits();
        let basis_view = basis.as_array();
        let mut basis_words = decode_basis(&basis_view, n_q)?;
        assert_basis_unique(&basis_words)?;
        let mut coeffs_vec = coeffs.as_slice()?.to_vec();
        if coeffs_vec.len() != basis_words.len() {
            return Err(PyValueError::new_err(format!(
                "coeffs has length {} but basis has {} rows",
                coeffs_vec.len(),
                basis_words.len()
            )));
        }
        let protected_words: Vec<Word> = if let Some(ref p) = protected {
            decode_basis(&p.as_array(), n_q)?
        } else {
            Vec::new()
        };
        self.inner
            .rk4_step(
                &mut basis_words,
                &mut coeffs_vec,
                dt,
                drop_tol,
                &protected_words,
                num_threads,
            )
            .map_err(map_err)?;
        let pairs: Vec<(Word, f64)> = basis_words.into_iter().zip(coeffs_vec).collect();
        pack_pauli_map(py, pairs, n_q)
    }

    /// Sparse generator matrix in COO form: `(rows, cols, vals)`.
    fn generator<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
    ) -> PyResult<PyCoo<'py>> {
        let n_q = self.inner.n_qubits();
        let basis_view = basis.as_array();
        let basis_words = decode_basis(&basis_view, n_q)?;
        assert_basis_unique(&basis_words)?;
        let triplets = self.inner.generator(&basis_words);
        let total = triplets.len();
        let mut rows = Vec::with_capacity(total);
        let mut cols = Vec::with_capacity(total);
        let mut vals = Vec::with_capacity(total);
        for (r, c, v) in triplets {
            rows.push(r as u64);
            cols.push(c as u64);
            vals.push(v);
        }
        Ok((
            rows.into_pyarray(py),
            cols.into_pyarray(py),
            vals.into_pyarray(py),
        ))
    }
}
