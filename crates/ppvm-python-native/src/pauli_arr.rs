// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Shared codec and argument validation for the `(N, n_qubits)` uint8
//! Pauli-basis array representation used across the Lindblad and
//! symmetry bindings.
//!
//! Every `*_arr` entry point decodes an incoming basis array into packed
//! [`Word`]s, validates the companion `coeffs` / `momentum` lengths, and
//! re-encodes the result on the way out. These are those three steps.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use ppvm_lindblad::{Word, codes_from_word, word_from_codes};
use pyo3::{exceptions::PyValueError, prelude::*};

use crate::lindblad::map_err;

/// Decode a `(N, n_qubits)` uint8 ndarray view into `N` packed [`Word`]s.
pub(crate) fn decode_basis(
    view: &numpy::ndarray::ArrayView2<u8>,
    n_qubits: usize,
) -> PyResult<Vec<Word>> {
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

/// Encode packed [`Word`]s back into an `(M, n_qubits)` uint8 array.
pub(crate) fn encode_basis<'py>(
    py: Python<'py>,
    words: &[Word],
    n_qubits: usize,
) -> PyResult<Bound<'py, PyArray2<u8>>> {
    let m = words.len();
    let mut flat = vec![0u8; m * n_qubits];
    for (i, w) in words.iter().enumerate() {
        codes_from_word(w, &mut flat[i * n_qubits..(i + 1) * n_qubits]);
    }
    flat.into_pyarray(py)
        .reshape([m, n_qubits])
        .map_err(|e| PyValueError::new_err(format!("reshape failed: {e}")))
}

/// Check the row width of a basis array against the qubit count a
/// [`crate::symmetry::TranslationGroup`] acts on. Reported separately from
/// [`decode_basis`]'s own width check so the error names the group rather
/// than the spec.
pub(crate) fn check_group_width(
    view: &numpy::ndarray::ArrayView2<u8>,
    n_qubits: usize,
) -> PyResult<()> {
    let width = view.shape().get(1).copied();
    if width != Some(n_qubits) {
        return Err(PyValueError::new_err(format!(
            "basis has {} qubits per row but group acts on {n_qubits}",
            width.unwrap_or(0)
        )));
    }
    Ok(())
}

/// Check that a [`crate::symmetry::TranslationGroup`] acts on the same
/// qubit count as the object being evolved. The core group routines
/// assert this on the first Pauli word they see, so without this the
/// mismatch surfaces as a panic from deep inside the step.
pub(crate) fn check_group_qubits(n_qubits: usize, group_n_qubits: usize) -> PyResult<()> {
    if n_qubits != group_n_qubits {
        return Err(PyValueError::new_err(format!(
            "spec has {n_qubits} qubits but the TranslationGroup acts on {group_n_qubits}"
        )));
    }
    Ok(())
}

/// Check that a coefficient vector has one entry per basis row.
pub(crate) fn check_coeffs_len(n_coeffs: usize, n_rows: usize) -> PyResult<()> {
    if n_coeffs != n_rows {
        return Err(PyValueError::new_err(format!(
            "coeffs has length {n_coeffs} but basis has {n_rows} rows"
        )));
    }
    Ok(())
}

/// Check that a momentum vector has one mode index per group generator.
pub(crate) fn check_momentum_len(n_modes: usize, n_generators: usize) -> PyResult<()> {
    if n_modes != n_generators {
        return Err(PyValueError::new_err(format!(
            "momentum has {n_modes} entries but group has {n_generators} generators"
        )));
    }
    Ok(())
}
