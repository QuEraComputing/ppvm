// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Python bindings for the symmetry-merging primitive.
//!
//! Exposes:
//! - [`TranslationGroup`] PyO3 class with constructors for 1D, 2D, 3D
//!   tori and multi-leg ladders, plus a generic generator-list path.
//! - [`canonicalize_basis_arr`] free function that takes the numpy
//!   `(basis_arr, coeffs)` representation used by `Lindbladian.pc_step_arr`
//!   /  `rk4_step_arr` and merges in place.

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use ppvm_lindblad::{Word, word_from_codes, codes_from_word};
use ppvm_runtime::symmetry as core_sym;
use pyo3::{exceptions::PyValueError, prelude::*};

type PyPauliMap<'py> = (Bound<'py, PyArray2<u8>>, Bound<'py, PyArray1<f64>>);

/// A finite abelian symmetry group acting on qubit positions by
/// permutations. Use this to merge translation-equivalent Pauli strings
/// in either the `Lindbladian.pc_step_arr` basis or the `PauliSum`
/// dictionary, reducing per-step memory by up to `|G|×`.
///
/// Build via the static methods:
/// - `TranslationGroup.chain_1d(n)` — 1D chain of `n` sites with PBC.
/// - `TranslationGroup.torus_2d(lx, ly)` — 2D torus; qubit `(i, j)` at
///   index `j*lx + i`.
/// - `TranslationGroup.torus_3d(lx, ly, lz)` — 3D torus; qubit
///   `(i, j, k)` at index `k*lx*ly + j*lx + i`.
/// - `TranslationGroup.ladder(l, n_legs)` — `n_legs`-leg ladder of `l`
///   sites, translation only along chain direction; qubit `(leg, j)` at
///   index `leg*l + j`.
/// - `TranslationGroup.from_generators(n_qubits, perms, orders)` —
///   arbitrary list of generator permutations + cyclic orders.
#[pyclass(frozen)]
pub struct TranslationGroup {
    inner: core_sym::TranslationGroup,
}

#[pymethods]
impl TranslationGroup {
    #[staticmethod]
    pub fn chain_1d(n: usize) -> Self {
        Self {
            inner: core_sym::TranslationGroup::chain_1d(n),
        }
    }

    #[staticmethod]
    pub fn torus_2d(lx: usize, ly: usize) -> Self {
        Self {
            inner: core_sym::TranslationGroup::torus_2d(lx, ly),
        }
    }

    #[staticmethod]
    pub fn torus_3d(lx: usize, ly: usize, lz: usize) -> Self {
        Self {
            inner: core_sym::TranslationGroup::torus_3d(lx, ly, lz),
        }
    }

    #[staticmethod]
    pub fn ladder(l: usize, n_legs: usize) -> Self {
        Self {
            inner: core_sym::TranslationGroup::ladder(l, n_legs),
        }
    }

    #[staticmethod]
    pub fn from_generators(
        n_qubits: usize,
        perms: Vec<Vec<u32>>,
        orders: Vec<u32>,
    ) -> PyResult<Self> {
        if perms.len() != orders.len() {
            return Err(PyValueError::new_err(format!(
                "perms ({} generators) and orders ({}) must have the same length",
                perms.len(),
                orders.len()
            )));
        }
        for (g, perm) in perms.iter().enumerate() {
            if perm.len() != n_qubits {
                return Err(PyValueError::new_err(format!(
                    "generator {g}: permutation length {} != n_qubits {n_qubits}",
                    perm.len()
                )));
            }
            let mut seen = vec![false; n_qubits];
            for &p in perm {
                let p = p as usize;
                if p >= n_qubits {
                    return Err(PyValueError::new_err(format!(
                        "generator {g}: target {p} out of range [0, {n_qubits})"
                    )));
                }
                if seen[p] {
                    return Err(PyValueError::new_err(format!(
                        "generator {g}: not a permutation (duplicate target {p})"
                    )));
                }
                seen[p] = true;
            }
        }
        Ok(Self {
            inner: core_sym::TranslationGroup::from_generators(n_qubits, perms, orders),
        })
    }

    /// Number of qubits this group acts on.
    #[getter]
    pub fn n_qubits(&self) -> usize {
        self.inner.n_qubits()
    }

    /// Number of generators (rank as an abelian product group).
    #[getter]
    pub fn n_generators(&self) -> usize {
        self.inner.n_generators()
    }

    /// Total group order: product of generator orders.
    #[getter]
    pub fn order(&self) -> usize {
        self.inner.order()
    }

    /// Return the canonical (lex-min) orbit representative of `pauli`.
    /// `pauli` is a length-`n_qubits` uint8 array with the encoding
    /// `I=0, X=1, Y=2, Z=3`. Result is the same shape.
    pub fn canonicalize<'py>(
        &self,
        py: Python<'py>,
        pauli: PyReadonlyArray1<'py, u8>,
    ) -> PyResult<Bound<'py, PyArray1<u8>>> {
        let codes = pauli.as_slice()?;
        if codes.len() != self.inner.n_qubits() {
            return Err(PyValueError::new_err(format!(
                "pauli has length {} but group expects {} qubits",
                codes.len(),
                self.inner.n_qubits()
            )));
        }
        let w = word_from_codes(codes).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let canon = self.inner.canonicalize(&w);
        let mut out = vec![0u8; codes.len()];
        codes_from_word(&canon, &mut out);
        Ok(out.into_pyarray(py))
    }
}

/// Merge a `(basis_arr, coeffs)` Pauli sum (the representation used by
/// `Lindbladian.pc_step_arr` / `rk4_step_arr`) into orbit-representative
/// form. Each row of `basis_arr` is replaced by its canonical
/// representative; coefficients of rows collapsing to the same rep are
/// summed.
///
/// Returns `(merged_basis_arr, merged_coeffs)`. Output length ≤ input
/// length.
///
/// For dynamics that commute with `group` and initial states that are
/// `group`-invariant, this preserves all `group`-invariant expectation
/// values (Theorem 1 of Teng et al., arXiv:2512.12094).
#[pyfunction]
pub fn canonicalize_basis_arr<'py>(
    py: Python<'py>,
    basis: PyReadonlyArray2<'py, u8>,
    coeffs: PyReadonlyArray1<'py, f64>,
    group: &TranslationGroup,
) -> PyResult<PyPauliMap<'py>> {
    let basis_view = basis.as_array();
    let n_q = group.inner.n_qubits();
    if basis_view.shape().get(1).copied() != Some(n_q) {
        return Err(PyValueError::new_err(format!(
            "basis has {} qubits per row but group acts on {n_q}",
            basis_view.shape().get(1).copied().unwrap_or(0)
        )));
    }
    let n = basis_view.shape()[0];
    let coeffs_slice = coeffs.as_slice()?;
    if coeffs_slice.len() != n {
        return Err(PyValueError::new_err(format!(
            "coeffs has length {} but basis has {} rows",
            coeffs_slice.len(),
            n
        )));
    }

    // Decode into Vec<Word>, run the merge, re-encode.
    let mut basis_words: Vec<Word> = Vec::with_capacity(n);
    for i in 0..n {
        let row = basis_view.row(i);
        let row_slice = row.as_slice().ok_or_else(|| {
            PyValueError::new_err("basis array rows are not contiguous; pass a C-order array")
        })?;
        basis_words.push(word_from_codes(row_slice).map_err(|e| PyValueError::new_err(e.to_string()))?);
    }
    let mut coeffs_vec = coeffs_slice.to_vec();

    core_sym::canonicalize_pauli_sum(&mut basis_words, &mut coeffs_vec, &group.inner);

    // Re-encode.
    let m = basis_words.len();
    let mut out_basis = vec![0u8; m * n_q];
    for (i, w) in basis_words.iter().enumerate() {
        codes_from_word(w, &mut out_basis[i * n_q..(i + 1) * n_q]);
    }
    let basis_arr = out_basis
        .into_pyarray(py)
        .reshape([m, n_q])
        .map_err(|e| PyValueError::new_err(format!("reshape failed: {e}")))?;
    Ok((basis_arr, coeffs_vec.into_pyarray(py)))
}
