// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Rust shim for direct Pauli-Lindbladian time evolution on an adaptive
//! Pauli-string basis. See `xy-experiments/main_k_adaptive.py` for the
//! Python driver that owns the linear algebra and basis adaptation.
//!
//! Given a Hermitian Pauli Hamiltonian `H = Σ c_i P_i` and Hermitian Pauli
//! jump operators `L_k` with rates `γ_k ≥ 0`, the adjoint Lindbladian acts
//! on an operator `p` (a single Pauli string) as
//!
//! ```text
//! L*(p) = i Σ_i c_i [P_i, p] + Σ_k γ_k (L_k p L_k − p).
//! ```
//!
//! For Hermitian Pauli `P` and Hermitian Pauli `p`, the product `P · p` has
//! a phase in `{±1, ±i}`. Real phase ⇔ commute (no contribution to the
//! commutator); imaginary phase ⇔ anti-commute, in which case
//! `[P, p] = 2 P p` and `i · 2 P p ∈ {+2 r, −2 r}` where `r = P ⊕ p` (xz
//! bits XORed). Hermitian Pauli jumps give `L p L = ±p`: the anti-commuting
//! case contributes `−2 γ_k · p` (diagonal).
//!
//! The hot path bypasses the higher-level `PauliWord` / `PhasedPauliWord`
//! types because their `set()`/`rehash()` calls dominate when this code is
//! invoked O(10^5) times per evolution step. Each Pauli string is held as
//! a packed pair of bit arrays (`[u8; W]` for X-bits and `[u8; W]` for
//! Z-bits) and we run a single fused multiply that produces both the
//! product phase and the XOR'd output in one byte-level pass.

use dashmap::DashMap;
use fxhash::{FxBuildHasher, FxHashMap};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{exceptions::PyValueError, prelude::*};
use rayon::prelude::*;

/// Bytes per X / Z bit array. 16 bytes = 128 bits, covers up to 128 qubits.
const W: usize = 16;

/// Packed key: first `W` bytes are the X-bit array, next `W` bytes are the
/// Z-bit array. Equality and hashing are over the 32 raw bytes.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct Word {
    x: [u8; W],
    z: [u8; W],
}

impl Word {
    #[inline(always)]
    fn zero() -> Self {
        Self {
            x: [0; W],
            z: [0; W],
        }
    }
}

/// Parsed Hamiltonian or jump term: the packed Pauli word and its real
/// coefficient (or rate, for a jump).
#[derive(Clone)]
struct Term {
    word: Word,
    coeff: f64,
}

/// `(out_basis, out_coeffs)` — packed Pauli rows + their real coefficients.
type PyPauliMap<'py> = (Bound<'py, PyArray2<u8>>, Bound<'py, PyArray1<f64>>);
/// `(rows, cols, vals)` — COO triplets for a real sparse matrix.
type PyCoo<'py> = (
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<u64>>,
    Bound<'py, PyArray1<f64>>,
);

/// Precompiled Lindbladian. Constructed once from string-form Hamiltonian
/// terms + jump operators; reused across many calls to `action`,
/// `leakage`, `generator`.
///
/// `action_cache` memoises the per-input contribution list, the Rust
/// analogue of the `@lru_cache(None)` in the original Python driver.
/// Without it, every `leakage`/`generator` call re-derives the same
/// `L*(p)` from scratch for the same Pauli strings.
#[pyclass]
pub struct LindbladSpec {
    n_qubits: usize,
    h_terms: Vec<Term>,
    j_terms: Vec<Term>,
    /// `h_support[q]` = indices of Hamiltonian terms acting on qubit `q`.
    h_support: Vec<Vec<u32>>,
    /// `j_support[q]` = indices of jump operators acting on qubit `q`.
    j_support: Vec<Vec<u32>>,
    /// `action_cache[p]` = unscaled list of `(output_word, contribution)`
    /// pairs from `L*(p)` excluding the input coefficient. Callers
    /// multiply through by `scale` at the use site.
    action_cache: DashMap<Word, Vec<(Word, f64)>, FxBuildHasher>,
}

// ────────────────── codec helpers ──────────────────

/// Set qubit `q` of `w` to the Pauli encoded by code `b ∈ {0,1,2,3}`.
/// 0 = I, 1 = X, 2 = Z, 3 = Y. (Bit layout: x = b & 1, z = (b >> 1) & 1.)
#[inline(always)]
fn set_code(w: &mut Word, q: usize, b: u8) {
    let byte = q >> 3;
    let mask = 1u8 << (q & 7);
    if b & 1 != 0 {
        w.x[byte] |= mask;
    }
    if b & 2 != 0 {
        w.z[byte] |= mask;
    }
}

/// Inverse of `set_code`: read qubit `q` of `w` and return the Pauli code.
#[inline(always)]
fn get_code(w: &Word, q: usize) -> u8 {
    let byte = q >> 3;
    let mask = 1u8 << (q & 7);
    let xb = ((w.x[byte] & mask) != 0) as u8;
    let zb = ((w.z[byte] & mask) != 0) as u8;
    xb | (zb << 1)
}

/// Decode a row of `n_qubits` Pauli codes into a `Word`.
fn word_from_codes(codes: &[u8], n_qubits: usize) -> PyResult<Word> {
    if codes.len() != n_qubits {
        return Err(PyValueError::new_err(format!(
            "Expected {n_qubits} Pauli codes per row, got {}",
            codes.len()
        )));
    }
    let mut w = Word::zero();
    for (q, &b) in codes.iter().enumerate() {
        if b > 3 {
            return Err(PyValueError::new_err(format!(
                "Pauli code must be 0 (I), 1 (X), 2 (Z), or 3 (Y); got {b}"
            )));
        }
        if b != 0 {
            set_code(&mut w, q, b);
        }
    }
    Ok(w)
}

/// Write `w` into `out` as `n_qubits` Pauli codes (one byte each).
fn codes_from_word(w: &Word, n_qubits: usize, out: &mut [u8]) {
    debug_assert_eq!(out.len(), n_qubits);
    for (q, slot) in out.iter_mut().take(n_qubits).enumerate() {
        *slot = get_code(w, q);
    }
}

/// Parse a string `"IXYZ..."` into a `Word` and the list of qubits where
/// the Pauli is non-identity (the term's support).
fn parse_term(s: &str, n_qubits: usize) -> PyResult<(Word, Vec<u32>)> {
    let chars: Vec<char> = s.chars().filter(|c| *c != '_').collect();
    if chars.len() != n_qubits {
        return Err(PyValueError::new_err(format!(
            "Pauli string \"{s}\" has length {} but n_qubits = {n_qubits}",
            chars.len()
        )));
    }
    if n_qubits > 8 * W {
        return Err(PyValueError::new_err(format!(
            "LindbladSpec supports n_qubits ≤ {}; got {n_qubits}",
            8 * W
        )));
    }
    let mut word = Word::zero();
    let mut support: Vec<u32> = Vec::new();
    for (q, c) in chars.into_iter().enumerate() {
        let code: u8 = match c {
            'I' => 0,
            'X' => 1,
            'Z' => 2,
            'Y' => 3,
            other => {
                return Err(PyValueError::new_err(format!(
                    "Invalid Pauli character '{other}' in \"{s}\""
                )));
            }
        };
        if code != 0 {
            set_code(&mut word, q, code);
            support.push(q as u32);
        }
    }
    Ok((word, support))
}

// ────────────────── algebra helpers ──────────────────

/// `true` if Pauli words `a` and `b` anti-commute.
///
/// Two Pauli strings anti-commute iff
/// `popcount(a.x & b.z) + popcount(a.z & b.x)` is odd.
#[inline(always)]
fn anti_commutes(a: &Word, b: &Word) -> bool {
    let mut bits: u32 = 0;
    for i in 0..W {
        bits += (a.x[i] & b.z[i]).count_ones();
        bits += (a.z[i] & b.x[i]).count_ones();
    }
    bits & 1 == 1
}

/// Fused product `h · p`: returns `(out, eps)` where `out.x = h.x ^ p.x`,
/// `out.z = h.z ^ p.z`, and
///
/// - `eps =  0` if `h` and `p` commute (caller should skip — `[h,p] = 0`),
/// - `eps = -2.0` if `h·p` has phase `+i` (so `i·[h,p] = -2·out`),
/// - `eps = +2.0` if `h·p` has phase `-i` (so `i·[h,p] = +2·out`).
///
/// Phase = `(2·sign_count + imag_count) mod 4`, where the per-byte XOR/AND
/// formulas are the same ones used by `PhasedPauliWord::mul_assign`.
#[inline(always)]
fn comm_product(h: &Word, p: &Word) -> (Word, f64) {
    let mut out = Word::zero();
    let mut sign_count: u32 = 0;
    let mut imag_count: u32 = 0;
    for i in 0..W {
        let a = h.x[i];
        let b = h.z[i];
        let c = p.x[i];
        let d = p.z[i];
        let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
        let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
        sign_count += sign.count_ones();
        imag_count += imag.count_ones();
        out.x[i] = a ^ c;
        out.z[i] = b ^ d;
    }
    let phase = (2 * sign_count + imag_count) & 3;
    let eps = match phase {
        1 => -2.0,
        3 => 2.0,
        _ => 0.0,
    };
    (out, eps)
}

/// Union of `index[q]` for each `q ∈ p_support`, deduped.
#[inline]
fn candidate_terms(p_support: &[u32], index: &[Vec<u32>], scratch: &mut Vec<u32>) {
    scratch.clear();
    for &q in p_support {
        scratch.extend_from_slice(&index[q as usize]);
    }
    scratch.sort_unstable();
    scratch.dedup();
}

/// Find qubits where `p` is non-identity (the support of `p`).
#[inline]
fn word_support(p: &Word, n_qubits: usize, out: &mut Vec<u32>) {
    out.clear();
    for q in 0..n_qubits {
        let byte = q >> 3;
        let mask = 1u8 << (q & 7);
        if (p.x[byte] | p.z[byte]) & mask != 0 {
            out.push(q as u32);
        }
    }
}

/// Compute the unscaled list of `(output, coefficient)` pairs that
/// `L*(p)` would contribute, without the input coefficient. Used by
/// `cached_action` to populate the cache. The scratch buffers are
/// allocated by the caller so this is allocation-light per p.
fn compute_action_terms(
    spec: &LindbladSpec,
    p: &Word,
    scratch_support: &mut Vec<u32>,
    scratch_cands: &mut Vec<u32>,
) -> Vec<(Word, f64)> {
    word_support(p, spec.n_qubits, scratch_support);
    // Aggregate within a small local map; for typical action sizes the
    // overhead of FxHashMap is dwarfed by the dedup it provides.
    let mut local: FxHashMap<Word, f64> = FxHashMap::default();

    // Hamiltonian (off-diagonal) terms
    candidate_terms(scratch_support, &spec.h_support, scratch_cands);
    for &i in scratch_cands.iter() {
        let h = &spec.h_terms[i as usize];
        let (r, eps) = comm_product(&h.word, p);
        if eps != 0.0 {
            *local.entry(r).or_insert(0.0) += h.coeff * eps;
        }
    }

    // Jump (diagonal) terms — single diagonal entry on `p`.
    candidate_terms(scratch_support, &spec.j_support, scratch_cands);
    let mut diag: f64 = 0.0;
    for &k in scratch_cands.iter() {
        let j = &spec.j_terms[k as usize];
        if anti_commutes(&j.word, p) {
            diag += -2.0 * j.coeff;
        }
    }
    if diag != 0.0 {
        *local.entry(*p).or_insert(0.0) += diag;
    }

    local.into_iter().collect()
}

/// Accumulate `scale · L*(p)` into `out`, hitting the cache when warm.
fn accumulate_action(
    spec: &LindbladSpec,
    p: &Word,
    scale: f64,
    out: &mut FxHashMap<Word, f64>,
    scratch_support: &mut Vec<u32>,
    scratch_cands: &mut Vec<u32>,
) {
    if let Some(r) = spec.action_cache.get(p) {
        for (w, c) in r.value().iter() {
            *out.entry(*w).or_insert(0.0) += scale * c;
        }
        return;
    }
    // Compute, install in cache, and accumulate.
    let terms = compute_action_terms(spec, p, scratch_support, scratch_cands);
    for (w, c) in terms.iter() {
        *out.entry(*w).or_insert(0.0) += scale * c;
    }
    // If a concurrent task raced and inserted first, the recomputed
    // value is identical so the overwrite is harmless.
    spec.action_cache.insert(*p, terms);
}

/// Decode a `(N, n_qubits)` uint8 ndarray view into `N` packed `Word`s.
fn decode_basis(view: &numpy::ndarray::ArrayView2<u8>, n_qubits: usize) -> PyResult<Vec<Word>> {
    let n_basis = view.shape()[0];
    let mut out: Vec<Word> = Vec::with_capacity(n_basis);
    for i in 0..n_basis {
        let row = view.row(i);
        let mut w = Word::zero();
        for q in 0..n_qubits {
            let b = row[q];
            if b > 3 {
                return Err(PyValueError::new_err(format!(
                    "Pauli code must be 0-3; got {b} at row {i}, col {q}"
                )));
            }
            if b != 0 {
                set_code(&mut w, q, b);
            }
        }
        out.push(w);
    }
    Ok(out)
}

// ────────────────── PyClass ──────────────────

#[pymethods]
impl LindbladSpec {
    #[new]
    #[pyo3(signature = (n_qubits, h_terms, h_coeffs, jump_terms, jump_rates))]
    fn new(
        n_qubits: usize,
        h_terms: Vec<String>,
        h_coeffs: Vec<f64>,
        jump_terms: Vec<String>,
        jump_rates: Vec<f64>,
    ) -> PyResult<Self> {
        if n_qubits > 8 * W {
            return Err(PyValueError::new_err(format!(
                "LindbladSpec supports n_qubits ≤ {}; got {n_qubits}",
                8 * W
            )));
        }
        if h_terms.len() != h_coeffs.len() {
            return Err(PyValueError::new_err(
                "h_terms and h_coeffs must have the same length",
            ));
        }
        if jump_terms.len() != jump_rates.len() {
            return Err(PyValueError::new_err(
                "jump_terms and jump_rates must have the same length",
            ));
        }

        let mut h: Vec<Term> = Vec::with_capacity(h_terms.len());
        let mut h_support_idx: Vec<Vec<u32>> = vec![Vec::new(); n_qubits];
        for (i, (s, c)) in h_terms.iter().zip(h_coeffs.iter()).enumerate() {
            let (word, support) = parse_term(s, n_qubits)?;
            for q in support {
                h_support_idx[q as usize].push(i as u32);
            }
            h.push(Term { word, coeff: *c });
        }

        let mut j: Vec<Term> = Vec::with_capacity(jump_terms.len());
        let mut j_support_idx: Vec<Vec<u32>> = vec![Vec::new(); n_qubits];
        for (k, (s, g)) in jump_terms.iter().zip(jump_rates.iter()).enumerate() {
            let (word, support) = parse_term(s, n_qubits)?;
            for q in support {
                j_support_idx[q as usize].push(k as u32);
            }
            j.push(Term { word, coeff: *g });
        }

        Ok(Self {
            n_qubits,
            h_terms: h,
            j_terms: j,
            h_support: h_support_idx,
            j_support: j_support_idx,
            action_cache: DashMap::with_hasher(FxBuildHasher::default()),
        })
    }

    /// Drop all memoised `L*(p)` entries. Useful between independent
    /// problems sharing the same spec but acting on disjoint Pauli sets.
    fn clear_cache(&self) {
        self.action_cache.clear();
    }

    #[getter]
    fn cache_size(&self) -> usize {
        self.action_cache.len()
    }

    #[getter]
    fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    #[getter]
    fn num_h_terms(&self) -> usize {
        self.h_terms.len()
    }

    #[getter]
    fn num_jump_terms(&self) -> usize {
        self.j_terms.len()
    }

    /// Apply `L*` to a single Pauli string `p`.
    fn action<'py>(
        &self,
        py: Python<'py>,
        p: PyReadonlyArray1<'py, u8>,
    ) -> PyResult<PyPauliMap<'py>> {
        let p_slice = p.as_slice()?;
        let p_word = word_from_codes(p_slice, self.n_qubits)?;
        let mut out: FxHashMap<Word, f64> = FxHashMap::default();
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        accumulate_action(self, &p_word, 1.0, &mut out, &mut s1, &mut s2);

        let m = out.len();
        let n = self.n_qubits;
        let mut basis = vec![0u8; m * n];
        let mut coeffs = vec![0f64; m];
        for (i, (w, c)) in out.into_iter().enumerate() {
            codes_from_word(&w, n, &mut basis[i * n..(i + 1) * n]);
            coeffs[i] = c;
        }
        let basis_arr = basis
            .into_pyarray(py)
            .reshape([m, n])
            .map_err(|e| PyValueError::new_err(format!("reshape failed: {e}")))?;
        Ok((basis_arr, coeffs.into_pyarray(py)))
    }

    /// Off-basis component of `L*( Σ_j coeffs[j] · basis[j] )`.
    ///
    /// `basis` is `(N, n_qubits)` uint8 Pauli codes; `coeffs` is length-N
    /// float64. Returns `(out_basis, out_coeffs)` with output Paulis NOT
    /// in `basis` (and not in `protected`, if given).
    #[pyo3(signature = (basis, coeffs, protected = None))]
    fn leakage<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
        coeffs: PyReadonlyArray1<'py, f64>,
        protected: Option<PyReadonlyArray2<'py, u8>>,
    ) -> PyResult<PyPauliMap<'py>> {
        let basis_view = basis.as_array();
        let coeffs_view = coeffs.as_slice()?;
        let (n_basis, n_q) = (basis_view.shape()[0], basis_view.shape()[1]);
        if n_q != self.n_qubits {
            return Err(PyValueError::new_err(format!(
                "basis has {n_q} columns but spec.n_qubits = {}",
                self.n_qubits
            )));
        }
        if coeffs_view.len() != n_basis {
            return Err(PyValueError::new_err(format!(
                "coeffs has length {} but basis has {n_basis} rows",
                coeffs_view.len()
            )));
        }

        let basis_words = decode_basis(&basis_view, self.n_qubits)?;
        let in_basis: FxHashMap<Word, ()> = basis_words.iter().map(|w| (*w, ())).collect();

        let protected_set: FxHashMap<Word, ()> = if let Some(ref prot) = protected {
            let pv = prot.as_array();
            if pv.shape()[1] != self.n_qubits {
                return Err(PyValueError::new_err(format!(
                    "protected has {} columns but spec.n_qubits = {}",
                    pv.shape()[1],
                    self.n_qubits
                )));
            }
            let words = decode_basis(&pv, self.n_qubits)?;
            words.into_iter().map(|w| (w, ())).collect()
        } else {
            FxHashMap::default()
        };

        let local: Vec<FxHashMap<Word, f64>> = basis_words
            .par_iter()
            .zip(coeffs_view.par_iter())
            .map(|(p, &c)| {
                let mut m: FxHashMap<Word, f64> = FxHashMap::default();
                let mut s1 = Vec::new();
                let mut s2 = Vec::new();
                accumulate_action(self, p, c, &mut m, &mut s1, &mut s2);
                m
            })
            .collect();

        let mut merged: FxHashMap<Word, f64> = FxHashMap::default();
        for m in local {
            for (k, v) in m {
                if in_basis.contains_key(&k) || protected_set.contains_key(&k) {
                    continue;
                }
                *merged.entry(k).or_insert(0.0) += v;
            }
        }

        let m = merged.len();
        let n = self.n_qubits;
        let mut out_basis = vec![0u8; m * n];
        let mut out_coeffs = vec![0f64; m];
        for (i, (w, c)) in merged.into_iter().enumerate() {
            codes_from_word(&w, n, &mut out_basis[i * n..(i + 1) * n]);
            out_coeffs[i] = c;
        }
        let basis_arr = out_basis
            .into_pyarray(py)
            .reshape([m, n])
            .map_err(|e| PyValueError::new_err(format!("reshape failed: {e}")))?;
        Ok((basis_arr, out_coeffs.into_pyarray(py)))
    }

    /// Sparse generator matrix in COO form: returns `(rows, cols, vals)`
    /// that the caller wraps with scipy.
    ///
    /// Row index = output Pauli's position in `basis`; col index = input
    /// Pauli's position. Output Paulis not in `basis` are silently
    /// dropped (handled by `leakage`).
    fn generator<'py>(
        &self,
        py: Python<'py>,
        basis: PyReadonlyArray2<'py, u8>,
    ) -> PyResult<PyCoo<'py>> {
        let basis_view = basis.as_array();
        let n_q = basis_view.shape()[1];
        if n_q != self.n_qubits {
            return Err(PyValueError::new_err(format!(
                "basis has {n_q} columns but spec.n_qubits = {}",
                self.n_qubits
            )));
        }

        let basis_words = decode_basis(&basis_view, self.n_qubits)?;
        let index: FxHashMap<Word, u32> = basis_words
            .iter()
            .enumerate()
            .map(|(i, w)| (*w, i as u32))
            .collect();

        let local: Vec<Vec<(u64, u64, f64)>> = basis_words
            .par_iter()
            .enumerate()
            .map(|(col, p)| {
                let mut m: FxHashMap<Word, f64> = FxHashMap::default();
                let mut s1 = Vec::new();
                let mut s2 = Vec::new();
                accumulate_action(self, p, 1.0, &mut m, &mut s1, &mut s2);
                let mut out = Vec::with_capacity(m.len());
                for (w, v) in m {
                    if let Some(&row) = index.get(&w) {
                        out.push((row as u64, col as u64, v));
                    }
                }
                out
            })
            .collect();

        let total: usize = local.iter().map(|v| v.len()).sum();
        let mut rows = Vec::with_capacity(total);
        let mut cols = Vec::with_capacity(total);
        let mut vals = Vec::with_capacity(total);
        for trips in local {
            for (r, c, v) in trips {
                rows.push(r);
                cols.push(c);
                vals.push(v);
            }
        }

        Ok((
            rows.into_pyarray(py),
            cols.into_pyarray(py),
            vals.into_pyarray(py),
        ))
    }
}
