// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-step orbit-representative evolution under translation symmetry.
//!
//! The state lives entirely in **orbit-rep form** throughout: `basis`
//! contains only canonical translation-orbit representatives, and
//! `coeffs` are complex (one per rep). The dynamics `L*` is computed
//! with **phase-aware action** — for each output Pauli `q`, we
//! canonicalize `q` to its orbit rep `r_q` with shift counter `cnt_q`,
//! and accumulate `χ_k(g_{cnt_q}) · v · c_r` (where `v` is the matrix
//! element of `L*` between input rep `r` and output `q`).
//!
//! This gives the user the per-step memory reduction promised by the
//! symmetry-merging paper: basis is ~|G|× smaller than the full-basis
//! representation, throughout the entire evolution.
//!
//! The phase-aware action is genuinely **complex** (because of the
//! `χ_k(g)` phase factors). Rather than materialise a CSR, the per-column
//! action — for each input rep, the list of `(row, χ_k·v)` pairs for the
//! in-basis outputs — is computed **once per expm call** (via
//! [`build_orbit_rep_cols`]) and then reused, CSC-style, across every
//! Krylov–Taylor matvec driving the external `quspin-expm` engine. No CSR
//! is ever materialised on the production path.
//!
//! ## Limitations
//!
//! - Caller is responsible for ensuring the input basis is in orbit-rep
//!   form (i.e. each entry is the canonical representative of its
//!   translation orbit). Use [`canonicalize_basis_to_rep`] if needed.
//! - The momentum sector `k_modes` is fixed for the duration of one
//!   pc_step call. To compute a full site-resolved profile, call
//!   `pc_step_orbit_rep` once per momentum mode and inverse-Fourier
//!   the results.

use crate::{Error, LindbladSpec};
use fxhash::{FxBuildHasher, FxHashMap};
use num::Complex;
use ppvm_pauli_sum::symmetry::TranslationGroup;
use quspin_expm::ExpmOp;
use quspin_types::{ExpmComputation, LinearOperator, QuSpinError};
use rayon::prelude::*;

// Word type re-exported from lib.rs.
use crate::Word;

/// Replace each entry of `basis` with its canonical orbit
/// representative under `group`. Pure rewrite; coefficients are
/// untouched. Useful to enforce the orbit-rep invariant before calling
/// [`pc_step_orbit_rep`].
///
/// Does NOT deduplicate — if multiple input entries collapse to the
/// same rep, both are kept (caller should run a merge afterwards).
pub fn canonicalize_basis_to_rep(basis: &mut [Word], group: &TranslationGroup) {
    for w in basis.iter_mut() {
        let canon = group.canonicalize(w);
        *w = canon;
    }
}

/// Build the per-column phase-aware action of the in-basis-restricted
/// orbit-rep generator `M` at momentum sector `k_modes`.
///
/// Returns, for each input rep `basis[c]` (column `c`), the list of
/// `(row, χ_k(g_{cnt_q}) · v_q)` pairs for every action output Pauli `q`
/// of `L*(basis[c])` whose orbit rep `r_q` is in `basis` at index `row`.
/// Outputs not in `basis` are dropped. This is the expensive part of the
/// orbit-rep dynamics (`compute_action_terms` + `canonicalize_with_shift`
/// + `character`); it is computed once and reused by the CSC-style matvec
/// in [`OrbitRepCscOp`] (which reads it directly).
pub(crate) fn build_orbit_rep_cols(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
    group: &TranslationGroup,
    k_modes: &[i32],
) -> Vec<Vec<(u32, Complex<f64>)>> {
    basis
        .par_iter()
        .map_init(
            || {
                (
                    Vec::<u32>::with_capacity(spec.n_qubits()),
                    Vec::<u32>::with_capacity(128),
                    FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                        128,
                        FxBuildHasher::default(),
                    ),
                )
            },
            |(s1, s2, lm), r| {
                let terms = spec.compute_action_terms(r, s1, s2, lm);
                let mut out = Vec::with_capacity(terms.len());
                for (q, v) in terms.iter() {
                    let (r_q, cnt_q) = group.canonicalize_with_shift(q);
                    if let Some(&row) = index.get(&r_q) {
                        let phase = group.character(k_modes, &cnt_q);
                        out.push((row, phase * (*v as f64)));
                    }
                }
                out
            },
        )
        .collect()
}

/// Borrowed CSC-style view of the in-basis-restricted orbit-rep generator
/// `M`, backed by the cached per-column action [`build_orbit_rep_cols`]
/// (computed once per expm call). The momentum-character phases make the
/// entries complex, so this implements `LinearOperator<Complex<f64>>`.
///
/// `cols[c]` holds `(row, χ_k·v)` pairs for column `c`; `dot` does a CSC
/// matvec `y = M·x` against the cached action, with no per-matvec action
/// recompute and no CSR materialisation.
pub(crate) struct OrbitRepCscOp<'a> {
    cols: &'a [Vec<(u32, Complex<f64>)>],
    dim: usize,
}

impl LinearOperator<Complex<f64>> for OrbitRepCscOp<'_> {
    fn dim(&self) -> usize {
        self.dim
    }

    fn parallel_hint(&self) -> bool {
        // `dot` parallelises internally over column chunks, and we drive
        // the sequential single-vector `apply` path; never let quspin run
        // its persistent-thread pool on top of our rayon parallelism.
        false
    }

    fn dot(
        &self,
        overwrite: bool,
        input: &[Complex<f64>],
        output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        let n = self.dim;
        let zero = Complex::new(0.0, 0.0);
        if n == 0 {
            return Ok(());
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);

        // Parallelise over column chunks: each thread accumulates into a
        // dense local `y` of length `dim`, reading the cached action.
        let partial_ys: Vec<Vec<Complex<f64>>> = self
            .cols
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![zero; n];
                for (c_local, col) in chunk.iter().enumerate() {
                    let c = c_offset + c_local;
                    let xc = input[c];
                    if xc == zero {
                        continue;
                    }
                    for &(row, val) in col.iter() {
                        y_local[row as usize] += val * xc;
                    }
                }
                y_local
            })
            .collect();

        if overwrite {
            output.fill(zero);
        }
        for partial in &partial_ys {
            for (oi, &pi) in output.iter_mut().zip(partial.iter()) {
                *oi += pi;
            }
        }
        Ok(())
    }

    fn trace(&self) -> Complex<f64> {
        // Computed eagerly in `expm_apply_orbit_rep_cached`; never reached
        // on the `from_parts` + single-vector `apply` path.
        unreachable!("OrbitRepCscOp::trace not used on the from_parts apply path")
    }

    fn onenorm(&self, _shift: Complex<f64>) -> f64 {
        unreachable!("OrbitRepCscOp::onenorm not used on the from_parts apply path")
    }

    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "OrbitRepCscOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }

    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, Complex<f64>>,
        _output: ndarray::ArrayViewMut2<'_, Complex<f64>>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "OrbitRepCscOp: dot_many not used on the from_parts apply path".into(),
        ))
    }

    fn dot_chunk(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output_chunk: &mut [Complex<f64>],
        _row_start: usize,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "OrbitRepCscOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }

    fn dot_transpose_chunk(
        &self,
        _input: &[Complex<f64>],
        _output: &[<Complex<f64> as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "OrbitRepCscOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Matrix-free (streaming) orbit-rep operator: recomputes the phase-aware
/// action per matvec instead of caching it (the `Complex` mirror of
/// [`crate::mf_expm`]'s streaming path). Holds NO `nnz`-sized cache — used for
/// the opt-in low-RAM mode (`PPVM_EXPM_STREAM`), at the cost of recomputing the
/// `compute_action_terms` + canonicalize + character work on every matvec.
pub(crate) struct StreamOrbitOp<'a> {
    spec: &'a LindbladSpec,
    basis: &'a [Word],
    index: &'a FxHashMap<Word, u32>,
    group: &'a TranslationGroup,
    k_modes: &'a [i32],
    dim: usize,
}

impl LinearOperator<Complex<f64>> for StreamOrbitOp<'_> {
    fn dim(&self) -> usize {
        self.dim
    }
    fn parallel_hint(&self) -> bool {
        false
    }
    fn dot(
        &self,
        overwrite: bool,
        input: &[Complex<f64>],
        output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        let n = self.dim;
        let zero = Complex::new(0.0, 0.0);
        if n == 0 {
            return Ok(());
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);
        let partial_ys: Vec<Vec<Complex<f64>>> = self
            .basis
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![zero; n];
                let mut s1 = Vec::<u32>::with_capacity(self.spec.n_qubits());
                let mut s2 = Vec::<u32>::with_capacity(128);
                let mut lm = FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                    128,
                    FxBuildHasher::default(),
                );
                for (c_local, r) in chunk.iter().enumerate() {
                    let xc = input[c_offset + c_local];
                    if xc == zero {
                        continue;
                    }
                    let terms = self.spec.compute_action_terms(r, &mut s1, &mut s2, &mut lm);
                    for (q, v) in terms.iter() {
                        let (r_q, cnt_q) = self.group.canonicalize_with_shift(q);
                        if let Some(&row) = self.index.get(&r_q) {
                            let phase = self.group.character(self.k_modes, &cnt_q);
                            y_local[row as usize] += phase * (*v as f64) * xc;
                        }
                    }
                }
                y_local
            })
            .collect();
        if overwrite {
            output.fill(zero);
        }
        for partial in &partial_ys {
            for (oi, &pi) in output.iter_mut().zip(partial.iter()) {
                *oi += pi;
            }
        }
        Ok(())
    }
    fn trace(&self) -> Complex<f64> {
        unreachable!("StreamOrbitOp::trace not used on the from_parts apply path")
    }
    fn onenorm(&self, _shift: Complex<f64>) -> f64 {
        unreachable!("StreamOrbitOp::onenorm not used on the from_parts apply path")
    }
    fn dot_transpose(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output: &mut [Complex<f64>],
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "StreamOrbitOp: dot_transpose not used on the from_parts apply path".into(),
        ))
    }
    fn dot_many(
        &self,
        _overwrite: bool,
        _input: ndarray::ArrayView2<'_, Complex<f64>>,
        _output: ndarray::ArrayViewMut2<'_, Complex<f64>>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "StreamOrbitOp: dot_many not used on the from_parts apply path".into(),
        ))
    }
    fn dot_chunk(
        &self,
        _overwrite: bool,
        _input: &[Complex<f64>],
        _output_chunk: &mut [Complex<f64>],
        _row_start: usize,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "StreamOrbitOp: dot_chunk not used on the from_parts apply path".into(),
        ))
    }
    fn dot_transpose_chunk(
        &self,
        _input: &[Complex<f64>],
        _output: &[<Complex<f64> as ExpmComputation>::Atomic],
        _rows: std::ops::Range<usize>,
    ) -> Result<(), QuSpinError> {
        Err(QuSpinError::RuntimeError(
            "StreamOrbitOp: dot_transpose_chunk not used on the from_parts apply path".into(),
        ))
    }
}

/// Light phase-aware pass computing ONLY `per_col[c] = (raw, diag)` for the
/// `μ`/1-norm selection, without materialising the cached columns. Used by the
/// streaming path so the norm pass costs no `nnz`-sized memory.
fn per_col_orbit_stream(
    spec: &LindbladSpec,
    basis: &[Word],
    index: &FxHashMap<Word, u32>,
    group: &TranslationGroup,
    k_modes: &[i32],
) -> Vec<(f64, Complex<f64>)> {
    basis
        .par_iter()
        .enumerate()
        .map_init(
            || {
                (
                    Vec::<u32>::with_capacity(spec.n_qubits()),
                    Vec::<u32>::with_capacity(128),
                    FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                        128,
                        FxBuildHasher::default(),
                    ),
                )
            },
            |(s1, s2, lm), (c, r)| {
                let terms = spec.compute_action_terms(r, s1, s2, lm);
                let mut raw = 0.0_f64;
                let mut diag = Complex::new(0.0, 0.0);
                for (q, v) in terms.iter() {
                    let (r_q, cnt_q) = group.canonicalize_with_shift(q);
                    if let Some(&row) = index.get(&r_q) {
                        let val = group.character(k_modes, &cnt_q) * (*v as f64);
                        raw += val.norm();
                        if row as usize == c {
                            diag += val;
                        }
                    }
                }
                (raw, diag)
            },
        )
        .collect()
}

/// Compute `exp(dt · M) · coeffs` for the in-basis-restricted orbit-rep
/// generator `M` at momentum sector `k_modes`, via `quspin-expm`. Returns
/// a fresh `Vec<Complex<f64>>` of length `basis.len()`. No CSR is
/// materialised.
///
/// The expensive phase-aware action is computed ONCE here (via
/// [`build_orbit_rep_cols`]) and reused, CSC-style, across every Krylov–
/// Taylor matvec (see [`OrbitRepCscOp`]). One pass over the cached columns
/// extracts the diagonal shift `μ = tr(M)/n` and a valid upper bound on
/// the column 1-norm of `M − μ·I`; from `‖dt·(M−μI)‖₁` we pick the Taylor
/// partition `(m*, s)` and hand everything to
/// [`quspin_expm::ExpmOp::from_parts`] (mirroring
/// [`crate::mf_expm::expm_apply_mf`]).
pub(crate) fn expm_apply_orbit_rep_cached(
    spec: &LindbladSpec,
    basis: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
    dt: f64,
    coeffs: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    let n = basis.len();
    if n == 0 {
        return Vec::new();
    }

    let index = crate::build_basis_index(basis);
    // Opt-in streaming (low-RAM) path: skip the phase-aware column cache and
    // recompute the action per matvec (StreamOrbitOp). Removes the dominant
    // peak-memory term (the nnz-sized complex cache on the doubly-enriched
    // corrector basis) at the cost of wall. Default (unset) = cached
    // OrbitRepCscOp (fast).
    let stream = std::env::var_os("PPVM_EXPM_STREAM").is_some();
    let cols_opt = if stream {
        None
    } else {
        Some(build_orbit_rep_cols(spec, basis, &index, group, k_modes))
    };

    // One pass for the per-column `(raw, diag)` used by the `μ`/1-norm
    // selection: `raw = Σ|val|` (upper bound on the absolute column sum),
    // `diag = M[c,c]`. From these: `trace = Σ diag`, `μ = trace/n`, and an
    // upper bound on the column 1-norm of `M − μ·I`: `raw − |diag| + |diag − μ|`.
    let per_col: Vec<(f64, Complex<f64>)> = match &cols_opt {
        Some(cols) => cols
            .par_iter()
            .enumerate()
            .map(|(c, col)| {
                let mut raw = 0.0_f64;
                let mut diag = Complex::new(0.0, 0.0);
                for &(row, val) in col.iter() {
                    raw += val.norm();
                    if row as usize == c {
                        diag += val;
                    }
                }
                (raw, diag)
            })
            .collect(),
        None => per_col_orbit_stream(spec, basis, &index, group, k_modes),
    };

    let trace: Complex<f64> = per_col.iter().map(|(_, d)| *d).sum();
    let mu = trace / n as f64;
    let onenorm = per_col
        .iter()
        .map(|(raw, diag)| raw - diag.norm() + (diag - mu).norm())
        .fold(0.0_f64, f64::max);

    let (m_star, s) = crate::expm::select_ms(dt.abs() * onenorm, None);

    let mut v = coeffs.to_vec();
    match &cols_opt {
        Some(cols) => {
            let op = OrbitRepCscOp { cols, dim: n };
            let e = ExpmOp::from_parts(
                op,
                Complex::new(dt, 0.0),
                mu,
                s as usize,
                m_star as usize,
                1e-12_f64,
            );
            e.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
                .expect("expm apply (orbit-rep cached)");
        }
        None => {
            let op = StreamOrbitOp {
                spec,
                basis,
                index: &index,
                group,
                k_modes,
                dim: n,
            };
            let e = ExpmOp::from_parts(
                op,
                Complex::new(dt, 0.0),
                mu,
                s as usize,
                m_star as usize,
                1e-12_f64,
            );
            e.apply(ndarray::ArrayViewMut1::from(v.as_mut_slice()))
                .expect("expm apply (orbit-rep stream)");
        }
    }
    v
}

/// Phase-aware leakage: out-of-basis component of `L*(O_k)` where `O_k`
/// is the operator represented by `basis` (orbit reps) and `coeffs`
/// (complex coefficients in momentum sector `k_modes`).
///
/// For each input rep `r` with coefficient `c_r`, and each output `q`
/// of `L*(r) = Σ_q v_q · q`:
/// 1. Canonicalize `q` → `(r_q, cnt_q)`.
/// 2. If `r_q` NOT in `basis` and NOT in `protected`:
///    `merged[r_q] += χ_k(g_{cnt_q}) · v_q · c_r`.
///
/// Returns `(r_q, sum)` pairs for all candidates with nonzero sum.
///
/// The live candidate map is capped to the *available room*
/// `room = max_basis − basis.len()` (the reps we could actually add),
/// applied during accumulation: input reps are processed in descending
/// `|c|` order and after each chunk only the `room` largest-`|sum|`
/// candidates are kept. A large `max_basis` (room ≥ all candidates)
/// disables the cap — the near-exact case.
#[allow(clippy::too_many_arguments)]
pub fn leakage_orbit_rep(
    spec: &LindbladSpec,
    basis: &[Word],
    coeffs: &[Complex<f64>],
    protected: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
    max_basis: usize,
) -> Result<Vec<(Word, Complex<f64>)>, Error> {
    if basis.len() != coeffs.len() {
        return Err(Error::LengthMismatch {
            what: "basis and coeffs",
            a: basis.len(),
            b: coeffs.len(),
        });
    }
    let in_basis: FxHashMap<&Word, ()> = basis.iter().map(|w| (w, ())).collect();
    let protected_set: FxHashMap<&Word, ()> = protected.iter().map(|w| (w, ())).collect();

    // Descending sort by |c|: process largest-magnitude contributors first
    // so the running room-cap keeps the right entries.
    let mut order: Vec<usize> = (0..basis.len()).collect();
    order.sort_by(|&a, &b| {
        coeffs[b]
            .norm()
            .partial_cmp(&coeffs[a].norm())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    const CHUNK_SIZE: usize = 4096;
    let room = max_basis.saturating_sub(basis.len());
    let mut merged: FxHashMap<Word, Complex<f64>> = FxHashMap::default();
    for chunk_indices in order.chunks(CHUNK_SIZE) {
        let local: Vec<Vec<(Word, Complex<f64>)>> = chunk_indices
            .par_iter()
            .map_init(
                || {
                    (
                        Vec::<u32>::with_capacity(spec.n_qubits()),
                        Vec::<u32>::with_capacity(128),
                        FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                            128,
                            FxBuildHasher::default(),
                        ),
                    )
                },
                |(s1, s2, lm), &i| {
                    let r = &basis[i];
                    let c_r = coeffs[i];
                    let terms = spec.compute_action_terms(r, s1, s2, lm);
                    let mut out = Vec::with_capacity(terms.len());
                    for (q, v) in terms.iter() {
                        let (r_q, cnt_q) = group.canonicalize_with_shift(q);
                        if !in_basis.contains_key(&r_q) && !protected_set.contains_key(&r_q) {
                            let phase = group.character(k_modes, &cnt_q);
                            out.push((r_q, phase * (*v as f64) * c_r));
                        }
                    }
                    out
                },
            )
            .collect();
        for v in local {
            for (k, val) in v {
                *merged.entry(k).or_insert(Complex::new(0.0, 0.0)) += val;
            }
        }

        // Room-cap: keep only the `room` largest-magnitude entries.
        if merged.len() > room {
            if room == 0 {
                merged.clear();
            } else {
                let mut mags: Vec<f64> = merged.values().map(|v| v.norm()).collect();
                let k = room.min(mags.len() - 1);
                mags.select_nth_unstable_by(k, |a, b| {
                    b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
                });
                let cutoff = mags[k];
                merged.retain(|_, &mut v| v.norm() >= cutoff);
            }
        }
    }
    Ok(merged.into_iter().filter(|(_, c)| c.norm() > 0.0).collect())
}

/// Per-step orbit-rep predictor-corrector evolution.
///
/// All state lives in orbit-rep form throughout. Each pc step does:
/// 1. Phase-aware leakage from `(basis, coeffs)`; append the largest
///    leakage reps, up to `room = max_basis − basis.len()`.
/// 2. Predictor: build complex orbit-rep CSR, run `expm_multiply_cx`.
/// 3. Phase-aware leakage from predicted state; append further reps.
/// 4. Corrector: rebuild CSR (basis grew), `expm_multiply_cx` from
///    pre-step coefficients (with the same trick as `pc_step_inner`).
/// 5. Prune `|c| < drop_tol`, then trim to the top-`max_basis` reps by
///    `|c|`; protected reps never dropped.
///
/// `max_basis` is a hard rank cap on the live orbit-rep basis: enrichment
/// adds at most `max_basis − basis.len()` of the largest leakage reps, the
/// leakage map is capped to the same room, and the post-step basis is
/// trimmed to the top-`max_basis` by `|c|`. Pass a large value (e.g.
/// `usize::MAX`) for the near-exact, uncapped case. `drop_tol` additionally
/// prunes by magnitude.
///
/// `basis` is assumed to contain only canonical orbit representatives.
/// If not, [`canonicalize_basis_to_rep`] should be called first.
#[allow(clippy::too_many_arguments)]
pub fn pc_step_orbit_rep(
    spec: &LindbladSpec,
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    dt: f64,
    max_basis: usize,
    drop_tol: f64,
    protected: &[Word],
    group: &TranslationGroup,
    k_modes: &[i32],
) -> Result<(), Error> {
    // Rate-based admission filter, same semantics as the real-space path:
    // a leakage rep is admitted only if its rate exceeds
    // `tau_add = K * drop_tol / dt` (K from PPVM_K_LEAKAGE, default 0 =>
    // admit everything, the historical behaviour).
    let tau_add = if dt > 0.0 { crate::k_leakage() * drop_tol / dt } else { 0.0 };
    // 1. First-hop phase-aware leakage.
    let mut leak = leakage_orbit_rep(spec, basis, coeffs, protected, group, k_modes, max_basis)?;
    if tau_add > 0.0 {
        leak.retain(|(_, c)| c.norm() > tau_add);
    }
    add_leakage_capped_complex(basis, coeffs, leak, max_basis);
    // 2. Predictor: cache-the-action expm (no CSR materialised; the
    //    phase-aware action is built once via `build_orbit_rep_cols`).
    let coeffs_predict = expm_apply_orbit_rep_cached(spec, basis, group, k_modes, dt, coeffs);
    // 3. Second-hop leakage from predicted state.
    let mut leak2 =
        leakage_orbit_rep(spec, basis, &coeffs_predict, protected, group, k_modes, max_basis)?;
    drop(coeffs_predict);
    if tau_add > 0.0 {
        leak2.retain(|(_, c)| c.norm() > tau_add);
    }
    add_leakage_capped_complex(basis, coeffs, leak2, max_basis);
    // 4. Corrector: cache-the-action expm from pre-step state (basis grew).
    *coeffs = expm_apply_orbit_rep_cached(spec, basis, group, k_modes, dt, coeffs);
    // 5. Prune by magnitude, then rank-cap to max_basis.
    if drop_tol > 0.0 {
        prune_basis_complex_local(basis, coeffs, drop_tol, protected);
    }
    cap_basis_complex(basis, coeffs, max_basis, protected);
    Ok(())
}

/// Complex analogue of `crate::add_leakage_capped`: add the largest leakage
/// reps to the basis, up to the available room `room = max_basis −
/// basis.len()`, so the in-step orbit-rep basis never exceeds `max_basis`.
/// New reps get coefficient 0; the surrounding expm fills them. No
/// magnitude filter — the top-`room` by `|leakage|` are added.
fn add_leakage_capped_complex(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    mut leak: Vec<(Word, Complex<f64>)>,
    max_basis: usize,
) {
    let room = max_basis.saturating_sub(basis.len());
    if leak.len() > room {
        if room > 0 {
            leak.select_nth_unstable_by(room - 1, |a, b| {
                b.1.norm()
                    .partial_cmp(&a.1.norm())
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        leak.truncate(room);
    }
    for (w, _) in leak {
        basis.push(w);
        coeffs.push(Complex::new(0.0, 0.0));
    }
}

/// Complex analogue of `crate::cap_basis`: keep only the `max_basis`
/// largest-`|c|` reps (protected reps always kept), dropping the rest.
/// A `max_basis` large enough to cover the whole basis is a no-op.
fn cap_basis_complex(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    max_basis: usize,
    protected: &[Word],
) {
    if basis.len() <= max_basis {
        return;
    }
    let protected_set: fxhash::FxHashSet<&Word> = protected.iter().collect();
    let n_prot = basis.iter().filter(|w| protected_set.contains(w)).count();
    let slots = max_basis.saturating_sub(n_prot);
    let mut mags: Vec<f64> = basis
        .iter()
        .zip(coeffs.iter())
        .filter(|(w, _)| !protected_set.contains(w))
        .map(|(_, c)| c.norm())
        .collect();
    let cutoff = if slots == 0 {
        f64::INFINITY
    } else if slots >= mags.len() {
        return;
    } else {
        let k = slots - 1;
        mags.select_nth_unstable_by(k, |a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        mags[k]
    };
    let mut write = 0;
    for read in 0..basis.len() {
        if protected_set.contains(&basis[read]) || coeffs[read].norm() >= cutoff {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}

fn prune_basis_complex_local(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    drop_tol: f64,
    protected: &[Word],
) {
    if drop_tol <= 0.0 {
        return;
    }
    let protected_set: fxhash::FxHashSet<&Word> = protected.iter().collect();
    let mut write = 0;
    for read in 0..basis.len() {
        if coeffs[read].norm() >= drop_tol || protected_set.contains(&basis[read]) {
            if write != read {
                basis.swap(write, read);
                coeffs.swap(write, read);
            }
            write += 1;
        }
    }
    basis.truncate(write);
    coeffs.truncate(write);
}
