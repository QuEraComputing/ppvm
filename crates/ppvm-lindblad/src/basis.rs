// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Basis-level `L*` operators: in-basis generator and off-basis leakage.

use crate::Error;
use crate::sector::Sector;
use crate::spec::LindbladSpec;
use crate::truncate::{cap_map_to_room, order_by_desc_mag};
use crate::word::{Word, word_hash};
use fxhash::{FxBuildHasher, FxHashMap, FxHashSet};
use num::Complex;
use rayon::prelude::*;

/// Chunk size for the leakage accumulation loops: candidates are folded
/// into the live map (and the room-cap applied) once per chunk.
const CHUNK_SIZE: usize = 4096;

/// Reusable scratch for adaptive real-coefficient evolution.
///
/// Keep one workspace per trajectory. It retains candidate payload and generation
/// buffers; hash-table growth starts from an empty allocation to preserve ordering.
#[derive(Default)]
pub struct PcStepWorkspace<const C: usize = { crate::W_CHUNKS }> {
    pub(crate) candidates: Option<crate::indexed_candidates::IndexedCandidates<C>>,
    outputs: Vec<Vec<(Word<C>, f64)>>,
    magnitude_bits: Vec<u64>,
}

pub(crate) enum LeakageCandidates<const C: usize> {
    Pairs(Vec<(Word<C>, f64)>),
    Indexed(crate::indexed_candidates::IndexedCandidates<C>, f64),
}

impl<const C: usize> LeakageCandidates<C> {
    fn into_vec(self) -> Vec<(Word<C>, f64)> {
        match self {
            Self::Pairs(pairs) => pairs,
            Self::Indexed(indexed, threshold) => indexed.into_vec(threshold),
        }
    }
    pub(crate) fn admit(
        self,
        basis: &mut Vec<Word<C>>,
        coefficients: &mut Vec<f64>,
        maximum: usize,
    ) -> Option<crate::indexed_candidates::IndexedCandidates<C>> {
        match self {
            Self::Pairs(pairs) => {
                crate::truncate::add_leakage_capped(basis, coefficients, pairs, maximum);
                None
            }
            Self::Indexed(mut indexed, threshold) => {
                indexed.admit(basis, coefficients, maximum, threshold);
                indexed.reset();
                Some(indexed)
            }
        }
    }
}

#[cfg(test)]
mod buffered_tests {
    use super::*;
    use crate::{LeakageBackend, word_from_codes};

    #[test]
    fn buffered_leakage_preserves_multiple_chunk_caps_and_source_order() {
        let n = 7;
        let mut h = Vec::new();
        for site in 0..n - 1 {
            for symbol in ['X', 'Y'] {
                let mut term = vec!['I'; n];
                term[site] = symbol;
                term[site + 1] = symbol;
                h.push((term.into_iter().collect(), 0.3));
            }
        }
        let reference = LindbladSpec::<2>::new(n, &h, &[]).unwrap();
        let buffered = LindbladSpec::<2>::new(n, &h, &[])
            .unwrap()
            .with_leakage_backend(LeakageBackend::Buffered);
        let indexed = LindbladSpec::<2>::new(n, &h, &[])
            .unwrap()
            .with_leakage_backend(LeakageBackend::Indexed);
        let basis: Vec<_> = (0..14000)
            .filter(|i| i % 3 != 0)
            .map(|i| {
                word_from_codes::<2>(
                    &(0..n)
                        .map(|q| ((i >> (2 * q)) & 3) as u8)
                        .collect::<Vec<_>>(),
                )
                .unwrap()
            })
            .collect();
        let coefficients: Vec<_> = (0..basis.len()).map(|i| (i as f64).sin()).collect();
        for room in [0, 100, 2000] {
            let expected = reference
                .leakage_with_prune(&basis, &coefficients, &[], basis.len() + room, 0.0)
                .unwrap();
            for threads in [1, 2, 6] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let actual = pool
                    .install(|| {
                        buffered.leakage_with_prune(
                            &basis,
                            &coefficients,
                            &[],
                            basis.len() + room,
                            0.0,
                        )
                    })
                    .unwrap();
                assert_eq!(actual, expected, "room={room}, threads={threads}");
                for implementation in [&indexed] {
                    let actual = pool
                        .install(|| {
                            implementation.leakage_with_prune(
                                &basis,
                                &coefficients,
                                &[],
                                basis.len() + room,
                                0.0,
                            )
                        })
                        .unwrap();
                    assert_eq!(actual, expected, "room={room}, threads={threads}");
                }
            }
        }
    }
}

/// Build a `word → row` map for a basis assumed to contain unique Pauli
/// words; debug-asserts the uniqueness invariant.
pub fn build_basis_index<const C: usize>(basis: &[Word<C>]) -> FxHashMap<Word<C>, u32> {
    let mut index: FxHashMap<Word<C>, u32> = FxHashMap::default();
    for (i, w) in basis.iter().enumerate() {
        let prev = index.insert(*w, i as u32);
        debug_assert!(
            prev.is_none(),
            "basis contains duplicate Pauli word at positions {} and {}",
            prev.unwrap(),
            i,
        );
    }
    index
}

impl<const C: usize> LindbladSpec<C> {
    /// Off-basis component of `L*( Σ_j coeffs[j] · basis[j] )`. Output
    /// strings that lie in `basis` or in `protected` are dropped.
    pub fn leakage(
        &self,
        basis: &[Word<C>],
        coeffs: &[f64],
        protected: &[Word<C>],
    ) -> Result<Vec<(Word<C>, f64)>, Error> {
        self.leakage_with_prune(basis, coeffs, protected, usize::MAX, 0.0)
    }

    /// Like [`Self::leakage`], but caps the live off-basis leakage map to
    /// the *available room* `room = max_basis − basis.len()` — only the
    /// strings we could actually add to the basis are worth keeping. The
    /// cap is applied during accumulation (after each chunk), keeping the
    /// `room` largest-magnitude entries.
    ///
    /// Basis indices are processed in descending-`|c|` order so the
    /// running cap keeps the entries that are most likely to be the true
    /// largest contributors. When `max_basis` is large enough that
    /// `room ≥ all candidates`, nothing is dropped — the near-exact case.
    pub fn leakage_with_prune(
        &self,
        basis: &[Word<C>],
        coeffs: &[f64],
        protected: &[Word<C>],
        max_basis: usize,
        tau_add: f64,
    ) -> Result<Vec<(Word<C>, f64)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        if self.leakage_backend != crate::LeakageBackend::Reference {
            return self
                .leakage_buffered(basis, coeffs, protected, max_basis, tau_add, None, None)
                .map(LeakageCandidates::into_vec);
        }
        // Hash-only membership tables: storing 8-byte `u64` keys instead
        // of 48-byte Words shrinks the in-basis structure ~6×, keeping it
        // in L3 (and often L2) at basis sizes where the full-Word version
        // would spill to DRAM.
        let in_basis: FxHashMap<u64, ()> = basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        let order = order_by_desc_mag(coeffs);
        let room = max_basis.saturating_sub(basis.len());
        let n_qubits = self.n_qubits();
        let mut merged: FxHashMap<Word<C>, f64> = FxHashMap::default();
        for chunk_indices in order.chunks(CHUNK_SIZE) {
            let local: Vec<Vec<(Word<C>, f64)>> = chunk_indices
                .par_iter()
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word<C>, Complex<f64>>::with_capacity_and_hasher(
                                128,
                                FxBuildHasher::default(),
                            ),
                        )
                    },
                    |(s1, s2, lm), &i| {
                        let p = &basis[i];
                        let c = coeffs[i];
                        let terms = self.compute_action_terms(p, s1, s2, lm);
                        let mut out = Vec::with_capacity(terms.len());
                        for (w, v) in terms.iter() {
                            let h = word_hash(w);
                            if !in_basis.contains_key(&h) && !protected_set.contains_key(&h) {
                                out.push((*w, c * *v));
                            }
                        }
                        out
                    },
                )
                .collect();
            for v in local {
                for (k, val) in v {
                    *merged.entry(k).or_insert(0.0) += val;
                }
            }
            cap_map_to_room(&mut merged, room);
        }
        // Rate-based admission: keep only candidates whose leakage rate
        // exceeds `tau_add`. `tau_add = 0` admits everything except exact
        // zeros.
        Ok(merged
            .into_iter()
            .filter(|(_, c)| c.abs() > tau_add)
            .collect())
    }

    pub(crate) fn leakage_for_step(
        &self,
        basis: &mut Vec<Word<C>>,
        coefficients: &[f64],
        protected: &[Word<C>],
        maximum: usize,
        threshold: f64,
        workspace: Option<&mut PcStepWorkspace<C>>,
    ) -> Result<LeakageCandidates<C>, Error> {
        if !matches!(
            self.leakage_backend,
            crate::LeakageBackend::Indexed | crate::LeakageBackend::Tail
        ) {
            return self
                .leakage_with_prune(basis, coefficients, protected, maximum, threshold)
                .map(LeakageCandidates::Pairs);
        }
        if basis.len() != coefficients.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coefficients.len(),
            });
        }
        let owned =
            (self.leakage_backend == crate::LeakageBackend::Tail).then(|| std::mem::take(basis));
        self.leakage_buffered(
            basis,
            coefficients,
            protected,
            maximum,
            threshold,
            workspace,
            owned,
        )
    }

    fn leakage_buffered(
        &self,
        basis: &[Word<C>],
        coeffs: &[f64],
        protected: &[Word<C>],
        max_basis: usize,
        tau_add: f64,
        mut workspace: Option<&mut PcStepWorkspace<C>>,
        owned_basis: Option<Vec<Word<C>>>,
    ) -> Result<LeakageCandidates<C>, Error> {
        let tail = owned_basis.is_some();
        let source = owned_basis.as_deref().unwrap_or(basis);
        let in_basis: FxHashMap<u64, ()> = source.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();
        let order = order_by_desc_mag(coeffs);
        let room = max_basis.saturating_sub(source.len());
        let mut merged = FxHashMap::default();
        let mut indexed = workspace
            .as_mut()
            .and_then(|w| w.candidates.take())
            .unwrap_or_else(crate::indexed_candidates::IndexedCandidates::new);
        if let Some(owned) = owned_basis {
            indexed.attach_basis(owned);
        }
        let mut magnitudes = Vec::new();
        let mut magnitude_bits = workspace
            .as_mut()
            .map(|w| std::mem::take(&mut w.magnitude_bits))
            .unwrap_or_default();
        let mut outputs = workspace
            .as_mut()
            .map(|w| std::mem::take(&mut w.outputs))
            .unwrap_or_default();
        outputs.resize_with(CHUNK_SIZE.div_ceil(64), Vec::new);
        for chunk in order.chunks(CHUNK_SIZE) {
            let basis = if tail { indexed.source_basis() } else { basis };
            outputs
                .par_iter_mut()
                .zip(chunk.par_chunks(64))
                .for_each_init(
                    crate::direct_expm::ActionScratch::new,
                    |scratch, (output, indices)| {
                        output.clear();
                        for &index in indices {
                            scratch.fill(self, &basis[index], true);
                            for &(word, value) in &scratch.terms {
                                let hash = word_hash(&word);
                                if !in_basis.contains_key(&hash)
                                    && !protected_set.contains_key(&hash)
                                {
                                    output.push((word, coeffs[index] * value));
                                }
                            }
                        }
                    },
                );
            // Serial source order preserves coherent sums and the reference map's tie ordering.
            for output in outputs.iter().take(chunk.len().div_ceil(64)) {
                for &(word, value) in output {
                    if matches!(
                        self.leakage_backend,
                        crate::LeakageBackend::Indexed | crate::LeakageBackend::Tail
                    ) {
                        indexed.add(word, value);
                    } else {
                        *merged.entry(word).or_insert(0.0) += value;
                    }
                }
            }
            if matches!(
                self.leakage_backend,
                crate::LeakageBackend::Indexed | crate::LeakageBackend::Tail
            ) {
                indexed.cap(room, &mut magnitude_bits);
            } else {
                crate::truncate::cap_real_map_with_scratch(&mut merged, room, &mut magnitudes);
            }
        }
        // Generation buffers are dead before the admission vector is materialized.
        if let Some(workspace) = workspace {
            for output in &mut outputs {
                output.clear();
            }
            workspace.outputs = outputs;
            magnitude_bits.clear();
            workspace.magnitude_bits = magnitude_bits;
        }
        drop(magnitudes);
        drop(order);
        drop(in_basis);
        drop(protected_set);
        if matches!(
            self.leakage_backend,
            crate::LeakageBackend::Indexed | crate::LeakageBackend::Tail
        ) {
            return Ok(LeakageCandidates::Indexed(indexed, tau_add));
        }
        Ok(LeakageCandidates::Pairs(
            merged
                .into_iter()
                .filter(|(_, value)| value.abs() > tau_add)
                .collect(),
        ))
    }

    /// Sparse generator matrix in COO form: returns `(row, col, val)`
    /// triplets. Row = output Pauli's position in `basis`; col = input
    /// Pauli's position. Output Paulis not in `basis` are silently dropped.
    ///
    /// Precondition: `basis` must not contain duplicate Pauli words
    /// (asserted in debug builds).
    pub fn generator(&self, basis: &[Word<C>]) -> Vec<(usize, usize, f64)> {
        let index = build_basis_index(basis);
        let n_qubits = self.n_qubits();

        // `compute_action_terms` returns a deduplicated `Vec<(Word, f64)>`,
        // so it can be scattered directly into COO triplets.
        let local: Vec<Vec<(usize, usize, f64)>> = basis
            .par_iter()
            .enumerate()
            .map_init(
                || {
                    (
                        Vec::<u32>::with_capacity(n_qubits),
                        Vec::<u32>::with_capacity(128),
                        FxHashMap::<Word<C>, Complex<f64>>::with_capacity_and_hasher(
                            128,
                            FxBuildHasher::default(),
                        ),
                    )
                },
                |(s1, s2, lm), (col, p)| {
                    let terms = self.compute_action_terms(p, s1, s2, lm);
                    let mut out = Vec::with_capacity(terms.len());
                    for (w, v) in terms.iter() {
                        if let Some(&row) = index.get(w) {
                            out.push((row as usize, col, *v));
                        }
                    }
                    out
                },
            )
            .collect();

        // Pre-allocate the flat output to avoid sequential push reallocation.
        let total: usize = local.iter().map(|v| v.len()).sum();
        let mut flat = Vec::with_capacity(total);
        for v in local {
            flat.extend(v);
        }
        flat
    }

    /// Complex-coefficient variant of [`Self::leakage`]: off-basis
    /// component of `L*( Σ_j coeffs[j] · basis[j] )` with complex `coeffs`.
    pub fn leakage_complex(
        &self,
        basis: &[Word<C>],
        coeffs: &[Complex<f64>],
        protected: &[Word<C>],
    ) -> Result<Vec<(Word<C>, Complex<f64>)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        let in_basis: FxHashMap<u64, ()> = basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        let n_qubits = self.n_qubits();
        let mut merged: FxHashMap<Word<C>, Complex<f64>> = FxHashMap::default();
        for chunk_start in (0..basis.len()).step_by(CHUNK_SIZE) {
            let chunk_end = (chunk_start + CHUNK_SIZE).min(basis.len());
            let chunk_basis = &basis[chunk_start..chunk_end];
            let chunk_coeffs = &coeffs[chunk_start..chunk_end];
            let local: Vec<Vec<(Word<C>, Complex<f64>)>> = chunk_basis
                .par_iter()
                .zip(chunk_coeffs.par_iter())
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word<C>, Complex<f64>>::with_capacity_and_hasher(
                                128,
                                FxBuildHasher::default(),
                            ),
                        )
                    },
                    |(s1, s2, lm), (p, &c)| {
                        let terms = self.compute_action_terms(p, s1, s2, lm);
                        let mut out = Vec::with_capacity(terms.len());
                        for (w, v) in terms.iter() {
                            let h = word_hash(w);
                            if !in_basis.contains_key(&h) && !protected_set.contains_key(&h) {
                                out.push((*w, c * *v));
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
        }
        Ok(merged.into_iter().filter(|(_, c)| c.norm() > 0.0).collect())
    }

    /// Phase-aware leakage: out-of-basis component of `L*(O_k)` where
    /// `O_k` is the operator represented by `basis` (orbit reps) and
    /// `coeffs` (complex coefficients in momentum `sector`).
    ///
    /// For each input rep `r` with coefficient `c_r`, and each output `q`
    /// of `L*(r) = Σ_q v_q · q`:
    /// 1. Canonicalize `q` → `(r_q, χ_k, |orbit_q|)` via
    ///    [`Sector::canonicalize_phase`].
    /// 2. If `r_q` NOT in `basis` and NOT in `protected`:
    ///    `merged[r_q] += χ_k · v_q · c_r · |orbit_r| / |orbit_q|`.
    ///
    /// The `|orbit_r| / |orbit_q|` factor is the convention conversion
    /// documented on [`crate::mf_expm`]'s `build_orbit_rep_cols`, so the
    /// admitted rates are comparable to the averaged-convention
    /// coefficients the caller holds.
    ///
    /// Returns `(r_q, sum)` pairs for all candidates with nonzero sum.
    ///
    /// This is the orbit-rep counterpart of [`Self::leakage_with_prune`],
    /// and caps the live candidate map the same way: to the *available
    /// room* `room = max_basis − basis.len()` (the reps we could actually
    /// add), applied during accumulation. A large `max_basis`
    /// (room ≥ all candidates) disables the cap — the near-exact case.
    pub fn leakage_orbit_rep(
        &self,
        basis: &[Word<C>],
        coeffs: &[Complex<f64>],
        protected: &[Word<C>],
        sector: &Sector<'_>,
        max_basis: usize,
    ) -> Result<Vec<(Word<C>, Complex<f64>)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        // Membership is tested on the canonical rep `r_q`, so unlike the
        // real path these are full-Word sets, not `word_hash` tables.
        let in_basis: FxHashSet<&Word<C>> = basis.iter().collect();
        let protected_set: FxHashSet<&Word<C>> = protected.iter().collect();

        let order = order_by_desc_mag(coeffs);
        let room = max_basis.saturating_sub(basis.len());
        let n_qubits = self.n_qubits();
        let mut merged: FxHashMap<Word<C>, Complex<f64>> = FxHashMap::default();
        for chunk_indices in order.chunks(CHUNK_SIZE) {
            let local: Vec<Vec<(Word<C>, Complex<f64>)>> = chunk_indices
                .par_iter()
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word<C>, Complex<f64>>::with_capacity_and_hasher(
                                128,
                                FxBuildHasher::default(),
                            ),
                        )
                    },
                    |(s1, s2, lm), &i| {
                        let r = &basis[i];
                        let c_r = coeffs[i];
                        // A rep that cannot carry the sector contributes
                        // nothing (its coefficient is identically zero).
                        let Some(orbit_in) = sector.orbit_size(r) else {
                            return Vec::new();
                        };
                        let terms = self.compute_action_terms(r, s1, s2, lm);
                        let mut out = Vec::with_capacity(terms.len());
                        for (q, v) in terms.iter() {
                            let Some((r_q, phase, orbit_out)) = sector.canonicalize_phase(q) else {
                                continue;
                            };
                            if !in_basis.contains(&r_q) && !protected_set.contains(&r_q) {
                                let rate = phase * *v * c_r * (orbit_in as f64 / orbit_out as f64);
                                out.push((r_q, rate));
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
            cap_map_to_room(&mut merged, room);
        }
        Ok(merged.into_iter().filter(|(_, c)| c.norm() > 0.0).collect())
    }
}
