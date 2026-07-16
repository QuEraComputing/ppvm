// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Direct Heisenberg-picture Lindbladian evolution on an adaptive
//! Pauli-string basis.
//!
//! For a Hermitian Pauli Hamiltonian `H = Σ c_i P_i` and jump operators
//! `L_k = Σ_a λ_{k,a} P_{k,a}` (each a Hermitian-Pauli linear combination
//! with possibly complex coefficients) with rates `γ_k ≥ 0`, the adjoint
//! Lindbladian acts on a single Pauli string `p` as
//!
//! ```text
//! L*(p) = i [H, p] + Σ_k γ_k ( L_k† p L_k − 1/2 {L_k† L_k, p} ).
//! ```
//!
//! Two jump shapes are supported with separate code paths:
//!
//! - **Hermitian Pauli** (`L = P`, `λ ∈ ℝ`): the dissipator collapses to a
//!   diagonal `-2γ` on Pauli strings that anti-commute with `P`. Same fast
//!   path used by every dephasing-style model.
//!
//! - **General** (complex `λ_a`, e.g. `σ± = (X ± iY)/2`): the dissipator
//!   becomes a double sum `Σ_{a,b} λ_a* λ_b P_a p P_b` plus a Pauli-
//!   linear-combination anti-commutator with `L†L`, which is precomputed
//!   once at construction. Intermediate coefficients are complex; the
//!   result is real because `L*` preserves Hermiticity, so we cast back
//!   to `f64` at the boundary (with a debug-only check that `|Im|` is at
//!   FP noise).
//!
//! Pauli strings are stored as [`ppvm_pauli_word::word::PauliWord`] backed by
//! `[u64; 2]` (≤128 qubits) with cached hashes for fast HashMap lookup. The
//! hot-path commutator/product loops bypass the higher-level word API and
//! operate directly on raw `u64` chunks for speed.

use fxhash::{FxBuildHasher, FxHashMap, FxHashSet};
use num::Complex;
use ppvm_traits::PauliWordTrait;
use ppvm_pauli_word::word::PauliWord;
use rayon::prelude::*;
use std::time::Instant;

pub mod config;
pub mod error;
pub(crate) mod expm;

/// Matrix-free / quspin-expm-backed `exp(dt·L*)·b` engine. See module docs.
pub(crate) mod mf_expm;

/// Per-step orbit-rep evolution under translation symmetry, with a
/// phase-aware complex action. See module docs.
pub mod orbit_rep;

/// Words pack up to 128 qubits.
const W_U64: usize = 2;

/// Maximum number of qubits supported by [`Word`] (= `8 · W_U64 · sizeof(u64)`).
pub const MAX_QUBITS: usize = 64 * W_U64;

/// The Pauli-word storage type used throughout this crate.
///
/// `[u64; 2]` covers up to 128 qubits; the `FxBuildHasher` matches the
/// hash used by the `FxHashMap` keys we wrap with; `REHASH=true` means
/// `set()` keeps the cached hash in sync.
pub type Word = PauliWord<[u64; W_U64], FxBuildHasher, true>;

pub use config::PcStepConfig;
pub use error::Error;

/// Per-phase timing breakdown (microseconds) returned by
/// [`LindbladSpec::pc_step_timed`].
#[derive(Default, Clone, Copy, Debug)]
pub struct PcStepTimings {
    pub leakage1_us: u64,
    pub expand1_us: u64,
    pub expm1_us: u64,
    pub leakage2_us: u64,
    pub expand2_us: u64,
    pub expm2_us: u64,
}

impl PcStepTimings {
    pub fn total_us(&self) -> u64 {
        self.leakage1_us
            + self.expand1_us
            + self.expm1_us
            + self.leakage2_us
            + self.expand2_us
            + self.expm2_us
    }
}

// ────────────────── codec helpers ──────────────────

/// Build a [`Word`] from a length-`n_qubits` slice of Pauli codes
/// (`0=I, 1=X, 2=Z, 3=Y`). Sets all bits and rehashes once.
pub fn word_from_codes(codes: &[u8]) -> Result<Word, Error> {
    let n_qubits = codes.len();
    if n_qubits > MAX_QUBITS {
        return Err(Error::TooManyQubits { got: n_qubits });
    }
    let mut w = Word::new(n_qubits);
    for (q, &b) in codes.iter().enumerate() {
        if b > 3 {
            return Err(Error::InvalidPauliCode { code: b });
        }
        if b & 1 != 0 {
            w.xbits.set(q, true);
        }
        if b & 2 != 0 {
            w.zbits.set(q, true);
        }
    }
    w.rehash();
    Ok(w)
}

/// Inverse of [`word_from_codes`]: write `n_qubits` Pauli codes into `out`.
pub fn codes_from_word(w: &Word, out: &mut [u8]) {
    debug_assert_eq!(out.len(), w.n_qubits());
    for (q, slot) in out.iter_mut().enumerate() {
        let xb = w.xbits[q] as u8;
        let zb = w.zbits[q] as u8;
        *slot = xb | (zb << 1);
    }
}

/// Parse a `"IXYZ..."` string into a [`Word`] together with the list of
/// qubits where the Pauli is non-identity (the term's support).
pub fn parse_pauli_string(s: &str, n_qubits: usize) -> Result<(Word, Vec<u32>), Error> {
    if n_qubits > MAX_QUBITS {
        return Err(Error::TooManyQubits { got: n_qubits });
    }
    let chars: Vec<char> = s.chars().filter(|c| *c != '_').collect();
    if chars.len() != n_qubits {
        return Err(Error::WrongLength {
            expected: n_qubits,
            got: chars.len(),
        });
    }
    let mut w = Word::new(n_qubits);
    let mut support = Vec::new();
    for (q, c) in chars.into_iter().enumerate() {
        match c {
            'I' => {}
            'X' => {
                w.xbits.set(q, true);
                support.push(q as u32);
            }
            'Z' => {
                w.zbits.set(q, true);
                support.push(q as u32);
            }
            'Y' => {
                w.xbits.set(q, true);
                w.zbits.set(q, true);
                support.push(q as u32);
            }
            other => return Err(Error::InvalidPauliChar { c: other }),
        }
    }
    w.rehash();
    Ok((w, support))
}

/// Compute the support (non-identity qubits) of `w`.
fn word_support(w: &Word, out: &mut Vec<u32>) {
    out.clear();
    for q in 0..w.n_qubits() {
        if w.xbits[q] || w.zbits[q] {
            out.push(q as u32);
        }
    }
}

// ────────────────── Pauli algebra ──────────────────
//
// Phase encoding: `Pauli product P·Q = ω · R` where `ω = i^phase` and
// `phase ∈ {0,1,2,3}` ↔ `ω ∈ {1, i, -1, -i}`. The per-byte XOR/AND
// formulas are the same ones used by `PhasedPauliWord::mul_assign` in
// `ppvm-runtime`.

#[inline(always)]
fn phase_factor(phase: u8) -> Complex<f64> {
    match phase & 3 {
        0 => Complex::new(1.0, 0.0),
        1 => Complex::new(0.0, 1.0),
        2 => Complex::new(-1.0, 0.0),
        _ => Complex::new(0.0, -1.0),
    }
}

/// `true` if Pauli words `a` and `b` anti-commute.
///
/// Two Pauli strings anti-commute iff
/// `popcount(a.x & b.z) + popcount(a.z & b.x)` is odd.
#[inline(always)]
fn anti_commutes(a: &Word, b: &Word) -> bool {
    let mut bits: u32 = 0;
    for i in 0..W_U64 {
        bits += (a.xbits.data[i] & b.zbits.data[i]).count_ones();
        bits += (a.zbits.data[i] & b.xbits.data[i]).count_ones();
    }
    bits & 1 == 1
}

/// Commutator product `h · p`: returns `(out, eps)` where `out = h ⊕ p` and
///
/// - `eps =  0` if `h` and `p` commute (caller should skip — `[h,p] = 0`),
/// - `eps = -2.0` if `h·p` has phase `+i` (so `i·[h,p] = -2·out`),
/// - `eps = +2.0` if `h·p` has phase `-i` (so `i·[h,p] = +2·out`).
#[inline(always)]
fn comm_product(h: &Word, p: &Word) -> (Word, f64) {
    let (out, phase) = pauli_mul(h, p);
    let eps = match phase {
        1 => -2.0,
        3 => 2.0,
        _ => 0.0,
    };
    (out, eps)
}

/// Full Pauli product `p · q`: returns `(out, phase)` where the product
/// is `ω · out` with `ω = i^phase`.
#[inline(always)]
fn pauli_mul(p: &Word, q: &Word) -> (Word, u8) {
    let mut out = Word::new(p.n_qubits());
    let mut sign_count: u32 = 0;
    let mut imag_count: u32 = 0;
    for i in 0..W_U64 {
        let a = p.xbits.data[i];
        let b = p.zbits.data[i];
        let c = q.xbits.data[i];
        let d = q.zbits.data[i];
        let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
        let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
        sign_count += sign.count_ones();
        imag_count += imag.count_ones();
        out.xbits.data[i] = a ^ c;
        out.zbits.data[i] = b ^ d;
    }
    out.rehash();
    (out, ((2 * sign_count + imag_count) & 3) as u8)
}

// ────────────────── Spec types ──────────────────

/// Parsed Hamiltonian term.
#[derive(Clone)]
struct HTerm {
    word: Word,
    coeff: f64,
}

/// One Pauli term in a complex linear combination (a single summand of
/// `L = Σ_a λ_a P_a` or of the precomputed `L†L`).
#[derive(Clone)]
struct PauliTerm {
    word: Word,
    coeff: Complex<f64>,
}

/// One jump operator `L_k` with rate `γ_k`. The `HermitianPauli` variant
/// is a fast path; `General` handles arbitrary complex Pauli sums.
#[derive(Clone)]
enum JumpKind {
    HermitianPauli {
        word: Word,
        rate: f64,
    },
    General {
        terms: Vec<PauliTerm>,         // L = Σ_a λ_a P_a
        dagger_dagger: Vec<PauliTerm>, // L†L = Σ_c μ_c P_c  (μ_c ∈ ℝ)
        rate: f64,
    },
}

/// Expand `L†L = (Σ_a λ_a P_a)† (Σ_b λ_b P_b) = Σ_{a,b} λ_a* λ_b P_a P_b`
/// as a Pauli linear combination, dropping FP-noise zeros. Coefficients are
/// real because `L†L` is Hermitian; we keep them complex for arithmetic
/// uniformity.
fn precompute_ldagger_l(terms: &[PauliTerm]) -> Vec<PauliTerm> {
    let zero = Complex::new(0.0, 0.0);
    let mut acc: FxHashMap<Word, Complex<f64>> = FxHashMap::default();
    for a in terms {
        for b in terms {
            let (word, phase) = pauli_mul(&a.word, &b.word);
            let coeff = a.coeff.conj() * b.coeff * phase_factor(phase);
            *acc.entry(word).or_insert(zero) += coeff;
        }
    }
    acc.into_iter()
        .filter(|(_, c)| c.norm() > 1e-14)
        .map(|(word, coeff)| PauliTerm { word, coeff })
        .collect()
}

/// Precompiled Lindbladian. Constructed once from string-form Hamiltonian
/// terms + jump operators; reused across many calls to [`Self::action`],
/// [`Self::leakage`], [`Self::generator`]. `L*(p)` is recomputed on every
/// call rather than cached: for sparse-local Hamiltonians a per-word cache
/// costs more than the recompute (hash lookup ≳ recompute) and its several
/// KB per cached word dominate memory at large basis sizes.
pub struct LindbladSpec {
    n_qubits: usize,
    h_terms: Vec<HTerm>,
    j_kinds: Vec<JumpKind>,
    /// `h_support[q]` = indices of Hamiltonian terms acting on qubit `q`.
    h_support: Vec<Vec<u32>>,
    /// `j_support[q]` = indices of jumps whose support contains qubit `q`.
    j_support: Vec<Vec<u32>>,
}

/// User-facing description of one jump operator: a complex Pauli linear
/// combination together with its rate.
#[derive(Clone, Debug)]
pub struct JumpInput {
    /// `(pauli_string, λ)` pairs forming `L_k = Σ_a λ_a P_a`.
    pub lincomb: Vec<(String, Complex<f64>)>,
    /// Non-negative GKSL rate `γ_k`.
    pub rate: f64,
}

impl LindbladSpec {
    /// Construct a Lindbladian spec from Hamiltonian terms and jump operators.
    ///
    /// `h_terms` are `(pauli_string, coefficient)` pairs forming the Hermitian
    /// Hamiltonian. Each jump operator is a complex Pauli linear combination;
    /// a length-1 jump with imaginary part `0` is routed to the Hermitian-Pauli
    /// fast path (with rate scaled by the squared real coefficient).
    pub fn new(
        n_qubits: usize,
        h_terms: &[(String, f64)],
        jumps: &[JumpInput],
    ) -> Result<Self, Error> {
        if n_qubits > MAX_QUBITS {
            return Err(Error::TooManyQubits { got: n_qubits });
        }

        let mut h_parsed: Vec<HTerm> = Vec::with_capacity(h_terms.len());
        let mut h_support_idx: Vec<Vec<u32>> = vec![Vec::new(); n_qubits];
        for (i, (s, c)) in h_terms.iter().enumerate() {
            let (word, support) = parse_pauli_string(s, n_qubits)?;
            for q in support {
                h_support_idx[q as usize].push(i as u32);
            }
            h_parsed.push(HTerm { word, coeff: *c });
        }

        let mut j_kinds: Vec<JumpKind> = Vec::with_capacity(jumps.len());
        let mut j_support_idx: Vec<Vec<u32>> = vec![Vec::new(); n_qubits];
        for (k, jump) in jumps.iter().enumerate() {
            if jump.rate < 0.0 {
                return Err(Error::NegativeRate {
                    index: k,
                    rate: jump.rate,
                });
            }
            if jump.lincomb.is_empty() {
                return Err(Error::EmptyLincomb { index: k });
            }

            // Fast path: single-term, purely real → Hermitian Pauli.
            if jump.lincomb.len() == 1 && jump.lincomb[0].1.im == 0.0 {
                let (s, c) = &jump.lincomb[0];
                let (word, support) = parse_pauli_string(s, n_qubits)?;
                for q in support {
                    j_support_idx[q as usize].push(k as u32);
                }
                j_kinds.push(JumpKind::HermitianPauli {
                    word,
                    rate: jump.rate * c.re * c.re,
                });
                continue;
            }

            // General path: parse all terms, precompute L†L, record union support.
            let mut terms: Vec<PauliTerm> = Vec::with_capacity(jump.lincomb.len());
            let mut union_support: std::collections::BTreeSet<u32> =
                std::collections::BTreeSet::new();
            for (s, c) in &jump.lincomb {
                let (word, support) = parse_pauli_string(s, n_qubits)?;
                for q in &support {
                    union_support.insert(*q);
                }
                terms.push(PauliTerm { word, coeff: *c });
            }
            for q in union_support {
                j_support_idx[q as usize].push(k as u32);
            }
            let dagger_dagger = precompute_ldagger_l(&terms);
            j_kinds.push(JumpKind::General {
                terms,
                dagger_dagger,
                rate: jump.rate,
            });
        }

        Ok(Self {
            n_qubits,
            h_terms: h_parsed,
            j_kinds,
            h_support: h_support_idx,
            j_support: j_support_idx,
        })
    }

    pub fn n_qubits(&self) -> usize {
        self.n_qubits
    }

    pub fn num_h_terms(&self) -> usize {
        self.h_terms.len()
    }

    pub fn num_jump_terms(&self) -> usize {
        self.j_kinds.len()
    }

    /// Apply `L*` to a single Pauli string `p`. Returns the output Pauli
    /// strings and their real coefficients (zero entries omitted).
    pub fn action(&self, p: &Word) -> Vec<(Word, f64)> {
        let mut out: FxHashMap<Word, f64> = FxHashMap::default();
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        self.accumulate_action(p, 1.0, &mut out, &mut s1, &mut s2);
        out.into_iter().filter(|(_, c)| *c != 0.0).collect()
    }

    /// Off-basis component of `L*( Σ_j coeffs[j] · basis[j] )`. Output
    /// strings that lie in `basis` or in `protected` are dropped.
    pub fn leakage(
        &self,
        basis: &[Word],
        coeffs: &[f64],
        protected: &[Word],
    ) -> Result<Vec<(Word, f64)>, Error> {
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
        basis: &[Word],
        coeffs: &[f64],
        protected: &[Word],
        max_basis: usize,
        tau_add: f64,
    ) -> Result<Vec<(Word, f64)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        // Hash-only membership tables: storing 8-byte `u64` keys instead
        // of 48-byte Words shrinks the in-basis structure ~6×, keeping it
        // in L3 (and often L2) at basis sizes where the full-Word version
        // would spill to DRAM.
        let in_basis: FxHashMap<u64, ()> = basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        // Descending sort by |c|: process largest-magnitude contributors
        // first so the running room-cap keeps the right entries.
        let mut order: Vec<usize> = (0..basis.len()).collect();
        order.sort_by(|&a, &b| {
            coeffs[b]
                .abs()
                .partial_cmp(&coeffs[a].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        const CHUNK_SIZE: usize = 4096;
        let room = max_basis.saturating_sub(basis.len());
        let mut merged: FxHashMap<Word, f64> = FxHashMap::default();
        for chunk_indices in order.chunks(CHUNK_SIZE) {
            let local: Vec<Vec<(Word, f64)>> = chunk_indices
                .par_iter()
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(self.n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
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

            // Room-cap: keep only the `room` largest-magnitude entries.
            if merged.len() > room {
                if room == 0 {
                    merged.clear();
                } else {
                    let mut mags: Vec<f64> = merged.values().map(|v| v.abs()).collect();
                    let k = room.min(mags.len() - 1);
                    mags.select_nth_unstable_by(k, |a, b| {
                        b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let cutoff = mags[k];
                    merged.retain(|_, &mut v| v.abs() >= cutoff);
                }
            }
        }
        // Rate-based admission: keep only candidates whose leakage rate
        // exceeds `tau_add`. `tau_add = 0` admits everything except exact
        // zeros.
        Ok(merged
            .into_iter()
            .filter(|(_, c)| c.abs() > tau_add)
            .collect())
    }

    /// Sparse generator matrix in COO form: returns `(row, col, val)`
    /// triplets. Row = output Pauli's position in `basis`; col = input
    /// Pauli's position. Output Paulis not in `basis` are silently dropped.
    ///
    /// Precondition: `basis` must not contain duplicate Pauli words
    /// (asserted in debug builds).
    pub fn generator(&self, basis: &[Word]) -> Vec<(usize, usize, f64)> {
        let index = build_basis_index(basis);

        // `compute_action_terms` returns a deduplicated `Vec<(Word, f64)>`,
        // so it can be scattered directly into COO triplets.
        let local: Vec<Vec<(usize, usize, f64)>> = basis
            .par_iter()
            .enumerate()
            .map_init(
                || {
                    (
                        Vec::<u32>::with_capacity(self.n_qubits),
                        Vec::<u32>::with_capacity(128),
                        FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
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
        basis: &[Word],
        coeffs: &[Complex<f64>],
        protected: &[Word],
    ) -> Result<Vec<(Word, Complex<f64>)>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        let in_basis: FxHashMap<u64, ()> =
            basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        const CHUNK_SIZE: usize = 4096;
        let mut merged: FxHashMap<Word, Complex<f64>> = FxHashMap::default();
        for chunk_start in (0..basis.len()).step_by(CHUNK_SIZE) {
            let chunk_end = (chunk_start + CHUNK_SIZE).min(basis.len());
            let chunk_basis = &basis[chunk_start..chunk_end];
            let chunk_coeffs = &coeffs[chunk_start..chunk_end];
            let local: Vec<Vec<(Word, Complex<f64>)>> = chunk_basis
                .par_iter()
                .zip(chunk_coeffs.par_iter())
                .map_init(
                    || {
                        (
                            Vec::<u32>::with_capacity(self.n_qubits),
                            Vec::<u32>::with_capacity(128),
                            FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
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

    /// One predictor-corrector step `O ← exp(dt·L*) O` in the adaptive
    /// real-coefficient Pauli basis: first-hop leakage admission, predictor
    /// exponential, second-hop admission from the predicted state, corrector
    /// exponential from the saved pre-step state, then truncation (prune +
    /// rank cap) per [`PcStepConfig`]. Exact in `dt` within the working
    /// basis — the only error is basis truncation.
    ///
    /// `protected` words are never dropped. All tuning knobs live in `cfg`.
    pub fn pc_step(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        if let Some(n) = cfg.num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| self.pc_step_inner(basis, coeffs, dt, protected, cfg))
        } else {
            self.pc_step_inner(basis, coeffs, dt, protected, cfg)
        }
    }

    /// Same as [`Self::pc_step`] but also returns a per-phase timing
    /// breakdown (microseconds), for profiling parallel scaling and hot
    /// spots.
    pub fn pc_step_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<PcStepTimings, Error> {
        if let Some(n) = cfg.num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| self.pc_step_inner_timed(basis, coeffs, dt, protected, cfg))
        } else {
            self.pc_step_inner_timed(basis, coeffs, dt, protected, cfg)
        }
    }

    fn pc_step_inner_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<PcStepTimings, Error> {
        let PcStepConfig { max_basis, admit_basis, drop_tol, tau_add, .. } = *cfg;
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);
        let mut t = PcStepTimings::default();

        let t0 = Instant::now();
        let leak = self.leakage_with_prune(basis, coeffs, protected, admit, tau_add)?;
        t.leakage1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak, admit);
        t.expand1_us = t0.elapsed().as_micros() as u64;

        // `coeffs` is read-only in the predictor, so it also serves as the
        // pre-step input for the corrector below — no clone needed.
        let t0 = Instant::now();
        let coeffs_predict = self.expm_step(basis, dt, coeffs, drop_tol);
        t.expm1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let leak2 = self.leakage_with_prune(basis, &coeffs_predict, protected, admit, tau_add)?;
        t.leakage2_us = t0.elapsed().as_micros() as u64;
        drop(coeffs_predict);

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak2, admit);
        t.expand2_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        *coeffs = self.expm_step(basis, dt, coeffs, drop_tol);
        t.expm2_us = t0.elapsed().as_micros() as u64;

        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);

        Ok(t)
    }

    fn pc_step_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        protected: &[Word],
        cfg: &PcStepConfig,
    ) -> Result<(), Error> {
        let PcStepConfig { max_basis, admit_basis, drop_tol, tau_add, .. } = *cfg;
        // Admission bound: enrichment may grow the live basis to `admit`
        // >= `max_basis`; the final `cap_basis` then keeps the top-
        // `max_basis` strings by evolved |coeff| over the whole union
        // (retained + admitted) — rank displacement. With `admit_basis =
        // None` admission is bounded by `max_basis` itself, `cap_basis` is
        // a no-op, and membership turnover requires `drop_tol > 0`.
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = tau_add.unwrap_or(0.0);
        // 1. First-hop expansion. After this, `coeffs` contains the pre-step
        // coefficients followed by zeros for the newly-added leakage strings.
        // We rely on `coeffs` itself as the pre-step buffer for the corrector
        // — no `.clone()` is needed because `expm_step` only borrows it.
        let leak = self.leakage_with_prune(basis, coeffs, protected, admit, tau_add)?;
        add_leakage_capped(basis, coeffs, leak, admit);
        // 2. Predictor: `expm_step` reads `coeffs` immutably and returns a
        // new owned vector with the predicted state.
        let coeffs_predict = self.expm_step(basis, dt, coeffs, drop_tol);
        // 3. Second-hop expansion from the predicted state. After leakage2
        // we no longer need `coeffs_predict`. Extend `coeffs` with zeros for
        // any newly-added second-hop strings so it remains a valid input
        // (pre-step state) for the corrector.
        let leak2 = self.leakage_with_prune(basis, &coeffs_predict, protected, admit, tau_add)?;
        drop(coeffs_predict);
        add_leakage_capped(basis, coeffs, leak2, admit);
        // 4. Corrector: redo from pre-step state on the doubly-enlarged basis.
        *coeffs = self.expm_step(basis, dt, coeffs, drop_tol);
        // 5. Prune basis entries below `drop_tol` (protected words never dropped).
        prune_basis(basis, coeffs, drop_tol, protected);
        cap_basis(basis, coeffs, max_basis, protected);
        Ok(())
    }

    /// Compute `exp(dt · M) · b` for the in-basis-restricted generator
    /// `M`, matrix-free, via `quspin-expm` (see [`mf_expm`]).
    fn expm_step(
        &self,
        basis: &[Word],
        dt: f64,
        b: &[f64],
        drop_tol: f64,
    ) -> Vec<f64> {
        mf_expm::expm_apply_mf(self, basis, dt, b, drop_tol)
    }

    /// Compute the unscaled list of `(output, coefficient)` pairs that
    /// `L*(p)` contributes (without the input coefficient).
    pub(crate) fn compute_action_terms(
        &self,
        p: &Word,
        scratch_support: &mut Vec<u32>,
        scratch_cands: &mut Vec<u32>,
        scratch_local: &mut FxHashMap<Word, Complex<f64>>,
    ) -> Vec<(Word, f64)> {
        word_support(p, scratch_support);
        let zero = Complex::new(0.0, 0.0);
        scratch_local.clear();
        let local = scratch_local;

        // ── i [H, p] ─────────────────────────────────────────────────
        candidate_terms(scratch_support, &self.h_support, scratch_cands);
        for &i in scratch_cands.iter() {
            let h = &self.h_terms[i as usize];
            let (r, eps) = comm_product(&h.word, p);
            if eps != 0.0 {
                *local.entry(r).or_insert(zero) += Complex::new(h.coeff * eps, 0.0);
            }
        }

        // ── dissipator ───────────────────────────────────────────────
        candidate_terms(scratch_support, &self.j_support, scratch_cands);
        for &k in scratch_cands.iter() {
            match &self.j_kinds[k as usize] {
                JumpKind::HermitianPauli { word, rate } => {
                    if anti_commutes(word, p) {
                        *local.entry(*p).or_insert(zero) += Complex::new(-2.0 * *rate, 0.0);
                    }
                }
                JumpKind::General {
                    terms,
                    dagger_dagger,
                    rate,
                } => {
                    let rate_c = Complex::new(*rate, 0.0);
                    // Sandwich: γ Σ_{a,b} λ_a* λ_b P_a p P_b.
                    for a in terms {
                        let (r_ap, phi1) = pauli_mul(&a.word, p);
                        for b in terms {
                            let (s, phi2) = pauli_mul(&r_ap, &b.word);
                            let coeff =
                                a.coeff.conj() * b.coeff * phase_factor(phi1 + phi2) * rate_c;
                            *local.entry(s).or_insert(zero) += coeff;
                        }
                    }
                    // -1/2 γ {L†L, p}. For Hermitian Pauli P_c and Pauli p,
                    // {P_c, p} = 2·sign·R if they commute (P_c·p = sign·R),
                    //         = 0          if they anti-commute.
                    for c_term in dagger_dagger {
                        let (r, phase) = pauli_mul(&c_term.word, p);
                        if phase & 1 == 0 {
                            let sign = if phase == 0 { 1.0 } else { -1.0 };
                            let coeff = -c_term.coeff * rate_c * Complex::new(sign, 0.0);
                            *local.entry(r).or_insert(zero) += coeff;
                        }
                    }
                }
            }
        }

        // L* preserves Hermiticity; imaginary parts must cancel to FP noise.
        // `drain()` empties `scratch_local` so its allocation can be reused
        // by the next call on the same thread (`Vec` keeps capacity).
        local
            .drain()
            .filter_map(|(w, c)| {
                debug_assert!(
                    c.im.abs() < 1e-9,
                    "L*(p) produced non-real coefficient {c}; bug in dissipator"
                );
                if c.re == 0.0 { None } else { Some((w, c.re)) }
            })
            .collect()
    }

    /// Accumulate `scale · L*(p)` into `out`.
    fn accumulate_action(
        &self,
        p: &Word,
        scale: f64,
        out: &mut FxHashMap<Word, f64>,
        scratch_support: &mut Vec<u32>,
        scratch_cands: &mut Vec<u32>,
    ) {
        let mut scratch_local = FxHashMap::default();
        let terms =
            self.compute_action_terms(p, scratch_support, scratch_cands, &mut scratch_local);
        for (w, c) in terms.iter() {
            *out.entry(*w).or_insert(0.0) += scale * c;
        }
    }
}

/// Compact 64-bit hash of a [`Word`], used as the key in cache-friendly
/// membership tables: an `FxHashMap<u64, ()>` over the basis has a working
/// set ~6× smaller than `FxHashMap<Word, ()>`. The hash mixes the word's
/// cached hash once through `FxHasher` and never touches the 32-byte
/// payload.
#[inline(always)]
fn word_hash(w: &Word) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = fxhash::FxHasher::default();
    w.hash(&mut h);
    h.finish()
}

/// Compact `basis` / `coeffs` in place: drop entries whose absolute
/// coefficient is below `drop_tol` unless the word appears in `protected`.
/// No-op when `drop_tol ≤ 0`.
fn prune_basis(basis: &mut Vec<Word>, coeffs: &mut Vec<f64>, drop_tol: f64, protected: &[Word]) {
    if drop_tol <= 0.0 {
        return;
    }
    debug_assert_eq!(basis.len(), coeffs.len());
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    let mut write = 0;
    for read in 0..basis.len() {
        if coeffs[read].abs() >= drop_tol || protected_set.contains(&basis[read]) {
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

/// Global max-basis cap (PauliStrings.jl-style top-M trim): keep only the
/// `max_basis` largest-|coeff| terms (protected strings always kept),
/// dropping the rest. Rank-based total-basis bound; dual of `drop_tol`.
/// A `max_basis` large enough to cover the whole basis is a no-op.
fn cap_basis(basis: &mut Vec<Word>, coeffs: &mut Vec<f64>, max_basis: usize, protected: &[Word]) {
    if basis.len() <= max_basis {
        return;
    }
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
    let n_prot = basis.iter().filter(|w| protected_set.contains(w)).count();
    let slots = max_basis.saturating_sub(n_prot);
    let mut mags: Vec<f64> = basis
        .iter()
        .zip(coeffs.iter())
        .filter(|(w, _)| !protected_set.contains(w))
        .map(|(_, c)| c.abs())
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
        if protected_set.contains(&basis[read]) || coeffs[read].abs() >= cutoff {
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

/// Add the largest leakage strings to the basis, up to the available room
/// `room = max_basis − basis.len()` — so the in-step basis (hence the
/// expm/leakage peak memory) never exceeds `max_basis`. New strings get
/// coefficient 0; the surrounding expm fills them. No magnitude filter: the
/// top-`room` by `|leakage|` are added (a large `max_basis` adds them all).
fn add_leakage_capped(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<f64>,
    mut leak: Vec<(Word, f64)>,
    max_basis: usize,
) {
    let room = max_basis.saturating_sub(basis.len());
    if leak.len() > room {
        if room > 0 {
            leak.select_nth_unstable_by(room - 1, |a, b| {
                b.1.abs().partial_cmp(&a.1.abs()).unwrap_or(std::cmp::Ordering::Equal)
            });
        }
        leak.truncate(room);
    }
    for (w, _) in leak {
        basis.push(w);
        coeffs.push(0.0);
    }
}

/// Build a `word → row` map for a basis assumed to contain unique Pauli
/// words; debug-asserts the uniqueness invariant.
pub fn build_basis_index(basis: &[Word]) -> FxHashMap<Word, u32> {
    let mut index: FxHashMap<Word, u32> = FxHashMap::default();
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only full-space complex predictor-corrector step, UNTRUNCATED:
    /// adds every nonzero leakage string (two hops) and applies the exact
    /// in-basis exponential to the complex coefficient vector. Reference
    /// bridge between the real `pc_step` and the orbit-rep path.
    fn pc_step_complex_full(
        spec: &LindbladSpec,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
    ) {
        let protected: Vec<Word> = Vec::new();
        let leak = spec.leakage_complex(basis, coeffs, &protected).unwrap();
        for (w, v) in leak {
            if v.norm() > 0.0 {
                basis.push(w);
                coeffs.push(Complex::new(0.0, 0.0));
            }
        }
        let predict = mf_expm::expm_apply_mf_cxvec(spec, basis, dt, coeffs, 0.0);
        let leak2 = spec.leakage_complex(basis, &predict, &protected).unwrap();
        drop(predict);
        for (w, v) in leak2 {
            if v.norm() > 0.0 {
                basis.push(w);
                coeffs.push(Complex::new(0.0, 0.0));
            }
        }
        *coeffs = mf_expm::expm_apply_mf_cxvec(spec, basis, dt, coeffs, 0.0);
    }

    fn jump_hpauli(s: &str, rate: f64) -> JumpInput {
        JumpInput {
            lincomb: vec![(s.to_string(), Complex::new(1.0, 0.0))],
            rate,
        }
    }

    #[test]
    fn z_dephasing_action_on_x() {
        // L = Z on a single qubit; L*(X) = γ(ZXZ - X) = γ(-X - X) = -2γ X.
        let spec = LindbladSpec::new(
            1,
            &[("X".to_string(), 0.0)], // no Hamiltonian
            &[jump_hpauli("Z", 0.5)],
        )
        .unwrap();
        let (x, _) = parse_pauli_string("X", 1).unwrap();
        let terms = spec.action(&x);
        assert_eq!(terms.len(), 1);
        assert!((terms[0].1 - (-1.0)).abs() < 1e-12); // -2·0.5 = -1
    }

    #[test]
    fn amplitude_damping_action_on_z() {
        // Single-qubit σ⁻ jump: L*(Z) = -γ(I + Z). With γ=1 we expect
        // I coefficient = -1, Z coefficient = -1.
        let sigma_minus = JumpInput {
            lincomb: vec![
                ("X".to_string(), Complex::new(0.5, 0.0)),
                ("Y".to_string(), Complex::new(0.0, -0.5)),
            ],
            rate: 1.0,
        };
        let spec = LindbladSpec::new(1, &[], &[sigma_minus]).unwrap();
        let (z, _) = parse_pauli_string("Z", 1).unwrap();
        let terms = spec.action(&z);
        let (i_word, _) = parse_pauli_string("I", 1).unwrap();
        let mut i_coeff = 0.0;
        let mut z_coeff = 0.0;
        for (w, c) in &terms {
            if w == &i_word {
                i_coeff = *c;
            } else if w == &z {
                z_coeff = *c;
            }
        }
        assert!((i_coeff - (-1.0)).abs() < 1e-10, "I coeff = {i_coeff}");
        assert!((z_coeff - (-1.0)).abs() < 1e-10, "Z coeff = {z_coeff}");
    }

    #[test]
    fn word_codec_roundtrip() {
        let codes = [0u8, 1, 2, 3, 1, 0, 3, 2];
        let w = word_from_codes(&codes).unwrap();
        let mut out = vec![0u8; codes.len()];
        codes_from_word(&w, &mut out);
        assert_eq!(out.as_slice(), &codes);
    }

    /// Per-step orbit-rep evolution gives the SAME final orbit-rep
    /// state as full-basis complex evolution followed by a single
    /// projection at the end. Validates that the phase-aware complex
    /// action machinery is consistent with the full-basis reference.
    #[test]
    fn pc_step_orbit_rep_matches_full_basis_projection() {
        use std::f64::consts::PI;
        use ppvm_pauli_sum::symmetry::canonicalize_pauli_sum_complex;
        let n = 4usize;
        let dt = 0.01f64;
        let n_steps = 3usize;
        let mut h_terms: Vec<(String, f64)> = Vec::new();
        for j in 0..n {
            let nxt = (j + 1) % n;
            for op in ["X", "Y"] {
                let mut s = vec!['I'; n];
                s[j] = op.chars().next().unwrap();
                s[nxt] = op.chars().next().unwrap();
                h_terms.push((s.into_iter().collect(), 1.0));
            }
        }
        let spec = LindbladSpec::new(n, &h_terms, &[]).unwrap();
        let group = ppvm_pauli_sum::symmetry::TranslationGroup::chain_1d(n);
        let k_mode: i32 = 1;
        let k = vec![k_mode];

        // Build the k=1 eigenstate in FULL basis form.
        let basis_full: Vec<Word> = (0..n)
            .map(|j| {
                let mut s = vec!['I'; n];
                s[j] = 'Z';
                let (w, _) = parse_pauli_string(&s.into_iter().collect::<String>(), n).unwrap();
                w
            })
            .collect();
        let coeffs_full: Vec<Complex<f64>> = (0..n as i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * (k_mode as f64) * (a as f64) / (n as f64)))
            .collect();

        // ----- Full-basis path -----
        let mut bf = basis_full.clone();
        let mut cf = coeffs_full.clone();
        let protected: Vec<Word> = Vec::new();
        for _ in 0..n_steps {
            // Full enrichment (tau_add = 0.0 adds every leakage string):
            // for a momentum eigenstate the leakage is pure-sector, so the
            // full-basis and orbit-rep paths build corresponding bases and
            // the projection theorem gives an exact match. The orbit-rep
            // side uses a large max_basis so its rank cap never binds.
            pc_step_complex_full(&spec, &mut bf, &mut cf, dt);
        }
        // Project at the end.
        canonicalize_pauli_sum_complex(&mut bf, &mut cf, &group, &k);

        // ----- Orbit-rep path -----
        // Initial orbit-rep form: project the full-basis input.
        let mut br = basis_full.clone();
        let mut cr = coeffs_full.clone();
        canonicalize_pauli_sum_complex(&mut br, &mut cr, &group, &k);
        // Evolve in orbit-rep form (max_basis large ⇒ full enrichment).
        for _ in 0..n_steps {
            orbit_rep::pc_step_orbit_rep(
                &spec,
                &mut br,
                &mut cr,
                dt,
                &protected,
                &group,
                &k,
                &PcStepConfig { max_basis: 10_000_000, ..Default::default() },
            )
            .unwrap();
        }

        // Compare.
        let mf: FxHashMap<Word, Complex<f64>> = bf.into_iter().zip(cf).collect();
        let mr: FxHashMap<Word, Complex<f64>> = br.into_iter().zip(cr).collect();
        assert_eq!(
            mf.len(),
            mr.len(),
            "orbit-rep ({}) and full-basis-projected ({}) basis sizes differ",
            mr.len(),
            mf.len()
        );
        let mut max_diff = 0.0_f64;
        for (w, cm) in &mr {
            let cf_val = mf
                .get(w)
                .copied()
                .unwrap_or_else(|| panic!("rep {:?} in orbit-rep but not in full-basis", w));
            max_diff = max_diff.max((cm - cf_val).norm());
        }
        assert!(
            max_diff < 1e-9,
            "orbit-rep diverged from full-basis: max |Δc| = {max_diff:e}"
        );
    }

    /// The full-space complex step at momentum k=0 must reproduce the real
    /// pc_step on the same trajectory exactly.
    #[test]
    fn complex_full_matches_real_at_kzero() {
        let n = 4usize;
        let dt = 0.01f64;
        let n_steps = 5usize;
        let mut h_terms: Vec<(String, f64)> = Vec::new();
        for j in 0..n {
            let nxt = (j + 1) % n;
            for op in ["X", "Y"] {
                let mut s = vec!['I'; n];
                s[j] = op.chars().next().unwrap();
                s[nxt] = op.chars().next().unwrap();
                h_terms.push((s.into_iter().collect(), 1.0));
            }
        }
        let spec = LindbladSpec::new(n, &h_terms, &[]).unwrap();

        let mut basis_r: Vec<Word> = (0..n)
            .map(|j| {
                let mut s = vec!['I'; n];
                s[j] = 'Z';
                let st: String = s.into_iter().collect();
                let (w, _) = parse_pauli_string(&st, n).unwrap();
                w
            })
            .collect();
        let mut coeffs_r: Vec<f64> = vec![1.0; n];

        let mut basis_c = basis_r.clone();
        let mut coeffs_c: Vec<Complex<f64>> =
            coeffs_r.iter().map(|&v| Complex::new(v, 0.0)).collect();

        let protected: Vec<Word> = Vec::new();
        for _ in 0..n_steps {
            // Large max_basis: rank cap never binds, so the real path
            // enriches fully (adds every leakage string). Match the
            // complex path by setting its tau_add=0.0 (also full
            // enrichment) so the two stay in lock-step at k=0.
            spec.pc_step(
                &mut basis_r,
                &mut coeffs_r,
                dt,
                &protected,
                &PcStepConfig { max_basis: 10_000_000, ..Default::default() },
            )
            .unwrap();
            pc_step_complex_full(&spec, &mut basis_c, &mut coeffs_c, dt);
        }
        // Match as (word → coeff) maps.
        let map_r: FxHashMap<Word, f64> =
            basis_r.into_iter().zip(coeffs_r).collect();
        let map_c: FxHashMap<Word, Complex<f64>> =
            basis_c.into_iter().zip(coeffs_c).collect();
        assert_eq!(
            map_r.len(),
            map_c.len(),
            "real and complex pc_step produced different basis sizes ({} vs {})",
            map_r.len(),
            map_c.len()
        );
        let mut max_diff = 0.0_f64;
        for (w, cr) in &map_r {
            let cc = map_c
                .get(w)
                .copied()
                .unwrap_or_else(|| panic!("word {:?} in real but not complex", w));
            assert!(cc.im.abs() < 1e-10, "expected zero imag at k=0, got {cc:?}");
            max_diff = max_diff.max((cr - cc.re).abs());
        }
        assert!(
            max_diff < 1e-10,
            "real vs complex pc_step diverged: max |Δc| = {max_diff:e}"
        );
    }

    /// Small-system end-to-end check that orbit-rep merging gives the
    /// same physics as standard evolution, when no truncation is applied.
    ///
    /// Setup: n=4 qubit chain, PBC, translation-invariant XY Hamiltonian
    /// `H = Σ_j (X_j X_{j+1} + Y_j Y_{j+1})`, no dissipation. Initial
    /// operator `O(0) = Σ_j Z_j` is translation-invariant (k=0 sector).
    ///
    /// Run 10 pc_step iterations with `drop_tol = 0` (no truncation):
    /// once without merging, once applying `canonicalize_pauli_sum`
    /// after each step. Canonicalize the un-merged final state once at
    /// the end. The two orbit-rep representations should be
    /// bit-identical up to FP noise.
    #[test]
    fn pc_step_matches_symmetry_merged_on_small_chain() {
        use ppvm_pauli_sum::symmetry::{TranslationGroup, canonicalize_pauli_sum};

        let n = 4usize;
        let dt = 0.05f64;
        let n_steps = 10usize;

        // Build XY-chain Hamiltonian with PBC. 8 terms (4 bonds × {XX, YY}).
        let mut h_terms: Vec<(String, f64)> = Vec::new();
        for j in 0..n {
            let nxt = (j + 1) % n;
            for op in ["X", "Y"] {
                let mut s = vec!['I'; n];
                s[j] = op.chars().next().unwrap();
                s[nxt] = op.chars().next().unwrap();
                h_terms.push((s.into_iter().collect(), 1.0));
            }
        }
        // No dissipation.
        let spec = LindbladSpec::new(n, &h_terms, &[]).unwrap();
        let group = TranslationGroup::chain_1d(n);

        // Initial: O(0) = Σ_j Z_j (translation-invariant).
        let mut basis_u: Vec<Word> = (0..n)
            .map(|j| {
                let mut s = vec!['I'; n];
                s[j] = 'Z';
                let st: String = s.into_iter().collect();
                let (w, _) = parse_pauli_string(&st, n).unwrap();
                w
            })
            .collect();
        let mut coeffs_u: Vec<f64> = vec![1.0; n];

        // Mirror state for the "with merging" run.
        let mut basis_m = basis_u.clone();
        let mut coeffs_m = coeffs_u.clone();

        let protected: Vec<Word> = Vec::new();
        for _ in 0..n_steps {
            // max_basis == current basis size → room = 0: no leakage
            // enrichment, only the expm step (the regime where merging
            // commutes with evolution). drop_tol = 0 → no truncation.
            let cfg_u = PcStepConfig { max_basis: basis_u.len(), ..Default::default() };
            spec.pc_step(&mut basis_u, &mut coeffs_u, dt, &protected, &cfg_u)
                .unwrap();

            let cfg_m = PcStepConfig { max_basis: basis_m.len(), ..Default::default() };
            spec.pc_step(&mut basis_m, &mut coeffs_m, dt, &protected, &cfg_m)
                .unwrap();
            // Apply symmetry merging on the "with merging" run only.
            canonicalize_pauli_sum(&mut basis_m, &mut coeffs_m, &group);
        }

        // Canonicalize the un-merged final state once.
        canonicalize_pauli_sum(&mut basis_u, &mut coeffs_u, &group);

        // Both representations should now be in orbit-rep form; compare
        // as (word → coeff) maps with FP tolerance.
        let map_u: FxHashMap<Word, f64> = basis_u.into_iter().zip(coeffs_u).collect();
        let map_m: FxHashMap<Word, f64> = basis_m.into_iter().zip(coeffs_m).collect();
        assert_eq!(
            map_u.len(),
            map_m.len(),
            "merged basis size {} != post-merged-unmerged basis size {}",
            map_m.len(),
            map_u.len()
        );
        let mut max_diff = 0.0f64;
        for (w, c_u) in &map_u {
            let c_m = map_m.get(w).copied().unwrap_or_else(|| {
                panic!(
                    "rep {:?} present in un-merged-then-canonicalized but not in merged",
                    w
                );
            });
            max_diff = max_diff.max((c_u - c_m).abs());
        }
        assert!(
            max_diff < 1e-9,
            "with-merging vs without-merging diverged: max |Δc| = {max_diff:e}"
        );
    }
}
