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
//! Pauli strings are stored as [`ppvm_runtime::word::PauliWord`] backed by
//! `[u64; 2]` (≤128 qubits) with cached hashes for fast HashMap lookup. The
//! hot-path commutator/product loops bypass the higher-level word API and
//! operate directly on raw `u64` chunks for speed.

use fxhash::{FxBuildHasher, FxHashMap, FxHashSet};
use num::Complex;
use ppvm_runtime::traits::PauliWordTrait;
use ppvm_runtime::word::PauliWord;
use rayon::prelude::*;
use std::fmt;
use std::time::Instant;

pub mod expm;
pub use expm::{Csr, ExpmOpts, csr_from_triplets, csr_one_norm, expm_multiply, spmv_parallel};

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

/// Errors raised when constructing a [`LindbladSpec`].
#[derive(Debug, Clone)]
pub enum Error {
    TooManyQubits {
        got: usize,
    },
    LengthMismatch {
        what: &'static str,
        a: usize,
        b: usize,
    },
    InvalidPauliCode {
        code: u8,
    },
    InvalidPauliChar {
        c: char,
    },
    WrongLength {
        expected: usize,
        got: usize,
    },
    NegativeRate {
        index: usize,
        rate: f64,
    },
    EmptyLincomb {
        index: usize,
    },
    Internal(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::TooManyQubits { got } => {
                write!(
                    f,
                    "LindbladSpec supports n_qubits ≤ {MAX_QUBITS}; got {got}"
                )
            }
            Error::LengthMismatch { what, a, b } => {
                write!(f, "{what}: expected matching lengths, got {a} and {b}")
            }
            Error::InvalidPauliCode { code } => write!(
                f,
                "Pauli code must be 0 (I), 1 (X), 2 (Z), or 3 (Y); got {code}"
            ),
            Error::InvalidPauliChar { c } => {
                write!(f, "invalid Pauli character '{c}'; expected I, X, Y, or Z")
            }
            Error::WrongLength { expected, got } => {
                write!(f, "Pauli string has length {got} but n_qubits = {expected}")
            }
            Error::NegativeRate { index, rate } => {
                write!(f, "jump rate must be non-negative; got γ_{index} = {rate}")
            }
            Error::EmptyLincomb { index } => {
                write!(
                    f,
                    "jump {index}: lincomb must contain at least one Pauli term"
                )
            }
            Error::Internal(msg) => write!(f, "internal error: {msg}"),
        }
    }
}

impl std::error::Error for Error {}

/// `Send + Sync` wrapper around `*mut T` for parallel scatter writes from
/// rayon tasks into pre-allocated arrays. Caller must guarantee that the
/// writes target disjoint indices.
///
/// The field is intentionally private so closures that use [`Self::write`]
/// capture the whole struct (which is `Send + Sync`) rather than the
/// inner `*mut T` (which is not) — Rust 2021's disjoint-field capture
/// would otherwise peel the wrapper.
#[derive(Clone, Copy)]
struct SendPtr<T>(*mut T);
unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}

impl<T> SendPtr<T> {
    /// SAFETY: caller guarantees `offset` is in bounds and that no other
    /// thread writes to the same offset.
    #[inline(always)]
    unsafe fn write(self, offset: usize, value: T) {
        unsafe { self.0.add(offset).write(value) }
    }
}

/// Per-phase timing breakdown (microseconds) returned by
/// [`LindbladSpec::pc_step_timed`].
#[derive(Default, Clone, Copy, Debug)]
pub struct PcStepTimings {
    pub leakage1_us: u64,
    pub expand1_us: u64,
    pub gencsr1_us: u64,
    pub expm1_us: u64,
    pub leakage2_us: u64,
    pub expand2_us: u64,
    pub gencsr2_us: u64,
    pub expm2_us: u64,
}

impl PcStepTimings {
    pub fn total_us(&self) -> u64 {
        self.leakage1_us
            + self.expand1_us
            + self.gencsr1_us
            + self.expm1_us
            + self.leakage2_us
            + self.expand2_us
            + self.gencsr2_us
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

/// Fused product `h · p`: returns `(out, eps)` where `out = h ⊕ p` and
///
/// - `eps =  0` if `h` and `p` commute (caller should skip — `[h,p] = 0`),
/// - `eps = -2.0` if `h·p` has phase `+i` (so `i·[h,p] = -2·out`),
/// - `eps = +2.0` if `h·p` has phase `-i` (so `i·[h,p] = +2·out`).
#[inline(always)]
fn comm_product(h: &Word, p: &Word) -> (Word, f64) {
    let mut out = Word::new(h.n_qubits());
    let mut sign_count: u32 = 0;
    let mut imag_count: u32 = 0;
    for i in 0..W_U64 {
        let a = h.xbits.data[i];
        let b = h.zbits.data[i];
        let c = p.xbits.data[i];
        let d = p.zbits.data[i];
        let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
        let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
        sign_count += sign.count_ones();
        imag_count += imag.count_ones();
        out.xbits.data[i] = a ^ c;
        out.zbits.data[i] = b ^ d;
    }
    out.rehash();
    let phase = (2 * sign_count + imag_count) & 3;
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
/// [`Self::leakage`], [`Self::generator`]. `L*(p)` is recomputed fresh on
/// every call — empirical benchmarks showed that the previous global
/// `action_cache` hurt wall time for sparse-local Hamiltonians (hash-map
/// lookup ≳ recompute) and consumed several KB per cached Pauli word, which
/// blocked us from reaching the basis sizes needed for L=41 sweeps.
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
        // would spill to DRAM. Each query then costs ~30 ns instead of
        // ~500 ns. False-positive risk from a 64-bit hash collision is
        // ≈ N/2⁶⁴ per lookup (e.g. ≈ 3·10⁻¹⁴ at N = 5·10⁵); the worst-case
        // cost is occasionally dropping a single leakage entry that should
        // have been kept, which is negligible against `tau_add` and other
        // approximations.
        let in_basis: FxHashMap<u64, ()> = basis.iter().map(|w| (word_hash(w), ())).collect();
        let protected_set: FxHashMap<u64, ()> =
            protected.iter().map(|w| (word_hash(w), ())).collect();

        // `map_init` gives each rayon worker a thread-local scratch tuple
        // that lives for the duration of the par_iter. We reuse the
        // hashmap + scratch vectors across all basis elements processed
        // by the same thread — the previous per-task allocations were a
        // real cost at large basis sizes.
        let local: Vec<Vec<(Word, f64)>> = basis
            .par_iter()
            .zip(coeffs.par_iter())
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
                            out.push((w.clone(), c * *v));
                        }
                    }
                    out
                },
            )
            .collect();

        let mut merged: FxHashMap<Word, f64> = FxHashMap::default();
        for v in local {
            for (k, val) in v {
                *merged.entry(k).or_insert(0.0) += val;
            }
        }
        Ok(merged.into_iter().filter(|(_, c)| *c != 0.0).collect())
    }

    /// Sparse generator matrix in COO form: returns `(row, col, val)`
    /// triplets. Row = output Pauli's position in `basis`; col = input
    /// Pauli's position. Output Paulis not in `basis` are silently dropped.
    ///
    /// Precondition: `basis` must not contain duplicate Pauli words
    /// (asserted in debug builds).
    pub fn generator(&self, basis: &[Word]) -> Vec<(usize, usize, f64)> {
        let index = build_basis_index(basis);

        // The cached `L*(p)` is already a deduplicated `Vec<(Word, f64)>`,
        // so we can iterate it directly without going through a per-task
        // `FxHashMap` accumulator (which was the previous hot spot).
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

    /// Build the generator restricted to `basis` directly in CSR form,
    /// skipping the intermediate `(row, col, val)` triplet materialisation.
    ///
    /// Algorithm: each task computes its column's contributions and writes
    /// them into a thread-local `(col_idx, values)` buffer plus a per-row
    /// count. Counts are then prefix-summed across threads into the final
    /// CSR `row_ptr`, and a second parallel pass scatters each task's data
    /// into the right CSR positions. This avoids both the
    /// `Vec<Vec<...>>` → `Vec<...>` flatten *and* the sequential
    /// `from_triplets` count-and-scatter that dominated the original
    /// `generator_csr` (~40% of PC-step wall time).
    pub fn generator_csr(&self, basis: &[Word]) -> Csr {
        let n = basis.len();
        if n == 0 {
            return Csr::new((0, 0), vec![0], Vec::new(), Vec::new());
        }

        let index = build_basis_index(basis);

        // Phase 1: per-column, collect `(row, value)` pairs (already
        // filtered to in-basis), in parallel with thread-local scratch.
        let cols: Vec<Vec<(u32, f64)>> = basis
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
                |(s1, s2, lm), p| {
                    let terms = self.compute_action_terms(p, s1, s2, lm);
                    let mut out = Vec::with_capacity(terms.len());
                    for (w, v) in terms.iter() {
                        if let Some(&row) = index.get(w) {
                            out.push((row, *v));
                        }
                    }
                    out
                },
            )
            .collect();

        // Phase 2: per-row counts. Each task scans its column and bumps
        // an atomic per-row counter; total cost is `O(nnz / threads)`.
        use std::sync::atomic::{AtomicU32, Ordering};
        let row_counts: Vec<AtomicU32> = (0..n).map(|_| AtomicU32::new(0)).collect();
        cols.par_iter().for_each(|col_data| {
            for &(row, _) in col_data.iter() {
                row_counts[row as usize].fetch_add(1, Ordering::Relaxed);
            }
        });

        // Phase 3: prefix-sum the counts into `row_ptr`. Sequential but
        // tiny (`O(n)`).
        let mut row_ptr = vec![0usize; n + 1];
        for i in 0..n {
            let c = row_counts[i].load(Ordering::Relaxed) as usize;
            row_ptr[i + 1] = row_ptr[i] + c;
        }
        let nnz = row_ptr[n];

        // Phase 4: scatter. Reuse `row_counts` as a per-row write offset
        // (atomic fetch-and-add); each task writes its column's entries
        // into the row's allocated slot. The increments target distinct
        // rows mostly, and even when they collide, atomics on `Relaxed`
        // are cheap on modern x86 / ARM.
        for c in row_counts.iter() {
            c.store(0, Ordering::Relaxed);
        }
        // We need parallel writes into shared `indices` and `data`.
        // Safety: each `(row, slot)` pair is unique by construction (every
        // (row, col) appears at most once per column, and atomic increment
        // gives a unique slot per (row, write) event). We use raw-pointer
        // writes inside `unsafe` to express this without RwLock overhead.
        let mut indices = vec![0u32; nnz];
        let mut data = vec![0f64; nnz];
        let indices_ptr = SendPtr(indices.as_mut_ptr());
        let data_ptr = SendPtr(data.as_mut_ptr());
        let row_ptr_ref: &[usize] = &row_ptr;
        let row_counts_ref: &[AtomicU32] = &row_counts;
        // `move` on the closure forces whole-struct capture of `SendPtr`
        // (rust 2021 edition does field-level capture by default, which
        // would peel the inner `*mut T` and leak the !Sync). All captures
        // are Copy (SendPtr) or `&` (refs), so `Fn` is still satisfied.
        cols.par_iter()
            .enumerate()
            .for_each(move |(col, col_data)| {
                for &(row, v) in col_data.iter() {
                    let slot_offset = row_counts_ref[row as usize].fetch_add(1, Ordering::Relaxed);
                    let pos = row_ptr_ref[row as usize] + slot_offset as usize;
                    // SAFETY: each `pos` is unique (atomic gives a unique slot
                    // per (row, write event)) and `pos < nnz` (counter bounded
                    // by precomputed row count).
                    unsafe {
                        indices_ptr.write(pos, col as u32);
                        data_ptr.write(pos, v);
                    }
                }
            });

        Csr::new_from_unsorted((n, n), row_ptr, indices, data)
            .map_err(|(_, _, _, e)| e)
            .expect("invalid CSR structure")
    }

    /// Predictor-corrector adaptive step.
    ///
    /// Mutates `basis` (may grow) and `coeffs` in place to reflect the
    /// state at time `t + dt`. The step is:
    ///
    /// 1. **First-hop expansion**: compute leakage from the current state
    ///    and append any leakage Pauli with `|coeff| > tau_add` to the
    ///    basis (with starting coefficient 0).
    /// 2. **Predictor**: apply `exp(dt · M)` to `coeffs` on the enlarged
    ///    basis, yielding a predicted state.
    /// 3. **Second-hop expansion**: compute leakage from the *predicted*
    ///    state and append any further leakage strings — these are the
    ///    second-hop Paulis the predictor flowed into but did not yet have
    ///    in basis.
    /// 4. **Corrector**: redo `exp(dt · M)` on the doubly-enlarged basis
    ///    starting from the saved pre-step coefficients.
    ///
    /// Lifts the per-step truncation error from `O(dt²)` (single-hop) to
    /// `O(dt³)`. Strings in `protected` are never added to the basis as
    /// leakage targets — typically the observable's support, which the
    /// caller wants tracked exactly.
    ///
    /// `num_threads`, when set, runs the entire step inside a freshly built
    /// rayon thread pool of that size — useful for benchmarking parallel
    /// scaling. When `None`, the global rayon pool is used.
    #[allow(clippy::too_many_arguments)]
    pub fn pc_step(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        opts: ExpmOpts,
        num_threads: Option<usize>,
    ) -> Result<(), Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| self.pc_step_inner(basis, coeffs, dt, tau_add, drop_tol, protected, opts))
        } else {
            self.pc_step_inner(basis, coeffs, dt, tau_add, drop_tol, protected, opts)
        }
    }

    /// Same as [`Self::pc_step`] but also returns a per-phase timing
    /// breakdown (microseconds), for profiling parallel scaling and hot
    /// spots. Output: `(leakage1, expand1, gencsr1, expm1, leakage2,
    /// expand2, gencsr2, expm2)`.
    #[allow(clippy::too_many_arguments)]
    pub fn pc_step_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        opts: ExpmOpts,
        num_threads: Option<usize>,
    ) -> Result<PcStepTimings, Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| self.pc_step_inner_timed(basis, coeffs, dt, tau_add, drop_tol, protected, opts))
        } else {
            self.pc_step_inner_timed(basis, coeffs, dt, tau_add, drop_tol, protected, opts)
        }
    }

    fn pc_step_inner_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        opts: ExpmOpts,
    ) -> Result<PcStepTimings, Error> {
        let mut t = PcStepTimings::default();

        let t0 = Instant::now();
        let leak = self.leakage(basis, coeffs, protected)?;
        t.leakage1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        for (w, v) in leak {
            if v.abs() > tau_add {
                basis.push(w);
                coeffs.push(0.0);
            }
        }
        t.expand1_us = t0.elapsed().as_micros() as u64;

        let coeffs_pre = coeffs.clone();

        let t0 = Instant::now();
        let csr = self.generator_csr(basis);
        t.gencsr1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let coeffs_predict = expm_multiply(&csr, dt, coeffs, opts);
        t.expm1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let leak2 = self.leakage(basis, &coeffs_predict, protected)?;
        t.leakage2_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let mut coeffs_pre_padded = coeffs_pre;
        for (w, v) in leak2 {
            if v.abs() > tau_add {
                basis.push(w);
                coeffs_pre_padded.push(0.0);
            }
        }
        t.expand2_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let csr = self.generator_csr(basis);
        t.gencsr2_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        *coeffs = expm_multiply(&csr, dt, &coeffs_pre_padded, opts);
        t.expm2_us = t0.elapsed().as_micros() as u64;

        prune_basis(basis, coeffs, drop_tol, protected);

        Ok(t)
    }

    fn pc_step_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        opts: ExpmOpts,
    ) -> Result<(), Error> {
        // 1. First-hop expansion.
        let leak = self.leakage(basis, coeffs, protected)?;
        for (w, v) in leak {
            if v.abs() > tau_add {
                basis.push(w);
                coeffs.push(0.0);
            }
        }
        // 2. Predictor.
        let coeffs_pre = coeffs.clone();
        let csr = self.generator_csr(basis);
        let coeffs_predict = expm_multiply(&csr, dt, coeffs, opts);
        // 3. Second-hop expansion from the predicted state.
        let leak2 = self.leakage(basis, &coeffs_predict, protected)?;
        let mut coeffs_pre_padded = coeffs_pre;
        for (w, v) in leak2 {
            if v.abs() > tau_add {
                basis.push(w);
                coeffs_pre_padded.push(0.0);
            }
        }
        // 4. Corrector: redo from pre-step state on the doubly-enlarged basis.
        let csr = self.generator_csr(basis);
        *coeffs = expm_multiply(&csr, dt, &coeffs_pre_padded, opts);
        // 5. Prune basis entries below `drop_tol` (protected words never dropped).
        prune_basis(basis, coeffs, drop_tol, protected);
        Ok(())
    }

    /// Compute the unscaled list of `(output, coefficient)` pairs that
    /// `L*(p)` contributes (without the input coefficient).
    fn compute_action_terms(
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
                        *local.entry(p.clone()).or_insert(zero) += Complex::new(-2.0 * *rate, 0.0);
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
            *out.entry(w.clone()).or_insert(0.0) += scale * c;
        }
    }
}

/// Compact 64-bit hash of a [`Word`], used as the key in cache-friendly
/// membership tables.
///
/// `PauliWord::hash` writes its cached `u64` into the supplied `Hasher`;
/// running it through `FxHasher` once mixes the bits enough that the
/// resulting `u64` is well-distributed when stored as a key in a downstream
/// `FxHashMap` (whose outer hash function then needs only a cheap multiply
/// on a `u64` key). The whole call is ~2-5 ns and never touches `Word`'s
/// 32-byte payload, which is the entire point: an `FxHashMap<u64, ()>`
/// over the basis has a working set ~6× smaller than `FxHashMap<Word, ()>`
/// and stays in L2/L3 well past basis 10⁶.
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

/// Build a `word → row` map for a basis assumed to contain unique Pauli
/// words; debug-asserts the uniqueness invariant.
fn build_basis_index(basis: &[Word]) -> FxHashMap<Word, u32> {
    let mut index: FxHashMap<Word, u32> = FxHashMap::default();
    for (i, w) in basis.iter().enumerate() {
        let prev = index.insert(w.clone(), i as u32);
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
}
