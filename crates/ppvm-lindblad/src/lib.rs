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
use ppvm_pauli_sum::symmetry::TranslationGroup;
use ppvm_traits::PauliWordTrait;
use ppvm_pauli_word::word::PauliWord;
use rayon::prelude::*;
use std::fmt;
use std::time::Instant;

pub mod expm;

/// Matrix-free / quspin-expm-backed `exp(dt·L*)·b` engine. See module docs.
pub(crate) mod mf_expm;

/// Per-step orbit-rep evolution under translation symmetry — phase-aware
/// action + complex CSR throughout. See module docs.
pub mod orbit_rep;

/// Words pack up to 128 qubits.
const W_U64: usize = 2;

/// Predictor-corrector leakage-admission shape factor `K` (dimensionless).
///
/// A leakage string is admitted only if its rate `|ell_q|` exceeds
/// `tau_add = K * drop_tol / dt`, i.e. if its acquired coefficient
/// `|ell_q|*dt` would exceed `K * drop_tol` — comfortably above the post-step
/// prune threshold. This stops the doubly-enriched transient from admitting the
/// (~90%+) strings that would be pruned anyway.
///
/// Read once from `PPVM_K_LEAKAGE` (default `0.0`). `K = 0` ⇒ `tau_add = 0` ⇒
/// admit every leakage string (the historical behaviour, and the behaviour at
/// `drop_tol = 0`), so the default changes nothing and exact-reference tests are
/// unaffected. Set `PPVM_K_LEAKAGE` to sweep `K` without rebuilding.
fn k_leakage() -> f64 {
    use std::sync::OnceLock;
    static K: OnceLock<f64> = OnceLock::new();
    *K.get_or_init(|| {
        std::env::var("PPVM_K_LEAKAGE")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|k| k.is_finite() && *k >= 0.0)
            .unwrap_or(0.0)
    })
}

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

    /// Conservative upper bound on `max_{p, w} |⟨w| L*(p) |p⟩|`, the
    /// largest coefficient any single input Pauli can produce in any
    /// single output Pauli through one action of `L*`. Useful for
    /// Cauchy-Schwarz style streaming bounds on accumulated leakage.
    ///
    /// The bound is `2 · max |h_i|` (commutator factor) plus a per-jump
    /// term. For the common case of Hermitian-Pauli H and Hermitian-Pauli
    /// jumps this collapses to `2 · max(|h_i|, γ_k)`.
    pub fn max_action_coef(&self) -> f64 {
        let mut m = 0.0_f64;
        for h in &self.h_terms {
            // `i [H_i, p]` has coefficient `±2 i · h_i · eps`; magnitude ≤ 2|h_i|.
            m = m.max(2.0 * h.coeff.abs());
        }
        for jk in &self.j_kinds {
            match jk {
                JumpKind::HermitianPauli { rate, .. } => {
                    // sandwich + anticommutator → max |2γ| on anticommuting p
                    m = m.max(2.0 * rate.abs());
                }
                JumpKind::General {
                    terms,
                    rate,
                    dagger_dagger,
                } => {
                    let lam_sum: f64 = terms.iter().map(|t| t.coeff.norm()).sum();
                    let dd_sum: f64 = dagger_dagger.iter().map(|t| t.coeff.norm()).sum();
                    // Sandwich `γ Σ_{a,b} λ_a* λ_b P_a p P_b`: per-output ≤ γ · (Σ|λ|)²
                    // L†L anticommutator: per-output ≤ γ · Σ|c_dd|
                    m = m.max(rate * (lam_sum * lam_sum + dd_sum));
                }
            }
        }
        m
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
        // Cap the live leakage map to the *available room* `max_basis −
        // basis.len()` — only the strings we could actually add are worth
        // keeping. When `room ≥ all candidates`, nothing is dropped (the
        // near-exact case). `usize::MAX` room ⇒ no cap.
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
                                out.push((w.clone(), c * *v));
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
        // exceeds `tau_add` (= K*drop_tol/dt). `tau_add = 0` admits everything
        // except exact zeros, matching the prior behaviour.
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


    /// Matrix-free SpMV: `y ← M · x` where `M` is the in-basis-restricted
    /// generator, never materialised. For each basis column `j` with `x[j] != 0`,
    /// compute `L*(basis[j])`, look up each output Pauli in `basis_index`,
    /// and accumulate `v · x[j]` into `y` at the matching row.
    ///
    /// `basis_index` must be the `Word → row` map for `basis` (use
    /// [`build_basis_index`]). It is the caller's responsibility to build
    /// it once and reuse across all SpMVs within an expm call.
    ///
    /// Threading uses `rayon::current_num_threads()` per-task dense
    /// accumulators, reduced into `y` at the end. Peak transient memory is
    /// roughly `T × n × 8 B` where `T` = thread count, `n` = basis size.
    pub fn spmv_matrix_free(
        &self,
        basis: &[Word],
        basis_index: &FxHashMap<Word, u32>,
        x: &[f64],
        y: &mut [f64],
    ) {
        debug_assert_eq!(basis.len(), x.len());
        debug_assert_eq!(basis.len(), y.len());
        let n = basis.len();
        if n == 0 {
            return;
        }
        let num_threads = rayon::current_num_threads().max(1);
        let chunk_size = n.div_ceil(num_threads);

        let partial_ys: Vec<Vec<f64>> = basis
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let c_offset = chunk_idx * chunk_size;
                let mut y_local = vec![0.0; n];
                let mut s1 = Vec::<u32>::with_capacity(self.n_qubits);
                let mut s2 = Vec::<u32>::with_capacity(128);
                let mut lm = FxHashMap::<Word, Complex<f64>>::with_capacity_and_hasher(
                    128,
                    FxBuildHasher::default(),
                );
                for (c_local, p) in chunk.iter().enumerate() {
                    let c = c_offset + c_local;
                    let xc = x[c];
                    if xc == 0.0 {
                        continue;
                    }
                    let terms = self.compute_action_terms(p, &mut s1, &mut s2, &mut lm);
                    for (w, v) in terms.iter() {
                        if let Some(&row) = basis_index.get(w) {
                            y_local[row as usize] += *v * xc;
                        }
                    }
                }
                y_local
            })
            .collect();

        // Sequential reduce. T × n adds; at T=8, n=10⁶ this is ~8M adds
        // (~10 ms), trivial vs the action evaluations above.
        y.fill(0.0);
        for partial in &partial_ys {
            for (yi, &pi) in y.iter_mut().zip(partial.iter()) {
                *yi += pi;
            }
        }
    }


    /// Apply `L*` to the Pauli sum `Σ_j coeffs[j] · basis[j]` and return
    /// the result as a (Word → real coefficient) map. The basis of the
    /// returned operator is determined by the action — entries appear
    /// for every Pauli that `L*` reaches with nonzero coefficient.
    ///
    /// Same parallel structure and chunked merge as [`Self::leakage`],
    /// without the in-basis / protected filter: this returns the full
    /// `L*(O)`, not just its off-basis component.
    pub fn compute_action_sum(
        &self,
        basis: &[Word],
        coeffs: &[f64],
    ) -> Result<FxHashMap<Word, f64>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
        const CHUNK_SIZE: usize = 4096;
        let mut merged: FxHashMap<Word, f64> = FxHashMap::default();
        for chunk_start in (0..basis.len()).step_by(CHUNK_SIZE) {
            let chunk_end = (chunk_start + CHUNK_SIZE).min(basis.len());
            let chunk_basis = &basis[chunk_start..chunk_end];
            let chunk_coeffs = &coeffs[chunk_start..chunk_end];
            let local: Vec<Vec<(Word, f64)>> = chunk_basis
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
                        terms.into_iter().map(|(w, v)| (w, c * v)).collect()
                    },
                )
                .collect();
            for v in local {
                for (k, val) in v {
                    *merged.entry(k).or_insert(0.0) += val;
                }
            }
        }
        Ok(merged)
    }

    /// Complex-coefficient variant of [`Self::compute_action_sum`]:
    /// `L*( Σ_j coeffs[j] · basis[j] )` with complex coefficients.
    /// Matrix elements of `L*` are real (Hermiticity-preserving), so
    /// each complex coefficient just multiplies into the real action
    /// outputs.
    pub fn compute_action_sum_complex(
        &self,
        basis: &[Word],
        coeffs: &[Complex<f64>],
    ) -> Result<FxHashMap<Word, Complex<f64>>, Error> {
        if basis.len() != coeffs.len() {
            return Err(Error::LengthMismatch {
                what: "basis and coeffs",
                a: basis.len(),
                b: coeffs.len(),
            });
        }
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
                        terms
                            .into_iter()
                            .map(|(w, v)| (w, c * v))
                            .collect()
                    },
                )
                .collect();
            for v in local {
                for (k, val) in v {
                    *merged.entry(k).or_insert(Complex::new(0.0, 0.0)) += val;
                }
            }
        }
        Ok(merged)
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
                                out.push((w.clone(), c * *v));
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

    /// Complex-coefficient predictor-corrector step.
    ///
    /// Same algorithm as [`Self::pc_step`]: first-hop leakage expansion,
    /// predictor `exp(dt·M)`, second-hop expansion from predicted state,
    /// corrector `exp(dt·M)` from saved pre-step state, then `drop_tol`
    /// prune. Coefficients are complex throughout; `L*` matrix elements
    /// are real, so the real matrix-free engine is applied to the real and
    /// imaginary parts separately (see [`Self::expm_step_complex`]).
    ///
    /// `drop_tol` is compared against `|c|` (complex magnitude). The
    /// `protected` set is honoured exactly as in the real path.
    ///
    /// When `sym_group` and `momentum` are both `Some`, the input
    /// `(basis, coeffs)` is checked against the momentum sector using
    /// [`ppvm_pauli_sum::symmetry::check_momentum_sector`] before any
    /// evolution. If the input is not in sector `momentum`, an error
    /// is returned and no evolution is performed. This is purely
    /// validation — pc_step itself evolves in the standard full-basis
    /// representation, and translation symmetry of `L*` keeps the
    /// state inside the input's momentum sector throughout.
    ///
    /// **For memory reduction**: this function does NOT merge into
    /// orbit-rep form during evolution. The correct workflow is to
    /// evolve as many steps as desired in full-basis form, then
    /// call [`ppvm_pauli_sum::symmetry::canonicalize_pauli_sum_complex`]
    /// at observable-readout / snapshot points. Per-step orbit-rep
    /// evolution would require phase-aware action / CSR (complex matrix
    /// elements, not just complex vectors) — a separate code path to
    /// be added later.
    #[allow(clippy::too_many_arguments)]
    pub fn pc_step_complex(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        num_threads: Option<usize>,
        sym_group: Option<&TranslationGroup>,
        momentum: Option<&[i32]>,
    ) -> Result<(), Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| {
                self.pc_step_complex_inner(
                    basis, coeffs, dt, tau_add, drop_tol, protected, sym_group, momentum,
                )
            })
        } else {
            self.pc_step_complex_inner(
                basis, coeffs, dt, tau_add, drop_tol, protected, sym_group, momentum,
            )
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn pc_step_complex_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<Complex<f64>>,
        dt: f64,
        tau_add: f64,
        drop_tol: f64,
        protected: &[Word],
        sym_group: Option<&TranslationGroup>,
        momentum: Option<&[i32]>,
    ) -> Result<(), Error> {
        // 0. Optional sector check on the input state.
        if let (Some(g), Some(k)) = (sym_group, momentum) {
            ppvm_pauli_sum::symmetry::check_momentum_sector(basis, coeffs, g, k, 1e-8)
                .map_err(|e| Error::Internal(format!(
                    "pc_step_complex: input is not in momentum sector {:?}: {}",
                    k, e
                )))?;
        }

        // 1. First-hop expansion (complex leakage).
        let leak = self.leakage_complex(basis, coeffs, protected)?;
        for (w, v) in leak {
            if v.norm() > tau_add {
                basis.push(w);
                coeffs.push(Complex::new(0.0, 0.0));
            }
        }
        // 2. Predictor.
        let coeffs_predict = self.expm_step_complex(basis, dt, coeffs, drop_tol);
        // 3. Second-hop expansion from the predicted state.
        let leak2 = self.leakage_complex(basis, &coeffs_predict, protected)?;
        drop(coeffs_predict);
        for (w, v) in leak2 {
            if v.norm() > tau_add {
                basis.push(w);
                coeffs.push(Complex::new(0.0, 0.0));
            }
        }
        // 4. Corrector.
        *coeffs = self.expm_step_complex(basis, dt, coeffs, drop_tol);
        // 5. Prune below drop_tol (complex magnitude).
        prune_basis_complex(basis, coeffs, drop_tol, protected);
        // NOTE: NO automatic symmetry-merge here. The orbit-rep
        // representation cannot be self-consistently evolved by the
        // standard L* action — applying L* to a single orbit-rep
        // would treat it as one Pauli string, not as the implicit
        // sum of orbit members with momentum-phases. Per-step
        // orbit-rep evolution requires phase-aware action with
        // complex CSR (to be added). The user is expected to call
        // `canonicalize_pauli_sum_complex` at observable-readout /
        // snapshot points instead.
        Ok(())
    }

    /// Complex-vector variant of [`Self::expm_step`]. The generator is the
    /// real in-basis matrix; only the vector is complex. Matrix-free.
    fn expm_step_complex(
        &self,
        basis: &[Word],
        dt: f64,
        b: &[Complex<f64>],
        drop_tol: f64,
    ) -> Vec<Complex<f64>> {
        // The generator is the REAL in-basis-restricted matrix; only the
        // vector `b` is complex. By linearity,
        // `exp(dt·M)·(re + i·im) = exp(dt·M)·re + i·exp(dt·M)·im`, so this is
        // two real matrix-free applies recombined. Fully matrix-free; no CSR.
        mf_expm::expm_apply_mf_cxvec(self, basis, dt, b, drop_tol)
    }

    /// One classical RK4 step on `O ← O + dt · L*(O)`, expanding the basis
    /// naturally as the action explores new strings. After the step, drops
    /// any string whose absolute coefficient is below `drop_tol` (protected
    /// words always kept). No predictor-corrector enrichment, no Krylov
    /// machinery, no CSR build — just four matrix-free action evaluations
    /// followed by a magnitude prune.
    ///
    /// Per-step local truncation error is `O(dt^5)` from the integrator.
    /// **Stability** requires `dt ≤ 2.78 / ‖L*‖` (RK4 absolute-stability
    /// boundary). For an n-qubit lattice Hamiltonian with bounded local
    /// terms this typically means `dt ≲ O(1) / n`; the dissipator further
    /// shrinks the bound at large `γ`. Violating the bound is **not
    /// signalled** — the trajectory will norm-conserve but individual Pauli
    /// coefficients diverge to oscillating ±large values that cancel; the
    /// observable looks fine, the basis still grows, but local quantities
    /// like MSD blow up. Always verify against a small in-band
    /// truncation case (e.g. against ED, or against [`Self::pc_step`]
    /// which is unconditionally stable) before trusting tight-`drop_tol`
    /// results at large `dt`.
    ///
    /// For stiff problems where the stability bound is restrictive, prefer
    /// [`Self::pc_step`], which integrates `exp(dt·L*)` via Krylov scaling-
    /// and-squaring and is unconditionally stable in `dt`.
    ///
    /// `num_threads`, when set, runs the entire step inside a freshly built
    /// rayon thread pool of that size.
    pub fn rk4_step(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        drop_tol: f64,
        protected: &[Word],
        num_threads: Option<usize>,
    ) -> Result<(), Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| self.rk4_step_inner(basis, coeffs, dt, drop_tol, protected))
        } else {
            self.rk4_step_inner(basis, coeffs, dt, drop_tol, protected)
        }
    }

    fn rk4_step_inner(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        drop_tol: f64,
        protected: &[Word],
    ) -> Result<(), Error> {
        // Helper: convert HashMap form back to (Vec<Word>, Vec<f64>) so we
        // can hand it to compute_action_sum, which expects slices.
        fn map_to_vecs(m: FxHashMap<Word, f64>) -> (Vec<Word>, Vec<f64>) {
            let n = m.len();
            let mut b = Vec::with_capacity(n);
            let mut c = Vec::with_capacity(n);
            for (w, v) in m {
                b.push(w);
                c.push(v);
            }
            (b, c)
        }

        // Helper: returns a fresh map representing O + α · K, where O is
        // (basis, coeffs) and K is a HashMap. We seed `out` with K (which
        // is typically smaller / has different support than O) so the
        // .entry() merges are fewer.
        fn scale_add(
            basis: &[Word],
            coeffs: &[f64],
            k: &FxHashMap<Word, f64>,
            alpha: f64,
        ) -> FxHashMap<Word, f64> {
            let mut out: FxHashMap<Word, f64> =
                FxHashMap::with_capacity_and_hasher(basis.len() + k.len(), FxBuildHasher::default());
            for (p, &c) in basis.iter().zip(coeffs.iter()) {
                out.insert(p.clone(), c);
            }
            for (w, v) in k {
                *out.entry(w.clone()).or_insert(0.0) += alpha * v;
            }
            out
        }

        // k1 = L*(O)
        let k1 = self.compute_action_sum(basis, coeffs)?;

        // k2 = L*(O + dt/2 · k1)
        let (b1, c1) = map_to_vecs(scale_add(basis, coeffs, &k1, dt / 2.0));
        let k2 = self.compute_action_sum(&b1, &c1)?;
        drop(b1);
        drop(c1);

        // k3 = L*(O + dt/2 · k2)
        let (b2, c2) = map_to_vecs(scale_add(basis, coeffs, &k2, dt / 2.0));
        let k3 = self.compute_action_sum(&b2, &c2)?;
        drop(b2);
        drop(c2);

        // k4 = L*(O + dt · k3)
        let (b3, c3) = map_to_vecs(scale_add(basis, coeffs, &k3, dt));
        let k4 = self.compute_action_sum(&b3, &c3)?;
        drop(b3);
        drop(c3);

        // O_new = O + dt/6 · (k1 + 2 k2 + 2 k3 + k4)
        let dt6 = dt / 6.0;
        let est_cap = basis.len() + k1.len() + k2.len() + k3.len() + k4.len();
        let mut out_map: FxHashMap<Word, f64> =
            FxHashMap::with_capacity_and_hasher(est_cap, FxBuildHasher::default());
        for (p, &c) in basis.iter().zip(coeffs.iter()) {
            out_map.insert(p.clone(), c);
        }
        for (w, v) in &k1 {
            *out_map.entry(w.clone()).or_insert(0.0) += dt6 * v;
        }
        for (w, v) in &k2 {
            *out_map.entry(w.clone()).or_insert(0.0) += dt6 * 2.0 * v;
        }
        for (w, v) in &k3 {
            *out_map.entry(w.clone()).or_insert(0.0) += dt6 * 2.0 * v;
        }
        for (w, v) in &k4 {
            *out_map.entry(w.clone()).or_insert(0.0) += dt6 * v;
        }
        drop(k1);
        drop(k2);
        drop(k3);
        drop(k4);

        // Prune below drop_tol; never drop a protected word, even if its
        // coefficient happens to land below threshold.
        let protected_hashes: FxHashSet<u64> = protected.iter().map(word_hash).collect();
        if drop_tol > 0.0 {
            out_map.retain(|w, &mut v| {
                v.abs() >= drop_tol || protected_hashes.contains(&word_hash(w))
            });
        }

        // Repack into the caller's Vec storage.
        basis.clear();
        coeffs.clear();
        basis.reserve(out_map.len());
        coeffs.reserve(out_map.len());
        for (w, v) in out_map {
            basis.push(w);
            coeffs.push(v);
        }
        Ok(())
    }

    /// Predictor-corrector adaptive step.
    ///
    /// Mutates `basis` (may grow) and `coeffs` in place to reflect the
    /// state at time `t + dt`. The step is:
    ///
    /// 1. **First-hop expansion**: compute leakage from the current state
    ///    and append the largest leakage Paulis to the basis (with starting
    ///    coefficient 0), up to `room = max_basis − basis.len()` new strings.
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
    /// `max_basis` is a hard rank cap on the live basis: enrichment adds at
    /// most `max_basis − basis.len()` of the largest leakage strings, the
    /// leakage map is capped to the same room, and the post-step basis is
    /// trimmed to the top-`max_basis` by `|coeff|` (protected words always
    /// kept). Pass a large value (e.g. `usize::MAX`) for the near-exact,
    /// uncapped case. `drop_tol` additionally prunes by magnitude.
    ///
    /// `num_threads`, when set, runs the entire step inside a freshly built
    /// rayon thread pool of that size — useful for benchmarking parallel
    /// scaling. When `None`, the global rayon pool is used.
    ///
    #[allow(clippy::too_many_arguments)]
    pub fn pc_step(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        max_basis: usize,
        drop_tol: f64,
        protected: &[Word],
        num_threads: Option<usize>,
        admit_basis: Option<usize>,
    ) -> Result<(), Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| {
                self.pc_step_inner(basis, coeffs, dt, max_basis, drop_tol, protected, admit_basis)
            })
        } else {
            self.pc_step_inner(basis, coeffs, dt, max_basis, drop_tol, protected, admit_basis)
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
        max_basis: usize,
        drop_tol: f64,
        protected: &[Word],
        num_threads: Option<usize>,
    ) -> Result<PcStepTimings, Error> {
        if let Some(n) = num_threads {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .map_err(|e| Error::Internal(format!("rayon pool build: {e}")))?;
            pool.install(|| {
                self.pc_step_inner_timed(
                    basis, coeffs, dt, max_basis, drop_tol, protected,
                )
            })
        } else {
            self.pc_step_inner_timed(
                basis, coeffs, dt, max_basis, drop_tol, protected,
            )
        }
    }

    fn pc_step_inner_timed(
        &self,
        basis: &mut Vec<Word>,
        coeffs: &mut Vec<f64>,
        dt: f64,
        max_basis: usize,
        drop_tol: f64,
        protected: &[Word],
    ) -> Result<PcStepTimings, Error> {
        let mut t = PcStepTimings::default();
        let tau_add = if dt > 0.0 { k_leakage() * drop_tol / dt } else { 0.0 };

        let t0 = Instant::now();
        let leak = self.leakage_with_prune(basis, coeffs, protected, max_basis, tau_add)?;
        t.leakage1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak, max_basis);
        t.expand1_us = t0.elapsed().as_micros() as u64;

        // Predictor: gencsr1 (or one_norm build for mf) + expm1.
        // For matrix-free, the "gencsr" timer captures the one-norm pass.
        // `coeffs` is read-only here, so we don't clone it — it serves as
        // the pre-step input for the corrector below as well.
        let t0 = Instant::now();
        let coeffs_predict = self.expm_step(basis, dt, coeffs, drop_tol);
        t.expm1_us = t0.elapsed().as_micros() as u64;

        let t0 = Instant::now();
        let leak2 = self.leakage_with_prune(basis, &coeffs_predict, protected, max_basis, tau_add)?;
        t.leakage2_us = t0.elapsed().as_micros() as u64;
        drop(coeffs_predict);

        let t0 = Instant::now();
        add_leakage_capped(basis, coeffs, leak2, max_basis);
        // (timed variant keeps the historical admission == max_basis path)
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
        max_basis: usize,
        drop_tol: f64,
        protected: &[Word],
        admit_basis: Option<usize>,
    ) -> Result<(), Error> {
        // Working-set (admission) bound: enrichment may grow the live basis
        // to `admit` >= `max_basis`; the final `cap_basis` then keeps the
        // top-`max_basis` strings by evolved |coeff| over the whole union
        // (retained + admitted), i.e. genuine rank displacement. With
        // `admit_basis = None` (the default) admission is bounded by
        // `max_basis` itself and `cap_basis` is a no-op — the historical
        // behaviour, where headroom must be created by `drop_tol` pruning.
        let admit = admit_basis.unwrap_or(max_basis).max(max_basis);
        let tau_add = if dt > 0.0 { k_leakage() * drop_tol / dt } else { 0.0 };
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

/// Complex-coefficient counterpart of [`prune_basis`]: drops entries
/// whose `|c|` (complex magnitude) is below `drop_tol`, never dropping
/// `protected` words.
fn prune_basis_complex(
    basis: &mut Vec<Word>,
    coeffs: &mut Vec<Complex<f64>>,
    drop_tol: f64,
    protected: &[Word],
) {
    if drop_tol <= 0.0 {
        return;
    }
    debug_assert_eq!(basis.len(), coeffs.len());
    let protected_set: FxHashSet<&Word> = protected.iter().collect();
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

    /// Per-step orbit-rep evolution gives the SAME final orbit-rep
    /// state as full-basis complex evolution followed by a single
    /// projection at the end. Validates that the phase-aware action
    /// + complex CSR machinery is consistent with the full-basis
    /// reference.
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
            spec.pc_step_complex(
                &mut bf, &mut cf, dt, 0.0, 0.0, &protected, None, None, None,
            )
            .unwrap();
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
                &spec, &mut br, &mut cr, dt, 10_000_000, 0.0, &protected, &group, &k,
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

    /// `pc_step_complex` correctly preserves the momentum sector of an
    /// initial k-eigenstate: evolve in full-basis form, project at the
    /// end, get a consistent orbit-rep representation.
    ///
    /// Also exercises the sector check on the input.
    #[test]
    fn pc_step_complex_preserves_momentum_sector() {
        use std::f64::consts::PI;
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

        // Momentum-k eigenstate.
        let mut basis: Vec<Word> = (0..n)
            .map(|j| {
                let mut s = vec!['I'; n];
                s[j] = 'Z';
                let (w, _) =
                    parse_pauli_string(&s.into_iter().collect::<String>(), n).unwrap();
                w
            })
            .collect();
        let mut coeffs: Vec<Complex<f64>> = (0..n as i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * (k_mode as f64) * (a as f64) / (n as f64)))
            .collect();

        let protected: Vec<Word> = Vec::new();

        // First step: pass sym_group + k to exercise the sector check.
        spec.pc_step_complex(
            &mut basis,
            &mut coeffs,
            dt,
            1.0,
            0.0,
            &protected,
            None,
            Some(&group),
            Some(&k),
        )
        .expect("sector check on k-eigenstate input should pass");
        // Subsequent steps without re-checking (would also pass).
        for _ in 1..n_steps {
            spec.pc_step_complex(
                &mut basis,
                &mut coeffs,
                dt,
                1.0,
                0.0,
                &protected,
                None,
                None,
                None,
            )
            .unwrap();
        }
        // After evolution, the state should STILL be in sector k.
        ppvm_pauli_sum::symmetry::check_momentum_sector(&basis, &coeffs, &group, &k, 1e-8)
            .expect("evolved state should remain in k-sector by translation-symmetry of L*");

        // Projecting to orbit-rep form should give a non-empty basis
        // with finite coefficients.
        let mut basis_proj = basis.clone();
        let mut coeffs_proj = coeffs.clone();
        ppvm_pauli_sum::symmetry::canonicalize_pauli_sum_complex(
            &mut basis_proj,
            &mut coeffs_proj,
            &group,
            &k,
        );
        assert!(basis_proj.len() > 0, "projected basis empty");
        assert!(basis_proj.len() <= basis.len(), "projection should shrink basis");
        assert!(
            coeffs_proj.iter().any(|c| c.norm() > 1e-10),
            "projected basis has only ~zero coefficients"
        );
    }

    /// Sector check rejects an input that mixes momentum components.
    #[test]
    fn pc_step_complex_rejects_wrong_sector() {
        let n = 4usize;
        let spec = LindbladSpec::new(n, &[("XIII".to_string(), 1.0)], &[]).unwrap();
        let group = ppvm_pauli_sum::symmetry::TranslationGroup::chain_1d(n);

        // A k=1 eigenstate (Z-sum with twist).
        let mut basis: Vec<Word> = (0..n)
            .map(|j| {
                let mut s = vec!['I'; n];
                s[j] = 'Z';
                let (w, _) =
                    parse_pauli_string(&s.into_iter().collect::<String>(), n).unwrap();
                w
            })
            .collect();
        use std::f64::consts::PI;
        let mut coeffs: Vec<Complex<f64>> = (0..n as i32)
            .map(|a| Complex::from_polar(1.0, -2.0 * PI * 1.0 * (a as f64) / (n as f64)))
            .collect();
        // Now try to evolve this k=1 state under the WRONG momentum (k=2).
        let wrong_k = vec![2_i32];
        let protected: Vec<Word> = Vec::new();
        let res = spec.pc_step_complex(
            &mut basis,
            &mut coeffs,
            0.01,
            1.0,
            0.0,
            &protected,
            None,
            Some(&group),
            Some(&wrong_k),
        );
        assert!(res.is_err(), "wrong-sector input should be rejected");
    }

    /// Complex pc_step at momentum k=0 (and no symmetry merging) must
    /// reproduce the real pc_step on the same trajectory exactly.
    #[test]
    fn pc_step_complex_matches_real_at_kzero_no_merge() {
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
            let tau_add = 0.0;
            spec.pc_step(
                &mut basis_r,
                &mut coeffs_r,
                dt,
                10_000_000,
                0.0,
                &protected,
                None,
                None,
            )
            .unwrap();
            spec.pc_step_complex(
                &mut basis_c,
                &mut coeffs_c,
                dt,
                tau_add,
                0.0,
                &protected,
                None,
                None, // no symmetry merging
                None,
            )
            .unwrap();
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
            // Set max_basis == current basis size so room = 0: NO leakage
            // enrichment, only the expm step (closest analog of standard
            // PP — the regime where merging commutes with evolution).
            // drop_tol = 0 → no truncation. This reproduces the old
            // tau_add = +∞ no-add behaviour under the rank-cap API.
            let max_basis_u = basis_u.len();
            spec.pc_step(
                &mut basis_u,
                &mut coeffs_u,
                dt,
                max_basis_u,
                0.0, // drop_tol = 0 → no truncation
                &protected,
                None,
                None,
            )
            .unwrap();

            let max_basis_m = basis_m.len();
            spec.pc_step(
                &mut basis_m,
                &mut coeffs_m,
                dt,
                max_basis_m,
                0.0,
                &protected,
                None,
                None,
            )
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
