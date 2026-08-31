// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Precompiled Lindbladian: construction and the single-Pauli `L*` kernel.

use crate::Error;
use crate::algebra::{anti_commutes, comm_product, pauli_mul, phase_factor};
use crate::word::{MAX_QUBITS, Word, parse_pauli_string, word_support};
use fxhash::FxHashMap;
use num::Complex;

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
