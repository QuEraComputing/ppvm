// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Pauli-string expectation values for [`GeneralizedTableau`].
//!
//! Two entry points:
//!
//! - [`GeneralizedTableau::expectation`] — single-Pauli `⟨ψ|P|ψ⟩` for any
//!   [`Word<Site = Pauli>`](ppvm_traits_2::Word). Conjugates `P` through the
//!   frame and overlaps the resulting Pauli with the amplitude vector using the
//!   *same* formulas as the measurement code.
//! - [`GeneralizedTableau::z_expectation`] — the single-qubit `⟨Z⟩` fast path
//!   the asymmetric loss channel calls.
//!
//! Both are **non-mutating**: they never collapse the state and never
//! normalize, so a zero-probability projection (e.g. `⟨YY⟩` on a Bell state)
//! returns `0.0` rather than panicking.
//!
//! Ported from `ppvm-tableau/src/expectation.rs`.

use fxhash::FxHashMap as HashMap;
use num::complex::Complex64;
use ppvm_pauli_sum_2::PauliPattern;
use ppvm_traits_2::{Pauli, Word};

use crate::data::{Bitstring, GeneralizedTableau, RowStorage};

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
    /// `⟨ψ|word|ψ⟩` for the multi-qubit Pauli `word`.
    ///
    /// Conjugates `word` through the Clifford frame (giving a Pauli on the
    /// canonical basis: an X-mask, a Z-mask and an `i^φ` phase), then sums
    /// `⟨α|P_conj|β⟩ c_α* c_β` over the amplitude vector. Always real (a
    /// Hermitian operator on a normalized state).
    pub fn expectation<W: Word<Site = Pauli>>(&self, word: &W) -> f64 {
        let (phase, stab_anticomm, destab_anticomm) = self.compute_decomposition_word(word);
        if stab_anticomm == I::zero() {
            let entries: Vec<(Complex64, I)> = self.coefficients.iter().copied().collect();
            Self::compute_overlap_case_b(&entries, phase, destab_anticomm)
        } else {
            let coeff_map: HashMap<I, Complex64> =
                self.coefficients.iter().map(|&(c, i)| (i, c)).collect();
            let odd_phase_mask = self.odd_phase_destabilizer_mask();
            Self::compute_overlap_case_a(
                &coeff_map,
                phase,
                destab_anticomm,
                stab_anticomm,
                odd_phase_mask,
            )
        }
    }

    /// `⟨Z⟩` on `qubit`, computed non-destructively (the state is not
    /// collapsed).
    ///
    /// Reuses the measurement overlap machinery; cost scales with the number of
    /// coefficients (and `n²` for the decomposition).
    pub fn z_expectation(&self, qubit: usize) -> f64 {
        let (phase_decomp, stab_anticomm_bits, destab_anticomm_bits) =
            self.compute_decomposition(qubit, Pauli::Z);

        if stab_anticomm_bits == I::zero() {
            // Case b: `Z` is a stabilizer — self-pairing overlap.
            let entries: Vec<(Complex64, I)> = self.coefficients.iter().copied().collect();
            Self::compute_overlap_case_b(&entries, phase_decomp, destab_anticomm_bits)
        } else {
            // Case a: cross-index pairing — clone the support into a map
            // (read-only, so unlike `measure` this does not drain).
            let coeff_map: HashMap<I, Complex64> =
                self.coefficients.iter().map(|&(c, i)| (i, c)).collect();
            let odd_phase_mask = self.odd_phase_destabilizer_mask();
            Self::compute_overlap_case_a(
                &coeff_map,
                phase_decomp,
                destab_anticomm_bits,
                stab_anticomm_bits,
                odd_phase_mask,
            )
        }
    }

    /// Sum `⟨ψ|P|ψ⟩` over every Pauli word accepted by `pattern`.
    ///
    /// This is old `GeneralizedTableau::trace`: pattern enumeration is
    /// exponential by definition, while each leaf delegates to the audited
    /// single-word [`expectation`](Self::expectation) kernel. Unbounded star
    /// patterns panic, matching the old enumerator.
    pub fn trace(&self, pattern: &PauliPattern) -> f64 {
        pattern
            .enumerate_matches::<A>(self.n_qubits())
            .map(|word| self.expectation(&word))
            .sum()
    }
}
