// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Pauli-string expectation values for `GeneralizedTableau`.
//!
//! Two entry points:
//!
//! - [`GeneralizedTableau::expectation`] — single-Pauli `⟨ψ|P|ψ⟩` for a
//!   `PauliWord`. Conjugates `P` through the tableau and overlaps the
//!   resulting Pauli with the sparse coefficient vector using the same
//!   formulas as the measurement code.
//! - [`GeneralizedTableau::trace`] — `Σ_{P matches pattern} ⟨ψ|P|ψ⟩` for a
//!   `PauliPattern`. Enumerates the matching Paulis and sums their
//!   expectations.
//!
//! Decision 9 in the multi-backend plan calls these out as the natural
//! primitive for the tableau backend; semantics intentionally diverge from
//! the PauliSum trace.

use crate::data::GeneralizedTableau;
use crate::prelude::*;
use bitvec::view::BitView;
use fxhash::FxHashMap as HashMap;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};
use ppvm_pauli_word::pattern::PauliPattern;
use std::fmt::Debug;

impl<T, I, C> GeneralizedTableau<T, I, C>
where
    T: Config,
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: SparseVector<Complex<T::Coeff>, I> + Debug,
    T::Coeff:
        One + Zero + Clone + num::Num + ToPrimitive + Debug + std::ops::Mul<f64> + PartialOrd<f64>,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: TableauIndex + Debug,
{
    /// `⟨ψ|word|ψ⟩` for the multi-qubit Pauli `word`.
    ///
    /// Conjugates `word` through the Clifford tableau (giving a Pauli on the
    /// canonical basis: an X-mask, Z-mask, and `i^φ` phase), then sums
    /// `⟨α|P_conj|β⟩ c_α* c_β` over the sparse coefficient vector. Always
    /// returns a real number (Hermitian operator on a normalized state).
    pub fn expectation<W: PauliWordTrait>(&self, word: &W) -> f64 {
        let (phase, stab_anticomm, destab_anticomm) = self.compute_decomposition_word(word);
        if stab_anticomm == I::zero() {
            let entries: Vec<(Complex<T::Coeff>, I)> = self.coefficients.iter().copied().collect();
            Self::compute_overlap_case_b(&entries, phase, destab_anticomm)
        } else {
            let coeff_map: HashMap<I, Complex<T::Coeff>> =
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

    /// `Σ_{P matches pattern} ⟨ψ|P|ψ⟩`.
    ///
    /// Enumerates every `PauliWord` accepted by `pattern` via
    /// [`PauliPattern::enumerate_matches`] and sums their expectations.
    /// Star quantifiers (`X*`) panic — the pattern must be bounded; use
    /// counted repetition (`Z?{n}`) or positional anchors instead.
    pub fn trace(&self, pattern: &PauliPattern) -> f64 {
        let mut sum = 0.0f64;
        for word in pattern.enumerate_matches::<T::Storage>(self.n_qubits()) {
            sum += self.expectation(&word);
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_pauli_sum::config::fxhash::ByteF64;

    type TestTableau = GeneralizedTableau<ByteF64<1>>;

    fn word(s: &str) -> PauliWord<u64> {
        s.into()
    }

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "expected {expected}, got {actual} (|Δ| = {})",
            (actual - expected).abs()
        );
    }

    // ─── Single-qubit expectations ──────────────────────────────────────

    #[test]
    fn expectation_z_on_zero_state_is_one() {
        let tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        assert_close(tab.expectation(&word("Z")), 1.0, 1e-12);
    }

    #[test]
    fn expectation_x_on_zero_state_is_zero() {
        let tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        assert_close(tab.expectation(&word("X")), 0.0, 1e-12);
    }

    #[test]
    fn expectation_identity_on_zero_state_is_one() {
        let tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        assert_close(tab.expectation(&word("I")), 1.0, 1e-12);
    }

    #[test]
    fn expectation_x_on_plus_state_is_one() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.h(0);
        assert_close(tab.expectation(&word("X")), 1.0, 1e-12);
        assert_close(tab.expectation(&word("Z")), 0.0, 1e-12);
    }

    // ─── Bell state ⟨Φ+|·|Φ+⟩ ────────────────────────────────────────────

    fn bell() -> TestTableau {
        let mut tab: TestTableau = GeneralizedTableau::new(2, 1e-12);
        tab.h(0);
        tab.cnot(0, 1);
        tab
    }

    #[test]
    fn bell_state_pauli_expectations() {
        let tab = bell();
        assert_close(tab.expectation(&word("II")), 1.0, 1e-12);
        assert_close(tab.expectation(&word("ZZ")), 1.0, 1e-12);
        assert_close(tab.expectation(&word("XX")), 1.0, 1e-12);
        assert_close(tab.expectation(&word("YY")), -1.0, 1e-12);
        // Cross terms vanish for the Bell state.
        assert_close(tab.expectation(&word("IZ")), 0.0, 1e-12);
        assert_close(tab.expectation(&word("ZI")), 0.0, 1e-12);
        assert_close(tab.expectation(&word("XZ")), 0.0, 1e-12);
        assert_close(tab.expectation(&word("YX")), 0.0, 1e-12);
    }

    // ─── GHZ state ────────────────────────────────────────────────────

    #[test]
    fn ghz_state_expectations() {
        let mut tab: TestTableau = GeneralizedTableau::new(3, 1e-12);
        tab.h(0);
        tab.cnot(0, 1);
        tab.cnot(1, 2);
        // GHZ = (|000⟩ + |111⟩)/√2. For Z^z: eigenvalue is (-1)^{popcount(z)·x}
        // on |xxx⟩, so on the two basis states it agrees iff popcount(z) is
        // even; the diagonal expectation is then +1, otherwise 0.
        assert_close(tab.expectation(&word("III")), 1.0, 1e-12); // popcount 0 → +1
        assert_close(tab.expectation(&word("ZZZ")), 0.0, 1e-12); // popcount 3 → 0
        assert_close(tab.expectation(&word("ZIZ")), 1.0, 1e-12); // popcount 2 → +1
        assert_close(tab.expectation(&word("ZZI")), 1.0, 1e-12); // popcount 2 → +1
        assert_close(tab.expectation(&word("IZI")), 0.0, 1e-12); // popcount 1 → 0
        // XXX flips |000⟩ ↔ |111⟩, both in the GHZ superposition → +1.
        assert_close(tab.expectation(&word("XXX")), 1.0, 1e-12);
        // Y has off-diagonal action with imaginary phase; YYY contributes 0.
        assert_close(tab.expectation(&word("YYY")), 0.0, 1e-12);
    }

    // ─── Single-qubit rotation: |ψ⟩ = RY(θ)|0⟩ ────────────────────────

    #[test]
    fn ry_rotation_z_expectation_is_cos_theta() {
        // RY(θ)|0⟩ = cos(θ/2)|0⟩ + sin(θ/2)|1⟩.  ⟨ψ|Z|ψ⟩ = cos(θ).
        for theta in [0.0, 0.3, 1.0, std::f64::consts::PI / 2.0] {
            let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
            tab.ry(0, theta);
            assert_close(tab.expectation(&word("Z")), theta.cos(), 1e-12);
            assert_close(tab.expectation(&word("X")), theta.sin(), 1e-12);
        }
    }

    // ─── X / Y on a non-Clifford superposition ───────────────────────
    //
    // After H(0); T(0) the state is |ψ⟩ = (|0⟩ + e^{iπ/4}|1⟩)/√2:
    //   ⟨ψ|X|ψ⟩ = cos(π/4) = √2/2,
    //   ⟨ψ|Y|ψ⟩ = sin(π/4) = √2/2.
    // The T gate populates two branches in the sparse coefficient vector,
    // so the overlap is a cross-product between them — case_a, not case_b.
    // Y in particular drives `phase_decomp` to an odd value (Y on a |+⟩-style
    // frame conjugates to -Y), forcing the `phase == 1 | 3` arms of
    // `overlap_case_a`'s `match`. A sign bug there would flip ⟨Y⟩.

    #[test]
    fn t_plus_x_expectation_is_cos_pi_over_4() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.h(0);
        tab.t(0);
        assert_close(
            tab.expectation(&word("X")),
            std::f64::consts::FRAC_1_SQRT_2,
            1e-12,
        );
    }

    #[test]
    fn t_plus_y_expectation_is_sin_pi_over_4() {
        let mut tab: TestTableau = GeneralizedTableau::new(1, 1e-12);
        tab.h(0);
        tab.t(0);
        assert_close(
            tab.expectation(&word("Y")),
            std::f64::consts::FRAC_1_SQRT_2,
            1e-12,
        );
    }

    // ─── trace(pattern) ───────────────────────────────────────────────

    #[test]
    fn trace_of_z_or_identity_pattern_on_bell_is_two() {
        // For |Φ+⟩ = (|00⟩+|11⟩)/√2:
        //   Σ_{P ∈ {I,Z}^2} ⟨Φ+|P|Φ+⟩ = ⟨II⟩ + ⟨IZ⟩ + ⟨ZI⟩ + ⟨ZZ⟩
        //                              = 1     + 0     + 0     + 1 = 2.
        //   Equivalently 2^n |⟨0…0|Φ+⟩|² = 4 · 1/2 = 2.
        let tab = bell();
        let pat = PauliPattern::parse("Z?{2}").expect("parse Z?{2}");
        assert_close(tab.trace(&pat), 2.0, 1e-12);
    }

    #[test]
    fn trace_of_y_or_identity_pattern_on_bell_is_zero() {
        // For |Φ+⟩, Σ_{P ∈ {I,Y}^2} ⟨Φ+|P|Φ+⟩ = ⟨II⟩ + ⟨IY⟩ + ⟨YI⟩ + ⟨YY⟩
        //                                       = 1     + 0     + 0     + (-1) = 0.
        // Equivalently 2^n |⟨+i+i|Φ+⟩|² = 4 · 0 = 0 — the projection onto
        // the all-|+i⟩ state has zero amplitude. `trace` does the
        // enumeration sum here without ever mutating state or calling
        // `normalize`, so the zero-probability case doesn't panic.
        let tab = bell();
        let pat = PauliPattern::parse("Y?{2}").expect("parse Y?{2}");
        assert_close(tab.trace(&pat), 0.0, 1e-12);
    }

    #[test]
    fn trace_of_positional_pattern_on_bell_matches_single_pauli() {
        // `Z0Z1` matches exactly the word ZZ; trace should equal ⟨ZZ⟩ = 1.
        let tab = bell();
        let pat = PauliPattern::parse("Z0Z1").expect("parse Z0Z1");
        assert_close(tab.trace(&pat), 1.0, 1e-12);
    }

    // ─── Cross-backend: forward tableau shots vs backward PauliSum ────────
    //
    // Run a noisy Clifford circuit many times on a tableau and average the
    // per-shot ⟨ψ|ZZ|ψ⟩. Independently, Heisenberg-propagate ZZ backward
    // through the same circuit via `PauliSum` and read off ⟨0…0|U†ZZ U|0…0⟩
    // as the sum of coefficients over Z/I-only Paulis. The Monte-Carlo
    // average and the deterministic value must agree within sampling error.
    //
    // `PauliSum::g(i)` performs O → g† O g, so the backward sweep applies
    // gates in reverse time order. Depolarize is self-dual under Heisenberg.
    #[test]
    fn forward_shots_match_backward_pauli_sum_under_depolarizing_noise() {
        use ppvm_pauli_sum::config::indexmap::ByteFxHashF64;

        let p = 0.05_f64;
        let n_shots: u64 = 4000;
        let n_qubits = 2;

        let mut sum = 0.0_f64;
        for shot in 0..n_shots {
            let mut tab: TestTableau = GeneralizedTableau::new_with_seed(n_qubits, 1e-12, shot);
            tab.h(0);
            tab.depolarize1(0, p);
            tab.cnot(0, 1);
            tab.depolarize1(0, p);
            tab.depolarize1(1, p);
            sum += tab.expectation(&word("ZZ"));
        }
        let avg = sum / (n_shots as f64);

        let mut ps: PauliSum<ByteFxHashF64<1>> = PauliSum::builder().n_qubits(n_qubits).build();
        ps += ("ZZ", 1.0);
        ps.depolarize1(1, p);
        ps.depolarize1(0, p);
        ps.cnot(0, 1);
        ps.depolarize1(0, p);
        ps.h(0);
        let z_or_i = PauliPattern::parse("Z?{2}").expect("parse Z?{2}");
        let exact = ps.trace(&z_or_i);

        // Per-shot |⟨ZZ⟩| ≤ 1 ⇒ σ_mean ≤ 1/√N; 5σ keeps this robust to RNG draws.
        let tol = 5.0 / (n_shots as f64).sqrt();
        assert!(
            (avg - exact).abs() < tol,
            "tableau avg {avg} vs PauliSum exact {exact}, |Δ|={} (tol {tol})",
            (avg - exact).abs()
        );
    }
}
