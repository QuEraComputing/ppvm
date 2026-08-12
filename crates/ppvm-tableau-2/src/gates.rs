// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The non-Clifford gate surface: the `T` gate, the single- and two-qubit Pauli
//! rotations, the in-plane rotation `R`, and `U3`.
//!
//! All of them branch the amplitude vector through
//! [`GeneralizedTableau::branch_with_coefficients`] (or, for `rotate_2`, the
//! bijective relabel), and therefore **auto-truncate**: the magnitude cutoff is
//! applied inline while the new vector is written. There is deliberately no
//! `truncate()` entry point on the tableau — unlike `PauliSum`, where truncation
//! is caller-driven.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits". The
//! branch each rotation stages — `c·P → cos·c·P + sin·c·(iGP)` — is a
//! norm-preserving, angle-additive 2-D rotation on the coefficient pair,
//! machine-checked in `lean/PPVM/Instantiations/Rotation.lean` (`rot_norm_sq`,
//! `rot_rot`), and its new key is genuinely new (`anticommute_new_key`).
//!
//! Ported from `ppvm-tableau/src/gates/{tgate,rot1,rot2,u3}.rs`.

use num::complex::{Complex, Complex64};
use ppvm_traits_2::{Pauli, RotXY, RotationOne, RotationTwo, TGate, U3Gate};

use crate::data::{Amplitudes, Bitstring, GeneralizedTableau};

/// `exp(iπ/8)·cos(π/8)`.
const COS_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.853_553_390_593_273_7,
    im: 0.353_553_390_593_273_8,
};
/// `-i·exp(iπ/8)·sin(π/8)`.
const ISIN_PI_OVER_8_TIMES_EXPIPI8: Complex64 = Complex {
    re: 0.146_446_609_406_726_24,
    im: -0.353_553_390_593_273_8,
};

impl<I: Bitstring, H> TGate for GeneralizedTableau<I, H> {
    fn t(&mut self, qubit: usize) {
        if self.is_lost[qubit] {
            return;
        }
        self.branch_with_coefficients(
            qubit,
            Pauli::Z,
            COS_PI_OVER_8_TIMES_EXPIPI8,
            ISIN_PI_OVER_8_TIMES_EXPIPI8,
        );
    }

    fn t_dag(&mut self, qubit: usize) {
        if self.is_lost[qubit] {
            return;
        }
        self.branch_with_coefficients(
            qubit,
            Pauli::Z,
            COS_PI_OVER_8_TIMES_EXPIPI8.conj(),
            ISIN_PI_OVER_8_TIMES_EXPIPI8.conj(),
        );
    }
}

impl<I: Bitstring, H> RotationOne<Complex64, f64> for GeneralizedTableau<I, H> {
    fn rotate_1(&mut self, axis: Pauli, qubit: usize, theta: f64) {
        if self.is_lost[qubit] {
            return;
        }
        let (sin, cos) = (theta * 0.5).sin_cos();
        let complex_cos = Complex64::new(cos, 0.0);
        let i_complex_sin = Complex64::new(0.0, -sin);
        self.branch_with_coefficients(qubit, axis, complex_cos, i_complex_sin);
    }
}

impl<I: Bitstring, H> RotXY<Complex64, f64> for GeneralizedTableau<I, H> {
    /// `R(axis_angle, θ) = RZ(axis_angle)·RX(θ)·RZ(−axis_angle)`.
    ///
    /// The tableau runs in the Schrödinger picture, so the sub-rotations are
    /// applied in forward order: `RZ(−axis_angle)`, then `RX(θ)`, then
    /// `RZ(axis_angle)`.
    fn r(&mut self, qubit: usize, axis_angle: f64, theta: f64) {
        self.rz(qubit, -axis_angle);
        self.rx(qubit, theta);
        self.rz(qubit, axis_angle);
    }
}

impl<I: Bitstring, H> U3Gate<Complex64, f64> for GeneralizedTableau<I, H> {
    /// `U3(θ, φ, λ) = RZ(φ)·RY(θ)·RZ(λ)`.
    fn u3(&mut self, qubit: usize, theta: f64, phi: f64, lambda: f64) {
        self.rz(qubit, lambda);
        self.ry(qubit, theta);
        self.rz(qubit, phi);
    }
}

/// Axis decoding for [`RotationTwo::rotate_2`]: `PAULIS[(axis_z << 1) | axis_x]`.
const PAULIS: [Pauli; 4] = [Pauli::I, Pauli::X, Pauli::Z, Pauli::Y];

impl<I: Bitstring, H> RotationTwo<Complex64, f64> for GeneralizedTableau<I, H> {
    /// `exp(−i·θ/2·P_a ⊗ P_b)`.
    ///
    /// # Loss
    ///
    /// If exactly one endpoint is lost the gate degrades to the **single-qubit**
    /// rotation on the surviving endpoint with the **same** angle (not halved,
    /// not otherwise adjusted). If both are lost the inner `rotate_1` no-ops via
    /// its own guard.
    ///
    /// # Ordering
    ///
    /// The branch stream is built by applying the Pauli on `b` **first**, then
    /// the Pauli on `a`, reproducing the old crate. That the composite of two
    /// single-site relabels *is* the frame-conjugated two-site Pauli — so this
    /// really is `exp(−i·θ/2·P_a ⊗ P_b)` — is `shiftOp_comp` in
    /// `lean/PPVM/Tableau/BranchPhase.lean`: the composite is again a phased XOR
    /// shift, with shift `L_a ⊕ L_b`, weight `w_a + w_b` and phase
    /// `pd_a + pd_b + 2⟨w_a, L_b⟩`. The `b`-before-`a` order turns out **not** to
    /// matter for the `ℤ/4` phase: `dot_crateWeight_order` isolates the whole
    /// order dependence in `⟨G_a, L_b⟩ + ⟨G_b, L_a⟩`, which
    /// `omega_eq_frame_coords` (`lean/PPVM/Tableau/Frame.lean`) identifies as
    /// `ω(P_a, P_b)`, zero for Paulis on distinct qubits
    /// (`omega_disjoint_support`, `rot2_order_irrelevant`;
    /// `crates/ppvm-conformance-2/tests/tableau_lean.rs` checks it on the real
    /// frames). The order is kept regardless — what it pins now is the float
    /// summation order. The merge keeps the absolute cutoff
    /// `|c| > |threshold|` (note: `abs`, not `norm_sqr` — the same *predicate*
    /// as the gates' `|c|² > threshold²` for a non-negative threshold, by
    /// `cutoff_abs_iff_sq` in `lean/PPVM/Algebra/Truncation.lean`, but a
    /// different float computation and a `hypot` per element; it is the
    /// measurement's rule that is genuinely different, being relative) and does
    /// **not** normalize.
    ///
    /// # Divergence from the old crate (adjudicated)
    ///
    /// Old merged into a `std::collections::HashMap` with the default
    /// `RandomState`, so the resulting amplitude **order** was randomized per
    /// process: two runs of the same seeded program produced different
    /// summation orders and therefore different float rounding in every later
    /// fold. There is no old order to preserve — it is not a function of the
    /// input — so this port merges deterministically (sort-merge, ascending by
    /// index, exactly as the `T` path) and pins the order. The support set and
    /// the per-index values are unchanged; only the previously-nondeterministic
    /// ordering becomes defined. `lean/PPVM/Algebra/GradedMap.lean` licenses the
    /// reordering over a *ring* (`accumulate_comm`/`accumulate_assoc`, and the
    /// batch-order-invariance `accumulateTerms_perm`), which `f64` is not —
    /// precisely the argument for pinning the order in the implementation rather
    /// than leaving it to the hasher.
    fn rotate_2(&mut self, axis_a: [u8; 2], axis_b: [u8; 2], a: usize, b: usize, theta: f64) {
        let [axis_a_x, axis_a_z] = axis_a;
        let [axis_b_x, axis_b_z] = axis_b;
        let pauli_a = PAULIS[(axis_a_z << 1 | axis_a_x) as usize];
        let pauli_b = PAULIS[(axis_b_z << 1 | axis_b_x) as usize];
        // NOTE: if both qubits are lost, the rot1 is a no-op.
        if self.is_lost[a] {
            return self.rotate_1(pauli_b, b, theta);
        } else if self.is_lost[b] {
            return self.rotate_1(pauli_a, a, theta);
        }

        let (sin, cos) = (theta * 0.5).sin_cos();
        let complex_cos = Complex64::new(cos, 0.0);
        let i_complex_sin = Complex64::new(0.0, -sin);

        let mut branch_coefficients = self.coefficients.clone();
        self.compute_coefficients_after_pauli_apply(&mut branch_coefficients, b, pauli_b);
        self.compute_coefficients_after_pauli_apply(&mut branch_coefficients, a, pauli_a);

        let old_coefficients = self.coefficients.take();

        // Deterministic sort-merge in place of old's `std::HashMap` (see the
        // divergence note above).
        let mut lhs: Vec<(I, Complex64)> = old_coefficients
            .into_iter()
            .map(|(coeff, idx)| (idx, complex_cos * coeff))
            .collect();
        let mut rhs: Vec<(I, Complex64)> = branch_coefficients
            .into_iter()
            .map(|(coeff, idx)| (idx, i_complex_sin * coeff))
            .collect();
        lhs.sort_unstable_by_key(|e| e.0);
        rhs.sort_unstable_by_key(|e| e.0);

        // `rotate_2`'s cutoff is absolute on the magnitude (not the squared
        // magnitude), and it does not normalize afterwards.
        let cutoff = Complex64::new(self.coefficient_threshold, 0.0).norm();
        self.coefficients = Amplitudes::new();
        self.coefficients.reserve(lhs.len() + rhs.len());
        let (mut i, mut j) = (0usize, 0usize);
        while i < lhs.len() && j < rhs.len() {
            match lhs[i].0.cmp(&rhs[j].0) {
                std::cmp::Ordering::Less => {
                    if lhs[i].1.norm() > cutoff {
                        self.coefficients.unsafe_insert(lhs[i].0, lhs[i].1);
                    }
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    if rhs[j].1.norm() > cutoff {
                        self.coefficients.unsafe_insert(rhs[j].0, rhs[j].1);
                    }
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    let sv = lhs[i].1 + rhs[j].1;
                    if sv.norm() > cutoff {
                        self.coefficients.unsafe_insert(lhs[i].0, sv);
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        while i < lhs.len() {
            if lhs[i].1.norm() > cutoff {
                self.coefficients.unsafe_insert(lhs[i].0, lhs[i].1);
            }
            i += 1;
        }
        while j < rhs.len() {
            if rhs[j].1.norm() > cutoff {
                self.coefficients.unsafe_insert(rhs[j].0, rhs[j].1);
            }
            j += 1;
        }
    }
}
