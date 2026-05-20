// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::config::Config;

/// Z-magnetization-conserving two-qubit gates.
///
/// These gates generate unitaries that commute with the total
/// Z-magnetization `Σ_i Z_i`, so they preserve the U(1) symmetry sector
/// of any observable that already commutes with `Σ_i Z_i`. They are the
/// natural building blocks for XY, Heisenberg, and related spin-model
/// dynamics.
///
/// Like every other Pauli-propagation gate in `ppvm-runtime`, these
/// methods act in the **Heisenberg picture** — they conjugate the
/// observable held in a [`PauliSum`](crate::sum::PauliSum) by the gate's
/// unitary. Compose them with the rest of the circuit in reverse order,
/// the same way `rxx` / `ryy` / `rzz` are composed.
///
/// The default implementations express each gate as a composition of
/// existing `rxx` / `ryy` / `rzz` calls so that any backend implementing
/// [`RotationTwo`](crate::traits::RotationTwo) automatically supports
/// U(1)-conserving dynamics. Backends are free to override with a fused
/// implementation when they can avoid the intermediate branching.
///
/// # Conservation under truncation
///
/// `exchange`, `xyzz`, and `rzz` all commute with the total Z
/// magnetization `Σ_k Z_k`. As a consequence, an observable built from
/// `{I, Z}`-only Pauli strings (`Σ_i Z_i`, `Σ_{i<j} Z_iZ_j`, products
/// thereof) is conserved under propagation **up to per-gate
/// floating-point precision** — typically `< 1e-14` per gate,
/// accumulating linearly over a circuit. The off-diagonal cross terms
/// that propagation introduces from individual `Z_i` summands cancel by
/// coefficient when the contributions are summed, so they never grow
/// large enough to be retained by a sensible truncation cutoff.
///
/// Standard truncation thresholds — anything well below the conserved
/// coefficient magnitude (`~1`) and well above the per-gate ε floor —
/// preserve conservation in practice. Setting
/// [`CoefficientThreshold`](crate::strategy::CoefficientThreshold) close
/// to 1, or running so many gates that the accumulated drift becomes
/// comparable to the threshold, can break this guarantee. For more
/// general U(1)-conserving Hermitian observables, conjugate-paired
/// Pauli strings (`Y_iX_j` ↔ `X_iY_j`, etc.) carry coefficients of
/// equal magnitude, so symmetric coefficient truncation keeps the
/// observable inside the U(1) sub-algebra even when individual
/// coefficients are perturbed.
pub trait U1Conserving<T: Config> {
    /// Fused XY exchange rotation:
    /// `exp(-i θ/2 (X_a X_b + Y_a Y_b))`.
    ///
    /// `X_a X_b` and `Y_a Y_b` commute, so this is mathematically
    /// equivalent to `rxx(a, b, theta)` followed by `ryy(a, b, theta)`.
    /// Exposed as a single routine so backends can minimize intermediate
    /// branching and so users can write U(1)-symmetric dynamics without
    /// re-deriving the XX / YY decomposition every time.
    fn exchange(&mut self, a: usize, b: usize, theta: impl Into<T::Coeff>);

    /// Combined XY + ZZ Heisenberg-style interaction:
    ///
    /// ```text
    /// exp(-i θ_xy/2 (X_a X_b + Y_a Y_b)) · exp(-i θ_zz/2 Z_a Z_b)
    /// ```
    ///
    /// The two generators commute, so the factorization is exact. With
    /// `theta_zz = 0` this reduces to [`exchange`](Self::exchange); with
    /// `theta_xy = 0` it reduces to `rzz`. With both non-zero it
    /// implements one Trotter slice of an XXZ-style Hamiltonian.
    fn xyzz(
        &mut self,
        a: usize,
        b: usize,
        theta_xy: impl Into<T::Coeff>,
        theta_zz: impl Into<T::Coeff>,
    );
}
