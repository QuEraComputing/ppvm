// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The shared behavioral gate/noise traits: [`Clifford`], [`RotationOne`],
//! [`PauliError`], and [`Measure`].
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Behavioral traits". These
//! describe operations, not representation layout; Clifford gates need no
//! coefficient parameter, numeric operations take the coefficient type directly.

use crate::coefficient::{Angle, Coefficient};

/// The Clifford gate set, applied in the Heisenberg picture.
///
/// `Clifford` is a *derived* behavioral trait: it is **not** implemented by
/// hand on each type but blanket-implemented once over the Pauli algebra
/// primitives ([`crate::pauli::SymplecticColumns`] + [`crate::pauli::PhaseTrack`],
/// see `pauli.rs`), and separately on `Sum` in terms of its key's `Clifford`.
///
/// Design: §"Behavioral traits" and §"Pauli algebra traits". Each generator is
/// an `Sp(2n, 2)` isometry on the symplectic bits — machine-checked per
/// generator in `lean/PPVM/Pauli/Symplectic.lean`
/// (`hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/`czAct_isometry`) — with
/// the sign action of `lean/PPVM/Pauli/Conjugation.lean` (`conjH_Y`: `HYH = −Y`,
/// etc.).
pub trait Clifford {
    /// Apply Pauli `X` to one qubit.
    fn x(&mut self, qubit: usize);
    /// Apply Pauli `Y` to one qubit.
    fn y(&mut self, qubit: usize);
    /// Apply Pauli `Z` to one qubit.
    fn z(&mut self, qubit: usize);
    /// Apply Hadamard `H` to one qubit.
    fn h(&mut self, qubit: usize);
    /// Apply the phase gate `S` to one qubit.
    fn s(&mut self, qubit: usize);
    /// Apply `CNOT` to one `(control, target)` pair.
    fn cnot(&mut self, control: usize, target: usize);
    /// Apply `CZ` to one qubit pair.
    fn cz(&mut self, qubit0: usize, qubit1: usize);

    /// stim alias for [`cnot`](Clifford::cnot).
    fn cx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cnot`](Clifford::cnot).
    fn zcx(&mut self, control: usize, target: usize) {
        self.cnot(control, target)
    }
    /// stim alias for [`cz`](Clifford::cz).
    fn zcz(&mut self, qubit0: usize, qubit1: usize) {
        self.cz(qubit0, qubit1)
    }
}

/// Single-qubit rotations parameterized by an angle domain `A` that yields
/// amplitudes in coefficient domain `C`.
///
/// The angle defaults to the coefficient (`A = C`), recovering today's
/// `rx(theta: C)` while permitting a symbolic/parametric angle over an
/// `f64`-coefficient sum.
///
/// Design: §"Behavioral traits" (`RotationOne`). The branch each rotation stages
/// — `c·P → cos·c·P + sin·c·(iGP)` — produces exactly one genuinely-new term and
/// is a norm-preserving, angle-additive 2-D rotation on the coefficient pair,
/// machine-checked in `lean/PPVM/Instantiations/Rotation.lean`
/// (`anticommute_new_key`, `rot_norm_sq`, `rot_rot`).
pub trait RotationOne<C: Coefficient, A: Angle<C> = C> {
    /// Rotate about `X` on `qubit` by `theta`.
    fn rx(&mut self, qubit: usize, theta: A);
    /// Rotate about `Y` on `qubit` by `theta`.
    fn ry(&mut self, qubit: usize, theta: A);
    /// Rotate about `Z` on `qubit` by `theta`.
    fn rz(&mut self, qubit: usize, theta: A);
}

/// A unital single-qubit Pauli error channel `P ↦ λ_P·P`.
///
/// Design: §"Behavioral traits" (`PauliError`). Acting diagonally in the Pauli
/// basis, its transfer eigenvalue collapses (using `Σ_Q p_Q = 1`) to
/// `λ_P = 1 − 2·Σ_{Q anticommutes with P} p_Q`, machine-checked in
/// `lean/PPVM/Algebra/Noise.lean` (`pauli_channel_eigenvalue`, and
/// `pauli_channel_eigenvalue_omega` tying anticommutation to
/// `PPVM.Symplectic.omega`).
pub trait PauliError<C: Coefficient> {
    /// Apply a single-qubit Pauli channel with `X`, `Y`, `Z` probabilities.
    fn pauli_error(&mut self, qubit: usize, probabilities: [C; 3]);
}

/// Loss-aware projective computational-basis measurement.
///
/// `Some(false)` and `Some(true)` denote the `|0⟩`/`|1⟩` outcomes and `None`
/// denotes a lost qubit — the former `Measure -> bool` / `LossyMeasure ->
/// Option<bool>` split is removed. Sharing the result type does not share the
/// algorithm: `Tableau` uses the pure Clifford procedure, `GeneralizedTableau`
/// the coefficient-aware `O(n²)` decomposition.
///
/// Design: §"Behavioral traits" (`Measure`). The deterministic-vs-random
/// dichotomy the pivot search rests on is machine-checked in
/// `lean/PPVM/Tableau/Frame.lean` (`measurement_dichotomy`,
/// `measure_deterministic_iff_xfree`).
pub trait Measure {
    /// Measure `qubit`; `None` if the qubit has been lost.
    fn measure(&mut self, qubit: usize) -> Option<bool>;

    /// Measure each target in order, one result per target.
    fn measure_many(&mut self, targets: &[usize]) -> Vec<Option<bool>> {
        targets.iter().map(|&q| self.measure(q)).collect()
    }
}
