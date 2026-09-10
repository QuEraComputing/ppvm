// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// A single-qubit Pauli symbol — the site alphabet of an ordinary packed Pauli
/// word (`Word<Site = Pauli>`).
///
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Pauli {
    /// Identity `I`.
    I,
    /// Pauli `X`.
    X,
    /// Pauli `Y`.
    Y,
    /// Pauli `Z`.
    Z,
}
/// `Sp`-part: bit-plane column algebra, written **once** and shared by
/// `PhasedPauliWord` (1-bit columns) and `Tableau` (SIMD blocks over its `2n`
/// rows). Same meaning, different width. No phase — this is the role-independent
/// symplectic action.
///
/// Design: §"Pauli algebra traits". The bit rules realize the per-generator
/// `Sp(2n, 2)` isometries of `lean/PPVM/Pauli/Symplectic.lean`
/// (`hAct_isometry`/`sAct_isometry`/`cnotAct_isometry`/`czAct_isometry`).
pub trait SymplecticColumns {
    /// Number of qubits (columns) this operator spans.
    fn n_qubits(&self) -> usize;

    /// `H` on `q`: swap the X and Z columns.
    fn swap_xz(&mut self, q: usize);

    /// `S` on `q`: `z_q ⊕= x_q` (maps `X → Y`).
    ///
    /// (Completes the design's abbreviated `// ...`; see the module friction
    /// note.)
    fn xor_z_from_x(&mut self, q: usize);

    /// `CNOT` bit rule, part one: `x_tgt ⊕= x_ctrl`.
    fn xor_x_col(&mut self, ctrl: usize, tgt: usize);

    /// `CNOT` bit rule, part two: `z_ctrl ⊕= z_tgt`.
    fn xor_z_col(&mut self, tgt: usize, ctrl: usize);

    /// `CZ` bit rule on `(a, b)`: `z_a ⊕= x_b` and `z_b ⊕= x_a`.
    ///
    /// (Completes the design's abbreviated `// ...`; see the module friction
    /// note.)
    fn cz_bits(&mut self, a: usize, b: usize);
}

/// Extension-part: the phase algebra. `ℤ₄` for a phased word, `ℤ₂` + the
/// Aaronson–Gottesman `g`-rule for a tableau. One phase delta per gate; the
/// role-dependent half of conjugation, written **per type**.
pub trait PhaseTrack {
    /// `H` phase delta: flip the sign of a component with both `x` and `z` set
    /// (`Y → −Y`).
    fn flip_phase_where_xz(&mut self, q: usize);

    /// `S` phase delta on `q`.
    fn s_phase(&mut self, q: usize);

    /// `CNOT` phase delta on `(ctrl, tgt)`.
    fn cnot_phase(&mut self, ctrl: usize, tgt: usize);

    /// `CZ` phase delta on `(a, b)`.
    fn cz_phase(&mut self, a: usize, b: usize);

    /// `X` phase delta on `q` (pure sign; no bit change).
    fn x_phase(&mut self, q: usize);

    /// `Y` phase delta on `q` (pure sign; no bit change).
    fn y_phase(&mut self, q: usize);

    /// `Z` phase delta on `q` (pure sign; no bit change).
    fn z_phase(&mut self, q: usize);
}

/// Role-*exclusive* operations that interpret the rows as a symplectic basis
/// rather than as independent operators. A tableau-only trait a word never
/// implements. Holds the frame **primitives**, not `measure` itself — the two
/// measurement algorithms are built *on* these.
pub trait StabilizerFrame {
    /// Find a generator that anticommutes with the measured Pauli (the pivot).
    fn anticommuting_pivot(&self, qubit: usize) -> Option<usize>;

    /// Multiply generator `src` into `dst` (uses the Aaronson–Gottesman
    /// `g`-rule).
    fn row_multiply(&mut self, src: usize, dst: usize);

    /// Restore canonical form after elimination.
    fn canonicalize(&mut self);
}

pub trait BlanketClifford {}
