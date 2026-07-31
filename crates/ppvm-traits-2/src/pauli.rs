// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The Pauli algebra primitives that the derived [`Clifford`] trait is
//! blanket-implemented over, plus the role-exclusive [`StabilizerFrame`].
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits:
//! symplectic structure and phase". Conjugation by a Clifford factors into a
//! symplectic map on the bit planes ([`SymplecticColumns`], the `Sp` part,
//! written once) and a phase update ([`PhaseTrack`], the extension part, written
//! per type). The blanket `Clifford` impl is the single audited copy of the
//! symplectic sign logic.
//!
//! # Friction: primitive surface beyond the two named methods
//!
//! The design sketches [`SymplecticColumns`] and [`PhaseTrack`] with `H` and
//! `CNOT` spelled out and the rest abbreviated (`swap_xz`/`xor_x_col`/
//! `xor_z_col`, `flip_phase_where_xz`/`cnot_phase`, then `// ...one phase delta
//! per gate`). To make the blanket `Clifford` impl type-check over the full
//! generator set the trait exposes (`x`/`y`/`z`/`h`/`s`/`cnot`/`cz`), this module
//! completes that abbreviation with the standard per-generator primitives:
//! `xor_z_from_x` (the `S` bit op `z ⊕= x`) and `cz_bits` (the `CZ` bit op) on
//! [`SymplecticColumns`], and one phase-delta method per remaining generator on
//! [`PhaseTrack`]. The Pauli gates `X`/`Y`/`Z` touch **no** symplectic bit (they
//! are pure signs on anticommuting components), so they appear only as phase
//! deltas. Each primitive follows the exact `Sp`-part/extension-part split the
//! design prescribes; only the method *count* is completed, not the shape.

use crate::gates::Clifford;

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
///
/// Design: §"Pauli algebra traits". The signs realize the conjugation identities
/// of `lean/PPVM/Pauli/Conjugation.lean` (`conjH_sign`, `conjS_sign`, …).
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
///
/// Design: §"Pauli algebra traits" (`StabilizerFrame`). The `2n` generators form
/// a genuine symplectic basis preserved by every Clifford generator, and the
/// pivot search rests on the measurement dichotomy — machine-checked in
/// `lean/PPVM/Tableau/Frame.lean` (`IsSymplecticFrame`, `frame_linearIndependent`,
/// `isSymplecticFrame_identity`, `isSymplecticFrame_hAct`/`sAct`/`cnotAct`/`czAct`,
/// `measurement_dichotomy`).
pub trait StabilizerFrame {
    /// Find a generator that anticommutes with the measured Pauli (the pivot).
    fn anticommuting_pivot(&self, qubit: usize) -> Option<usize>;

    /// Multiply generator `src` into `dst` (uses the Aaronson–Gottesman
    /// `g`-rule).
    fn row_multiply(&mut self, src: usize, dst: usize);

    /// Restore canonical form after elimination.
    fn canonicalize(&mut self);
}

/// The derived [`Clifford`] behavior, blanket-implemented once. The *sequence*
/// of primitives per gate is identical across roles even though the phase
/// primitive it calls is not — the single audited copy of the symplectic sign
/// logic that would otherwise be duplicated and drift.
///
/// Design: §"Pauli algebra traits" (the `impl<T: SymplecticColumns + PhaseTrack>
/// Clifford for T` block).
impl<T: SymplecticColumns + PhaseTrack> Clifford for T {
    #[inline]
    fn x(&mut self, q: usize) {
        self.x_phase(q);
    }

    #[inline]
    fn y(&mut self, q: usize) {
        self.y_phase(q);
    }

    #[inline]
    fn z(&mut self, q: usize) {
        self.z_phase(q);
    }

    #[inline]
    fn h(&mut self, q: usize) {
        self.flip_phase_where_xz(q);
        self.swap_xz(q);
    }

    #[inline]
    fn s(&mut self, q: usize) {
        self.s_phase(q);
        self.xor_z_from_x(q);
    }

    #[inline]
    fn cnot(&mut self, c: usize, t: usize) {
        self.cnot_phase(c, t);
        self.xor_x_col(c, t);
        self.xor_z_col(t, c);
    }

    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        self.cz_phase(a, b);
        self.cz_bits(a, b);
    }
}
