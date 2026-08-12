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

use crate::gates::{Clifford, CliffordExtensions};

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

/// Opt-in marker selecting the blanket [`Clifford`] impl below. A type that
/// implements [`SymplecticColumns`] + [`PhaseTrack`] receives the shared,
/// single-audited blanket `Clifford` **only** if it also implements this marker.
///
/// # Why an opt-in marker (coherence)
///
/// On stable Rust two `impl Clifford for _` blocks may not overlap (E0119). A
/// type that wants a *fused* `Clifford` — one that reads each symplectic bit once
/// and folds in the phase in the same pass, rather than running the column and
/// phase primitives as separate steps — must provide its own `impl Clifford`.
/// That is illegal while an *unconditional* blanket `impl<T: SymplecticColumns +
/// PhaseTrack> Clifford for T` exists and the fused type satisfies those bounds.
/// Gating the blanket on this empty marker resolves the overlap without giving up
/// the shared audited copy: the standard word types (`PauliWord`,
/// `LossyPauliWord`, the future `Tableau`) opt in and get the blanket, while
/// `Phased<W>` stays out and supplies its own fused impl (see
/// `ppvm-phased-pauli-word-2`).
///
/// Design: §"Pauli algebra traits".
pub trait BlanketClifford {}

/// The derived [`Clifford`] behavior, blanket-implemented once for every type
/// that opts in via [`BlanketClifford`]. The *sequence* of primitives per gate is
/// identical across roles even though the phase primitive it calls is not — the
/// single audited copy of the symplectic sign logic that would otherwise be
/// duplicated and drift.
///
/// Design: §"Pauli algebra traits" (the `impl<T: SymplecticColumns + PhaseTrack +
/// BlanketClifford> Clifford for T` block).
impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> Clifford for T {
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

/// The derived [`CliffordExtensions`] behavior, blanket-implemented for the same
/// [`BlanketClifford`] opt-ins — the counterpart of the old crate's
/// `impl<T: PauliWordTrait> CliffordExtensions for T`.
///
/// # Why generator products rather than new primitives
///
/// The old blanket could write each extension gate as a raw bit rule because a
/// bare `PauliWordTrait` carries no phase, so the *only* content was the
/// `Sp(2n,2)` action. Here the blanket must also be correct for phase-carrying
/// opt-ins (`Tableau`), so each gate is expressed as a product of the audited
/// [`Clifford`] generators instead: the bit rule and the sign then follow from
/// generators whose signs are already machine-checked, rather than from six new
/// hand-written [`PhaseTrack`] deltas that would have to be re-proved (and
/// re-implemented by every word crate).
///
/// Calls compose in the *backward* Heisenberg convention this crate uses
/// (`P ↦ U†PU`; `lean/PPVM/Pauli/Conjugation.lean`, `conjSdag`): applying `A`
/// then `B` conjugates by the operator product `A·B`. With that, and using
/// `S³ = S†`, `SZ = S³`:
///
/// | gate | operator identity | call sequence |
/// |:---:|:---|:---|
/// | `s_dag` | `S† = S·Z` | `s`, `z` |
/// | `sqrt_x` | `√X ≃ H·S·H` | `h`, `s`, `h` |
/// | `sqrt_x_dag` | `√X† ≃ H·S†·H` | `h`, `s_dag`, `h` |
/// | `sqrt_y` | `√Y ≃ H·Z` | `h`, `z` |
/// | `sqrt_y_dag` | `√Y† ≃ Z·H` | `z`, `h` |
/// | `cy` | `CY = (I⊗S)·CNOT·(I⊗S†)` | `s(t)`, `cnot(c,t)`, `s_dag(t)` |
///
/// Each row reproduces the old crate's bit rule **and** the old phased word's
/// sign formula exactly — pinned by `tests/phase1_gate_surface.rs`
/// (`blanket_clifford_extensions_match_old_conjugation_table`,
/// `blanket_cy_matches_old_two_qubit_table`), which replays the gates on a
/// ℤ₄-phased stub and checks the full `s`/`s_dag`/`√X`/`√X†`/`√Y`/`√Y†`
/// conjugation table plus the 16-entry `CY` table of `ppvm-traits`; by the
/// `ppvm-conformance-2` differential suites against the old crate; and, for the
/// gate identities themselves, by the `ℤ[i]` matrix oracle
/// (`phased_pauli_word_lean.rs`).
///
/// The composition itself — the step the stub tests cannot check for a
/// phase-*carrying* opt-in — is machine-checked in
/// `lean/PPVM/Pauli/Conjugation.lean`, where each row above is *defined* as the
/// product of the audited generator homs (the crate's backward `s` is `conjSdag`
/// there) and is therefore a `MonoidHom` for free: `extSdag`/`extSqrtX`/
/// `extSqrtXdag`/`extSqrtY`/`extSqrtYdag` (+ `extSdagHom`… `MonoidHom.comp`s),
/// with the tables as `extSdag_eq_conjS`, `extSqrtX_X`/`_Y`/`_Z`, …, the
/// dagger-inverse pairs as `extSqrtXdag_extSqrtX`/`extSqrtYdag_extSqrtY`/
/// `extSdag_conjSdag`, and `extSqrtX_sq`/`extSqrtY_sq` for `√X² = X`, `√Y² = Y`.
/// `cy` is the same on `𝒫₂`: `conjCY` (= `conjST ∘ conjCNOT ∘ conjSdagT`, with
/// `conjCY_calls` collapsing the literal four-primitive sequence and
/// `conjCYHom` the hom), whose `conjCY_bits` + `conjCY_sign` *are* the old
/// 16-entry table. Corollary — every composite delta is still real
/// (`extSqrtX_isRealPhase`, …, `conjCY_isRealPhase`), so the `±1` drain in
/// `ppvm-pauli-sum-2` stays total on the extension gates too.
///
/// # Loss guard
///
/// A lossy word implements the guard inside its column primitives ("a gate
/// touching a lost qubit is a no-op"), so the single-qubit rows above inherit it
/// unchanged. `cy` is the one case worth stating: with a **lost control and a
/// present target**, the old whole-gate skip did nothing, while the decomposition
/// still runs `s(t)` and `s_dag(t)`. Those two share the `z ⊕= x` bit map and are
/// inverse conjugations, so they cancel exactly and the word is left untouched —
/// verified against the old reference over the full 25-word lossy alphabet in
/// `ppvm-conformance-2::lossy_pauli_word_diff`, and proven on every loss
/// configuration by `sActL_cnotActL_sActL_eq_cyActL`
/// (`lean/PPVM/Pauli/Symplectic.lean`): the guarded composite equals the old
/// crate's atomic whole-gate skip `cyActL`. The phase half of the same
/// cancellation, for a phase-carrying opt-in, is `conjS_conjSdag`/
/// `conjSdag_conjS` (`lean/PPVM/Pauli/Conjugation.lean`).
///
/// A concrete type that wants the old fused single-pass cost (the `Tableau`'s
/// per-gate bit-plane sweep) opts out of [`BlanketClifford`] and writes its own
/// `impl CliffordExtensions`, exactly as `Phased<W>` does for [`Clifford`].
impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> CliffordExtensions for T {
    #[inline]
    fn s_dag(&mut self, q: usize) {
        self.s(q);
        self.z(q);
    }

    #[inline]
    fn sqrt_x(&mut self, q: usize) {
        self.h(q);
        self.s(q);
        self.h(q);
    }

    #[inline]
    fn sqrt_x_dag(&mut self, q: usize) {
        self.h(q);
        self.s_dag(q);
        self.h(q);
    }

    #[inline]
    fn sqrt_y(&mut self, q: usize) {
        self.h(q);
        self.z(q);
    }

    #[inline]
    fn sqrt_y_dag(&mut self, q: usize) {
        self.z(q);
        self.h(q);
    }

    #[inline]
    fn cy(&mut self, control: usize, target: usize) {
        self.s(target);
        self.cnot(control, target);
        self.s_dag(target);
    }
}
