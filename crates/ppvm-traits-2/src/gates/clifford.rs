// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::pauli::{BlanketClifford, PhaseTrack, SymplecticColumns};

/// The Clifford gate set, applied in the Heisenberg picture.
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

/// Additional Clifford gates beyond the minimal set: `S†`, `√X`, `√X†`, `√Y`,
/// `√Y†`, and `CY`.
pub trait CliffordExtensions: Clifford {
    /// Apply `S†` to one qubit.
    fn s_dag(&mut self, qubit: usize);
    /// Apply `√X` to one qubit.
    fn sqrt_x(&mut self, qubit: usize);
    /// Apply `(√X)†` to one qubit.
    fn sqrt_x_dag(&mut self, qubit: usize);
    /// Apply `√Y` to one qubit.
    fn sqrt_y(&mut self, qubit: usize);
    /// Apply `(√Y)†` to one qubit.
    fn sqrt_y_dag(&mut self, qubit: usize);
    /// Apply `CY` to one `(control, target)` pair.
    fn cy(&mut self, control: usize, target: usize);
    /// stim alias for [`cy`](CliffordExtensions::cy).
    fn zcy(&mut self, control: usize, target: usize) {
        self.cy(control, target)
    }
}

/// Batched Clifford gates: apply the same gate to many qubits in one call.
pub trait CliffordBatch: Clifford {
    /// Apply Pauli `X` to every qubit in `indices`.
    fn x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.x(q);
        }
    }
    /// Apply Pauli `Y` to every qubit in `indices`.
    fn y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.y(q);
        }
    }
    /// Apply Pauli `Z` to every qubit in `indices`.
    fn z_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.z(q);
        }
    }
    /// Apply Hadamard `H` to every qubit in `indices`.
    fn h_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.h(q);
        }
    }
    /// Apply the phase gate `S` to every qubit in `indices`.
    fn s_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s(q);
        }
    }
    /// Apply `CNOT` to every `(control, target)` pair.
    fn cnot_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cnot(c, t);
        }
    }
    /// Apply `CZ` to every `(control, target)` pair.
    fn cz_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cz(c, t);
        }
    }
}

pub trait CliffordExtensionsBatch: CliffordExtensions + CliffordBatch {
    /// Apply `S†` to every qubit in `indices`.
    fn s_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.s_dag(q);
        }
    }
    /// Apply `√X` to every qubit in `indices`.
    fn sqrt_x_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x(q);
        }
    }
    /// Apply `(√X)†` to every qubit in `indices`.
    fn sqrt_x_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_x_dag(q);
        }
    }
    /// Apply `√Y` to every qubit in `indices`.
    fn sqrt_y_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y(q);
        }
    }
    /// Apply `(√Y)†` to every qubit in `indices`.
    fn sqrt_y_dag_many(&mut self, indices: &[usize]) {
        for &q in indices {
            self.sqrt_y_dag(q);
        }
    }
    /// Apply `CY` to every `(control, target)` pair.
    fn cy_many(&mut self, pairs: &[(usize, usize)]) {
        for &(c, t) in pairs {
            self.cy(c, t);
        }
    }
}

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
