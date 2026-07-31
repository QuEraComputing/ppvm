// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Clifford conjugation for [`Phased`]: a **fused** `ℤ₄` implementation that reads
//! each inner symplectic bit once, computes the conjugation sign from those
//! reads, applies the bit update reusing them, and folds the sign into the stored
//! phase.
//!
//! # Why a fused impl rather than the blanket
//!
//! The standard word types get their [`Clifford`] from the single audited blanket
//! `impl<T: SymplecticColumns + PhaseTrack + BlanketClifford> Clifford for T` in
//! `ppvm-traits-2`, which runs the phase primitive and the column primitives as
//! **separate** steps. For a phased word that split means reading each inner X/Z
//! bit twice — once to compute the sign, once to apply the bit op — which
//! benchmarked ~1.6–1.8× slower than the old fused `PhasedPauliWord::cnot`.
//! `Phased<W>` therefore deliberately **does not** implement
//! [`BlanketClifford`](ppvm_traits_2::BlanketClifford); it supplies the fused
//! impl below, which reads each bit exactly once. The blanket stays the shared,
//! single-audited copy for the phaseless types (`PauliWord`, `LossyPauliWord`,
//! future `Tableau`).
//!
//! # Correctness
//!
//! The per-gate signs and bit updates are byte-for-byte the old fused kernel
//! `ppvm-pauli-word/src/phase/clifford.rs` (the `add_phase(_ << 1)` formulas plus
//! bit ops, including the lost-qubit guard): the sign is computed from the
//! *pre-mutation* bits, so reading them up front and writing afterward is
//! observationally identical. For a bare [`PauliWord`], [`PauliBits::is_lost`] is
//! a constant `false`, so the loss guard compiles away.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Pauli algebra traits:
//! symplectic structure and phase". Lean spec: the exact conjugation signs are
//! machine-checked in `lean/PPVM/Pauli/Conjugation.lean` — `conjH_sign`
//! (`H`: `+2` iff `x∧z`, i.e. `HYH = −Y`, `conjH_Y`), `conjSdag_sign`
//! (`S`: sign iff `x∧¬z`; this crate runs the backward `S†PS` direction — see the
//! `s` note on the `S`/`S†` convention), and the two-qubit public sign theorems
//! `conjCNOT_sign` / `conjCZ_sign`; the underlying bit maps are the `Sp(2n, 2)`
//! isometries of `lean/PPVM/Pauli/Symplectic.lean`.

use ppvm_traits_2::{Clifford, PauliBits, Phase};

use crate::data::Phased;

impl<W: PauliBits> Phased<W> {
    /// Fold a `−1` (`i²`) into the stored phase when `cond` holds. Every Clifford
    /// conjugation sign is a pure `±1` (a `<<1` phase delta in the old kernel), so
    /// each gate reduces its sign to this guarded `Neg1` compose.
    #[inline]
    fn flip_sign_if(&mut self, cond: bool) {
        if cond {
            self.phase *= Phase::Neg1;
        }
    }
}

/// The fused `ℤ₄` conjugation. Each gate reads the inner word's X/Z bits **once**,
/// computes the sign from those reads, applies the bit update reusing them, then
/// folds the sign into the stored phase. Ported from
/// `ppvm-pauli-word/src/phase/clifford.rs`, including the loss guard — for a bare
/// [`PauliWord`], [`PauliBits::is_lost`] is a constant `false`, so the guard is
/// compiled away.
impl<W: PauliBits> Clifford for Phased<W> {
    /// `X` conjugation (pure sign, no bit change): `XPX = (−1)^{z} P`.
    #[inline]
    fn x(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.z_bit(q));
    }

    /// `Y` conjugation (pure sign): `YPY = (−1)^{x ⊕ z} P`.
    #[inline]
    fn y(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.x_bit(q) ^ self.word.z_bit(q));
    }

    /// `Z` conjugation (pure sign): `ZPZ = (−1)^{x} P`.
    #[inline]
    fn z(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        self.flip_sign_if(self.word.x_bit(q));
    }

    /// `H` on `q`: swap the X/Z bits and pick up `Y → −Y`. Sign iff the component
    /// has both `x` and `z` set. Lean: `conjH_sign` (`+2` iff `x && z`), `conjH_Y`
    /// (`HYH = −Y`).
    #[inline]
    fn h(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_x_bit(q, z);
        self.word.set_z_bit(q, x);
        self.flip_sign_if(x && z);
    }

    /// `S` on `q`: `z ⊕= x`, sign iff `x ∧ ¬z`. In this crate's backward
    /// Heisenberg convention (`S†PS`; the old kernel's `test_s`) `S` on `X` yields
    /// `−Y`. Lean: `conjSdag_sign` pins this exact `+2`-iff-`x∧¬z` delta
    /// (`conjSdag_X`: `S†XS = −Y`), the backward counterpart of the forward
    /// `conjS_sign` (`x∧z`); `conjS_conjSdag` proves the two are inverse
    /// conjugations.
    #[inline]
    fn s(&mut self, q: usize) {
        if self.word.is_lost(q) {
            return;
        }
        let x = self.word.x_bit(q);
        let z = self.word.z_bit(q);
        self.word.set_z_bit(q, z ^ x);
        self.flip_sign_if(x && !z);
    }

    /// `CNOT` on `(ctrl, tgt)`: `x_tgt ⊕= x_ctrl` and `z_ctrl ⊕= z_tgt`, sign iff
    /// `x_ctrl ∧ z_tgt ∧ (x_tgt = z_ctrl)`. Lean: `conjCNOT_sign`.
    #[inline]
    fn cnot(&mut self, ctrl: usize, tgt: usize) {
        if self.word.is_lost(ctrl) || self.word.is_lost(tgt) {
            return;
        }
        let xc = self.word.x_bit(ctrl);
        let zc = self.word.z_bit(ctrl);
        let xt = self.word.x_bit(tgt);
        let zt = self.word.z_bit(tgt);
        self.word.set_x_bit(tgt, xt ^ xc);
        self.word.set_z_bit(ctrl, zc ^ zt);
        self.flip_sign_if(xc && zt && (xt == zc));
    }

    /// `CZ` on `(a, b)`: `z_a ⊕= x_b` and `z_b ⊕= x_a`, sign iff
    /// `x_a ∧ x_b ∧ (z_a ⊕ z_b)`. Lean: `conjCZ_sign`.
    #[inline]
    fn cz(&mut self, a: usize, b: usize) {
        if self.word.is_lost(a) || self.word.is_lost(b) {
            return;
        }
        let xa = self.word.x_bit(a);
        let za = self.word.z_bit(a);
        let xb = self.word.x_bit(b);
        let zb = self.word.z_bit(b);
        self.word.set_z_bit(a, za ^ xb);
        self.word.set_z_bit(b, zb ^ xa);
        self.flip_sign_if(xa && xb && (za ^ zb));
    }
}

#[cfg(test)]
mod tests {
    use crate::PhasedPauliWord;
    use ppvm_traits_2::Clifford;

    fn conj(input: &str, gate: impl Fn(&mut PhasedPauliWord)) -> String {
        let mut w: PhasedPauliWord = input.into();
        gate(&mut w);
        w.to_string()
    }

    #[test]
    fn single_qubit_gates_track_sign() {
        // (gate, input → target) ported from the old phase/clifford.rs tests.
        let x = [("+I", "+I"), ("+X", "+X"), ("+Y", "-Y"), ("+Z", "-Z")];
        let y = [("+I", "+I"), ("+X", "-X"), ("+Y", "+Y"), ("+Z", "-Z")];
        let z = [("+I", "+I"), ("+X", "-X"), ("+Y", "-Y"), ("+Z", "+Z")];
        let h = [("+I", "+I"), ("+X", "+Z"), ("+Y", "-Y"), ("+Z", "+X")];
        let s = [("+I", "+I"), ("+X", "-Y"), ("+Y", "+X"), ("+Z", "+Z")];
        for (i, t) in x {
            assert_eq!(conj(i, |w| w.x(0)), t, "X {i}");
        }
        for (i, t) in y {
            assert_eq!(conj(i, |w| w.y(0)), t, "Y {i}");
        }
        for (i, t) in z {
            assert_eq!(conj(i, |w| w.z(0)), t, "Z {i}");
        }
        for (i, t) in h {
            assert_eq!(conj(i, |w| w.h(0)), t, "H {i}");
        }
        for (i, t) in s {
            assert_eq!(conj(i, |w| w.s(0)), t, "S {i}");
        }
    }

    #[test]
    fn cnot_tracks_sign() {
        for (i, t) in [
            ("+II", "+II"),
            ("+IX", "+IX"),
            ("+IZ", "+ZZ"),
            ("+IY", "+ZY"),
            ("+XI", "+XX"),
            ("+XX", "+XI"),
            ("+XY", "+YZ"),
            ("+XZ", "-YY"),
            ("+ZI", "+ZI"),
            ("+ZX", "+ZX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+IZ"),
            ("+YI", "+YX"),
            ("+YX", "+YI"),
            ("+YY", "-XZ"),
            ("+YZ", "+XY"),
        ] {
            assert_eq!(conj(i, |w| w.cnot(0, 1)), t, "CNOT {i}");
        }
    }

    #[test]
    fn cz_tracks_sign() {
        for (i, t) in [
            ("+II", "+II"),
            ("+IX", "+ZX"),
            ("+IY", "+ZY"),
            ("+IZ", "+IZ"),
            ("+XI", "+XZ"),
            ("+XX", "+YY"),
            ("+XY", "-YX"),
            ("+XZ", "+XI"),
            ("+ZI", "+ZI"),
            ("+ZX", "+IX"),
            ("+ZY", "+IY"),
            ("+ZZ", "+ZZ"),
            ("+YI", "+YZ"),
            ("+YX", "-XY"),
            ("+YY", "+XX"),
            ("+YZ", "+YI"),
        ] {
            assert_eq!(conj(i, |w| w.cz(0, 1)), t, "CZ {i}");
        }
    }

    #[test]
    fn phase_carries_through_gate() {
        // A pre-existing +i factor is preserved under a signed conjugation:
        // H (−iY) = −i·(−Y) = +iY.
        assert_eq!(conj("-iY", |w| w.h(0)), "+iY");
    }
}
