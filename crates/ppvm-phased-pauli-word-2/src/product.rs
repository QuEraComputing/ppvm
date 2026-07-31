// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The phased product for [`Phased`]`<W>` where `W: KeyProduct`.
//!
//! The product **reuses** the base word's [`KeyProduct::key_mul`] kernel; it does
//! not re-implement the `phaseExp` boolean formula. Writing `(φ, v)` and `(ψ, w)`
//! for the two phased words and `(v·w, iᵏ)` for `key_mul(v, w)` (result bits and
//! emitted residual phase),
//!
//! ```text
//! (φ, v) · (ψ, w) = (φ · ψ · iᵏ, v·w)
//! ```
//!
//! i.e. the two explicit prefactors and the residual `iᵏ` compose in the `ℤ₄`
//! phase group while the bits are exactly the base word's twisted product. This
//! is the full (untwisted) group product of the phased Pauli group `𝒫₁`: the
//! `KeyProduct` cocycle on the base is absorbed by carrying the phase explicitly.
//!
//! A lossy base word has **no** `KeyProduct` (loss breaks the twisted-product
//! group — a lost qubit has no well-defined product), so `Phased<LossyPauliWord>`
//! gets Clifford conjugation but no product; only `Phased<PauliWord>` and other
//! `KeyProduct` bases get this impl.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`" (`KeyProduct`; keys form a group only up to phase). Lean spec:
//! the phased product is the group operation of `𝒫₁` in
//! `lean/PPVM/Pauli/Phase.lean`, associative because the residual `phaseExp` is a
//! 2-cocycle (`phaseExp_cocycle`; n-qubit `lean/PPVM/Pauli/Word.lean`
//! `phaseExpN_cocycle`), with `P·P = +I` (`phaseExpN_self`).

use std::ops::{Mul, MulAssign};

use ppvm_traits_2::KeyProduct;

use crate::data::Phased;

impl<W: KeyProduct> Phased<W> {
    /// The phased product `(φ, v) · (ψ, w) = (φ · ψ · iᵏ, v·w)`, reusing the base
    /// word's [`KeyProduct::key_mul`] for both the result bits `v·w` and the
    /// residual phase `iᵏ`. The single product kernel; the [`Mul`] / [`MulAssign`]
    /// operators all route through it.
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let (word, emitted) = self.word.key_mul(&rhs.word);
        Self {
            word,
            phase: self.phase * rhs.phase * emitted,
        }
    }
}

impl<W: KeyProduct> Mul for Phased<W> {
    type Output = Phased<W>;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Phased::mul(&self, &rhs)
    }
}

impl<W: KeyProduct> Mul<&Phased<W>> for &Phased<W> {
    type Output = Phased<W>;

    #[inline]
    fn mul(self, rhs: &Phased<W>) -> Phased<W> {
        Phased::mul(self, rhs)
    }
}

impl<W: KeyProduct> MulAssign<&Phased<W>> for Phased<W> {
    #[inline]
    fn mul_assign(&mut self, rhs: &Phased<W>) {
        *self = Phased::mul(self, rhs);
    }
}

impl<W: KeyProduct> MulAssign for Phased<W> {
    #[inline]
    fn mul_assign(&mut self, rhs: Phased<W>) {
        *self = Phased::mul(self, &rhs);
    }
}

#[cfg(test)]
mod tests {
    use crate::PhasedPauliWord;

    fn product(lhs: &str, rhs: &str) -> String {
        let x: PhasedPauliWord = lhs.into();
        let y: PhasedPauliWord = rhs.into();
        (x * y).to_string()
    }

    #[test]
    fn single_qubit_products() {
        // Ported from the old phase/mul.rs tests.
        for (lhs, rhs, ans) in [("+X", "+X", "+I"), ("+X", "+Y", "+iZ"), ("+X", "+Z", "-iY")] {
            assert_eq!(product(lhs, rhs), ans, "{lhs}*{rhs}");
        }
    }

    #[test]
    fn multi_qubit_and_phase_accumulation() {
        for (lhs, rhs, ans) in [
            ("+ZI", "-ZI", "-II"),
            ("+II", "-ZI", "-ZI"),
            ("+XI", "+iXI", "+iII"),
            ("-XX", "-XX", "+II"),
        ] {
            assert_eq!(product(lhs, rhs), ans, "{lhs}*{rhs}");
        }
    }

    #[test]
    fn operator_forms_agree() {
        let a: PhasedPauliWord = "+iXY".into();
        let b: PhasedPauliWord = "-ZX".into();
        let by_ref = (&a * &b).to_string();
        let mut acc = a.clone();
        acc *= &b;
        assert_eq!(by_ref, acc.to_string());
        assert_eq!(by_ref, a.mul(&b).to_string());
        assert_eq!(by_ref, (a * b).to_string());
    }

    #[test]
    fn square_is_identity_up_to_phase() {
        // `phaseExpN_self`: the residual of g(w)·g(w) is +1 for every Pauli word
        // w, so a `+`-phased word squares to `+I…I` (the explicit prefactors are
        // both +1 here, isolating the residual).
        for s in ["+XYZI", "+YYXZ", "+ZZZZ", "+IXYX"] {
            let w: PhasedPauliWord = s.into();
            let sq = w.mul(&w);
            assert_eq!(sq.to_string(), "+IIII", "{s}²");
        }
    }
}
