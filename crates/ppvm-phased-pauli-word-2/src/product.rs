// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Multiplication for [`Phased<W>`] where `W: KeyProduct`.
//! Reuses [`KeyProduct::key_mul`] and combines its residual phase with both prefactors.
//! `Mul` and `MulAssign` delegate to the same product implementation.
//! Lossy words do not implement `KeyProduct`, so they have no phased product.

use std::ops::{Mul, MulAssign};

use ppvm_traits_2::KeyProduct;

use crate::data::Phased;

impl<W: KeyProduct> Phased<W> {
    /// Multiply the base words and combine both prefactors with the emitted phase.
    /// `(φ, v) · (ψ, w) = (φ · ψ · iᵏ, v·w)`, where `key_mul` emits `iᵏ`.
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
        // Every Pauli word with prefactor +1 squares to +I on all sites.
        for s in ["+XYZI", "+YYXZ", "+ZZZZ", "+IXYX"] {
            let w: PhasedPauliWord = s.into();
            let sq = w.mul(&w);
            assert_eq!(sq.to_string(), "+IIII", "{s}²");
        }
    }
}
