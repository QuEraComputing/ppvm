// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli multiplication: XOR the X/Z planes and return the residual phase.
//! One pass counts sign and imaginary masks to compute `i^k`, with `k` modulo four.
//! The coefficient absorbs the phase; zero unused bits remain zero after XOR.
//! Phase identities are verified in `lean/PPVM/Pauli/{Phase,Word}.lean`.

use bitvec::view::BitView;
use num::PrimInt;
use ppvm_traits_2::{KeyProduct, Phase};
use std::hash::BuildHasher;

use crate::{HashFinalize, PauliStorage, PauliWord};

impl<A, H> KeyProduct for PauliWord<A, H>
where
    A: PauliStorage,
    <A as BitView>::Store: PrimInt,
    H: BuildHasher + Default + HashFinalize,
    // `Eq` and `Clone` add no further bounds on `H`.
{
    /// Returns the XORed word and its residual phase for the coefficient to absorb.
    /// Both inputs must have equal width; zero unused bits stay zero.
    fn key_mul(&self, other: &Self) -> (Self, Phase) {
        debug_assert_eq!(
            self.nqubits, other.nqubits,
            "twisted product requires equal-width words",
        );

        let lhs_x = self.xbits.data.as_raw_slice();
        let lhs_z = self.zbits.data.as_raw_slice();
        let rhs_x = other.xbits.data.as_raw_slice();
        let rhs_z = other.zbits.data.as_raw_slice();

        let mut product = Self::new(self.nqubits);
        let product_x = product.xbits.data.as_raw_mut_slice();
        let product_z = product.zbits.data.as_raw_mut_slice();
        let mut sign_count = 0u32;
        let mut imaginary_count = 0u32;

        // Each iteration multiplies a packed group of qubits in parallel.
        for chunk in 0..lhs_x.len() {
            let (left_x, left_z) = (lhs_x[chunk], lhs_z[chunk]);
            let (right_x, right_z) = (rhs_x[chunk], rhs_z[chunk]);

            // Each set bit contributes a factor of -1 (i²) or i, respectively.
            let sign_mask = (left_x & left_z & right_x & !right_z)
                | (left_x & !left_z & !right_x & right_z)
                | (!left_x & left_z & right_x & right_z);
            let imaginary_mask = (left_x & !left_z & right_z)
                | (left_x & !right_x & right_z)
                | (!left_x & left_z & right_x)
                | (left_z & right_x & !right_z);
            sign_count += sign_mask.count_ones();
            imaginary_count += imaginary_mask.count_ones();

            // XOR gives the product's Pauli letters; the masks account for phase.
            product_x[chunk] = left_x ^ right_x;
            product_z[chunk] = left_z ^ right_z;
        }

        product.refresh_hash();
        let phase_exponent = ((2 * sign_count + imaginary_count) % 4) as u8;
        (product, Phase::from_exponent(phase_exponent))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiplication_rejects_mismatched_widths() {
        for (lhs, rhs) in [("I", "IX"), ("IX", "I"), ("", "I"), ("I", "")] {
            let left = PauliWord::<u64>::from(lhs);
            let right = PauliWord::<u64>::from(rhs);
            assert!(
                std::panic::catch_unwind(|| left.key_mul(&right)).is_err(),
                "accepted mismatched widths: {lhs:?} * {rhs:?}"
            );
        }
    }

    /// A `Phase` rendered as a `+ / +i / - / -i` prefix, for readable asserts.
    fn phase_str(p: Phase) -> &'static str {
        match p {
            Phase::Pos1 => "+",
            Phase::PosI => "+i",
            Phase::Neg1 => "-",
            Phase::NegI => "-i",
        }
    }

    fn product_str(lhs: &str, rhs: &str) -> String {
        let x: PauliWord = lhs.into();
        let y: PauliWord = rhs.into();
        let (w, p) = x.key_mul(&y);
        format!("{}{}", phase_str(p), w)
    }

    #[test]
    fn single_qubit_products() {
        // Ported from `phase/mul.rs` tests; identity/phaseless cases added.
        for (lhs, rhs, ans) in [
            ("X", "X", "+I"),
            ("X", "Y", "+iZ"),
            ("X", "Z", "-iY"),
            ("Y", "Z", "+iX"),
            ("Z", "X", "+iY"),
            ("Y", "X", "-iZ"),
            ("I", "Y", "+Y"),
            ("Z", "Z", "+I"),
        ] {
            assert_eq!(product_str(lhs, rhs), ans, "{lhs}*{rhs}");
        }
    }

    #[test]
    fn multi_qubit_products() {
        for (lhs, rhs, ans) in [
            ("ZI", "ZI", "+II"),
            ("II", "ZI", "+ZI"),
            ("XX", "XX", "+II"),
        ] {
            assert_eq!(product_str(lhs, rhs), ans, "{lhs}*{rhs}");
        }
    }

    #[test]
    fn square_is_identity_up_to_phase() {
        // `phaseExpN_self`: P·P = +I (each Pauli squares to +I).
        for s in ["XYZI", "YYXZ", "ZZZZ", "IXYX"] {
            let (w, p) = PauliWord::<u64>::from(s).key_mul(&s.into());
            assert_eq!(w, PauliWord::<u64>::new(4), "{s}² word");
            assert_eq!(p, Phase::Pos1, "{s}² phase");
        }
    }

    #[test]
    fn commutation_sign_law() {
        // `phaseExpN_sub_comm`: P·Q = (−1)^{ω(P,Q)} Q·P. Compare the two orders'
        // phases: they are equal (commute) or differ by −1 (anticommute).
        for (l, r) in [("XY", "ZX"), ("XZ", "ZX"), ("YI", "IY"), ("XX", "ZZ")] {
            let (_, pq) = PauliWord::<u64>::from(l).key_mul(&r.into());
            let (_, qp) = PauliWord::<u64>::from(r).key_mul(&l.into());
            let ratio = pq.compose(qp.inverse());
            assert!(
                ratio == Phase::Pos1 || ratio == Phase::Neg1,
                "{l},{r}: phase ratio {ratio:?} is not ±1",
            );
        }
    }

    #[test]
    fn associativity() {
        // `tmul_assoc`: the twisted product is associative once the emitted
        // phase is folded onto a commutative coefficient (`Complex<f64>` here).
        let u: PauliWord = "XYZ".into();
        let v: PauliWord = "ZXY".into();
        let w: PauliWord = "YZX".into();

        let (uv, p_uv) = u.key_mul(&v);
        let (uv_w, p2) = uv.key_mul(&w);
        let left_word = uv_w;
        let left_phase = p_uv.compose(p2);

        let (vw, p_vw) = v.key_mul(&w);
        let (u_vw, p3) = u.key_mul(&vw);
        let right_word = u_vw;
        let right_phase = p_vw.compose(p3);

        assert_eq!(left_word, right_word);
        // Fold both onto a coefficient and compare.
        let one = num::Complex::new(1.0, 0.0);
        assert_eq!(left_phase.apply(&one), right_phase.apply(&one));
    }
}
