// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The twisted (2-cocycle) key product [`KeyProduct`] for [`PauliWord`]:
//! `v·w = iᵏ (v⊕w)` with `k = phaseExp(v,w) ∈ ℤ/4`.
//!
//! Ported from `ppvm-pauli-word/src/phase/mul.rs` (the `MulAssign` kernel) to
//! keep the hot path at parity: one pass over the packed machine words computing
//! the result planes `x = a⊕c`, `z = b⊕d` and the packed boolean `sign`/`imag`
//! masks whose popcounts give `k = (2·sign_count + imag_count) mod 4`.
//!
//! Design: `traits-2-configuration-and-hashing.md` §"The map is a graded algebra
//! over `C[K]`" (`KeyProduct`; keys form a group only up to phase). Lean spec:
//! the packed boolean `2·sign + imag` formula equals the ℤ[i] matrix-product
//! exponent (`lean/PPVM/Pauli/Phase.lean` `phaseExp_eq_ref`, grounded in
//! `PauliMatrix.pauliMat_mul`), the exponent is a genuine 2-cocycle so the
//! product is associative (`phaseExp_cocycle`; n-qubit `Word.lean`
//! `phaseExpN_cocycle`), `P·P = +I` (`phaseExpN_self`), and
//! `P·Q = (−1)^{ω} Q·P` (`phaseExpN_sub_comm`).

use bitvec::view::BitView;
use num::PrimInt;
use ppvm_traits_2::{KeyProduct, Phase};
use std::hash::BuildHasher;

use crate::data::PauliWord;
use crate::storage::{HashFinalize, PauliStorage};

impl<A, H> KeyProduct for PauliWord<A, H>
where
    A: PauliStorage,
    <A as BitView>::Store: PrimInt,
    H: BuildHasher + Default + HashFinalize,
    // `KeyProduct: Eq + Clone` — both hold for `PauliWord<A, H>` with no `H`
    // bound (see `data.rs`).
{
    /// `self · other = iᵏ (self ⊕ other)`, returning the residual [`Phase`]
    /// `iᵏ` for the coefficient to absorb. The result inherits the
    /// canonical-unused-bits invariant: both operands have zero high bits, and
    /// `0 ⊕ 0 = 0`, so the XORed planes stay canonical.
    fn key_mul(&self, other: &Self) -> (Self, Phase) {
        debug_assert_eq!(
            self.nqubits, other.nqubits,
            "twisted product requires equal-width words",
        );

        let lhs_x = self.xbits.data.as_raw_slice();
        let lhs_z = self.zbits.data.as_raw_slice();
        let rhs_x = other.xbits.data.as_raw_slice();
        let rhs_z = other.zbits.data.as_raw_slice();

        let mut out = Self::new(self.nqubits);
        let mut sign_count = 0u32;
        let mut imag_count = 0u32;
        {
            let out_x = out.xbits.data.as_raw_mut_slice();
            let out_z = out.zbits.data.as_raw_mut_slice();
            for i in 0..lhs_x.len() {
                let a = lhs_x[i];
                let b = lhs_z[i];
                let c = rhs_x[i];
                let d = rhs_z[i];
                // `sign` contributes i² = −1, `imag` contributes i.
                let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
                let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
                sign_count += sign.count_ones();
                imag_count += imag.count_ones();
                out_x[i] = a ^ c;
                out_z[i] = b ^ d;
            }
        }
        out.invalidate_hash();
        let k = ((2 * sign_count + imag_count) % 4) as u8;
        (out, Phase::from_exponent(k))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
