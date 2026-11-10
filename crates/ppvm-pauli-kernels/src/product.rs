// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The packed Pauli-word product and its phase.
//!
//! For words `P = ⊗ p_k`, `Q = ⊗ q_k` stored as X/Z bit planes,
//! `P · Q = iᵏ · ⊗ (p_k ⊕ q_k)`. Each site contributes `2·sign + imag` to `k`
//! (`site_product` in `lean/PPVM/Kernels/RowProduct.lean`), so `k` is
//! `2·popcount(sign) + popcount(imag)` mod 4 (`popcount_exponent`). Unused
//! (padding) bits must be zero; zero sites contribute nothing and stay zero.

use crate::Bits;

/// Rows up to this many words use an iterator loop, longer rows an indexed one.
/// Measured on arm64: the iterator loop is faster on short rows, the indexed
/// loop (which LLVM unrolls further) on long rows.
const SHORT_ROW: usize = 8;

/// The `(sign, imag)` masks for `p = (a, b)`, `q = (c, d)`.
#[inline(always)]
fn masks<W: Bits>(a: W, b: W, c: W, d: W) -> (W, W) {
    let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
    let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
    (sign, imag)
}

/// The exponent `k mod 4` from the mask counts. 4 divides 2³², so wrapping is exact.
#[inline(always)]
fn exponent(sign_count: u32, imag_count: u32) -> u8 {
    (sign_count.wrapping_mul(2).wrapping_add(imag_count) % 4) as u8
}

/// One storage word of the product `dst · src`.
///
/// Returns the product's bits `[dst_x ⊕ src_x, dst_z ⊕ src_z]` and this word's
/// contribution `2·popcount(sign) + popcount(imag)` to the phase exponent. A row's
/// exponent is the wrapping sum of its words' contributions, mod 4.
///
/// For containers that write their own loop around the product, e.g. to follow
/// it with a hash of the same words: a loop the compiler sees inline lets that
/// work reuse the words still in registers, which a call to [`row_product`] did
/// not (measured on `PauliWord`).
#[inline(always)]
pub fn word_product<W: Bits>(dst_x: W, dst_z: W, src_x: W, src_z: W) -> ([W; 2], u32) {
    let (sign, imag) = masks(dst_x, dst_z, src_x, src_z);
    let contribution = sign
        .popcount()
        .wrapping_mul(2)
        .wrapping_add(imag.popcount());
    ([dst_x ^ src_x, dst_z ^ src_z], contribution)
}

/// `dst ← dst · src`. Returns `k ∈ 0..4` with `dst · src = iᵏ · (dst ⊕ src)`.
///
/// All four slices must have the same length.
#[inline]
pub fn row_product<W: Bits>(dst_x: &mut [W], dst_z: &mut [W], src_x: &[W], src_z: &[W]) -> u8 {
    debug_assert!(
        dst_x.len() == dst_z.len() && src_x.len() == dst_x.len() && src_z.len() == dst_x.len()
    );
    let mut sign_count = 0u32;
    let mut imag_count = 0u32;
    if dst_x.len() <= SHORT_ROW {
        for (((dx, dz), &sx), &sz) in dst_x.iter_mut().zip(dst_z.iter_mut()).zip(src_x).zip(src_z) {
            let (sign, imag) = masks(*dx, *dz, sx, sz);
            sign_count = sign_count.wrapping_add(sign.popcount());
            imag_count = imag_count.wrapping_add(imag.popcount());
            *dx = *dx ^ sx;
            *dz = *dz ^ sz;
        }
    } else {
        for i in 0..dst_x.len() {
            let (a, b, c, d) = (dst_x[i], dst_z[i], src_x[i], src_z[i]);
            let (sign, imag) = masks(a, b, c, d);
            sign_count = sign_count.wrapping_add(sign.popcount());
            imag_count = imag_count.wrapping_add(imag.popcount());
            dst_x[i] = a ^ c;
            dst_z[i] = b ^ d;
        }
    }
    exponent(sign_count, imag_count)
}

/// The phase exponent of `a · b` without computing the product's bits.
///
/// All four slices must have the same length.
#[inline]
pub fn row_product_phase<W: Bits>(a_x: &[W], a_z: &[W], b_x: &[W], b_z: &[W]) -> u8 {
    debug_assert!(a_x.len() == a_z.len() && b_x.len() == a_x.len() && b_z.len() == a_x.len());
    let mut sign_count = 0u32;
    let mut imag_count = 0u32;
    if a_x.len() <= SHORT_ROW {
        for (((&ax, &az), &bx), &bz) in a_x.iter().zip(a_z).zip(b_x).zip(b_z) {
            let (sign, imag) = masks(ax, az, bx, bz);
            sign_count = sign_count.wrapping_add(sign.popcount());
            imag_count = imag_count.wrapping_add(imag.popcount());
        }
    } else {
        for i in 0..a_x.len() {
            let (sign, imag) = masks(a_x[i], a_z[i], b_x[i], b_z[i]);
            sign_count = sign_count.wrapping_add(sign.popcount());
            imag_count = imag_count.wrapping_add(imag.popcount());
        }
    }
    exponent(sign_count, imag_count)
}
