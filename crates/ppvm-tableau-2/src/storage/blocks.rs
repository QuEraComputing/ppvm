// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The word kernels: one contiguous sweep per Clifford generator.
//!
//! Each function takes `stride`-word slices out of
//! [`TableauData`](super::TableauData) — the X plane and Z plane of the
//! qubit column(s) the gate touches, plus the sign plane of the same half —
//! and applies the gate to every generator at once. In the canonical
//! column-major orientation a single-qubit gate is therefore two such calls
//! (destabilizer half, stabilizer half) over `n.div_ceil(64)` words, where the
//! replaced row-major layout walked all `2n` generators to touch one word each.
//!
//! # Where the formulas come from
//!
//! Every body here is the *same predicate* as the fused row kernel it replaces
//! (`clifford.rs`), lifted from one bit to a whole word: the old
//! `pw.phase ^= (((x & z) & mask) != 0) as u8) << 1` on a single masked bit
//! becomes `ph ^= x & z` on the plane. Because the `+2` phase delta is
//! `ℤ/2`-valued the per-generator XOR *is* the `ℤ/4` update — the same identity
//! (`two_mul_natCast` in `lean/PPVM/Tableau/Batch.lean`) that licensed the old
//! batched sweeps. The low phase plane is untouched: no Clifford generator ever
//! changes a generator's `i`-parity, which is fixed by
//! `phase % 2 == popcount(x & z) % 2` for a Hermitian generator, and the bit
//! maps below preserve `popcount(x & z)` parity per generator.
//!
//! Each generator is an `Sp(2n, 2)` isometry with the sign action of
//! `lean/PPVM/Pauli/Conjugation.lean` (`conjH_sign`, `conjS_sign`, …); the
//! per-gate tables are unchanged from `clifford.rs`, only their access pattern is.
//!
//! # Padding
//!
//! Every phase predicate is an `AND` against an X or Z plane and every bit map
//! is elementwise, so words and bits that are zero on entry are zero on exit.
//! That is what keeps the "padding is zero" invariant the bulk equality and
//! hashing paths rely on.

/// `X`: bit-preserving, sign flips where the generator has a `Z`.
#[inline]
pub(crate) fn pauli_x(z: &[u64], ph: &mut [u64]) {
    for (p, &zw) in ph.iter_mut().zip(z) {
        *p ^= zw;
    }
}

/// `Y`: bit-preserving, sign flips where the generator has an `X` or a `Z` but
/// not both.
#[inline]
pub(crate) fn pauli_y(x: &[u64], z: &[u64], ph: &mut [u64]) {
    for ((p, &xw), &zw) in ph.iter_mut().zip(x).zip(z) {
        *p ^= xw ^ zw;
    }
}

/// `Z`: bit-preserving, sign flips where the generator has an `X`.
#[inline]
pub(crate) fn pauli_z(x: &[u64], ph: &mut [u64]) {
    for (p, &xw) in ph.iter_mut().zip(x) {
        *p ^= xw;
    }
}

/// `H`: swap the X and Z bits, sign flips where both were set.
#[inline]
pub(crate) fn h(x: &mut [u64], z: &mut [u64], ph: &mut [u64]) {
    for ((xw, zw), p) in x.iter_mut().zip(z.iter_mut()).zip(ph.iter_mut()) {
        let (a, b) = (*xw, *zw);
        *p ^= a & b;
        *xw = b;
        *zw = a;
    }
}

/// Backward `S`: `z ^= x`, sign flips where `x & z`.
#[inline]
pub(crate) fn s(x: &[u64], z: &mut [u64], ph: &mut [u64]) {
    for ((&a, zw), p) in x.iter().zip(z.iter_mut()).zip(ph.iter_mut()) {
        *p ^= a & *zw;
        *zw ^= a;
    }
}

/// `S†`: the `S` bit map, sign flips where `x & !z`.
#[inline]
pub(crate) fn s_dag(x: &[u64], z: &mut [u64], ph: &mut [u64]) {
    for ((&a, zw), p) in x.iter().zip(z.iter_mut()).zip(ph.iter_mut()) {
        *p ^= a & !*zw;
        *zw ^= a;
    }
}

/// `√X`: `x ^= z`, sign flips where `z & !x`.
#[inline]
pub(crate) fn sqrt_x(x: &mut [u64], z: &[u64], ph: &mut [u64]) {
    for ((xw, &b), p) in x.iter_mut().zip(z).zip(ph.iter_mut()) {
        *p ^= b & !*xw;
        *xw ^= b;
    }
}

/// `(√X)†`: the `√X` bit map, sign flips where `x & z`.
#[inline]
pub(crate) fn sqrt_x_dag(x: &mut [u64], z: &[u64], ph: &mut [u64]) {
    for ((xw, &b), p) in x.iter_mut().zip(z).zip(ph.iter_mut()) {
        *p ^= *xw & b;
        *xw ^= b;
    }
}

/// `√Y`: the `H` bit map, sign flips where `x & !z`.
#[inline]
pub(crate) fn sqrt_y(x: &mut [u64], z: &mut [u64], ph: &mut [u64]) {
    for ((xw, zw), p) in x.iter_mut().zip(z.iter_mut()).zip(ph.iter_mut()) {
        let (a, b) = (*xw, *zw);
        *p ^= a & !b;
        *xw = b;
        *zw = a;
    }
}

/// `(√Y)†`: the `H` bit map, sign flips where `z & !x`.
#[inline]
pub(crate) fn sqrt_y_dag(x: &mut [u64], z: &mut [u64], ph: &mut [u64]) {
    for ((xw, zw), p) in x.iter_mut().zip(z.iter_mut()).zip(ph.iter_mut()) {
        let (a, b) = (*xw, *zw);
        *p ^= b & !a;
        *xw = b;
        *zw = a;
    }
}

/// `CNOT`: `x_t ^= x_c`, `z_c ^= z_t`, sign flips where
/// `x_c & z_t & !(x_t ^ z_c)` — all read before either write.
#[inline]
pub(crate) fn cnot(xc: &[u64], zc: &mut [u64], xt: &mut [u64], zt: &[u64], ph: &mut [u64]) {
    for i in 0..ph.len() {
        let (a, b, c, d) = (xc[i], zc[i], xt[i], zt[i]);
        ph[i] ^= a & d & !(c ^ b);
        zc[i] = b ^ d;
        xt[i] = c ^ a;
    }
}

/// `CZ`: `z_a ^= x_b`, `z_b ^= x_a`, sign flips where `x_a & x_b & (z_a ^ z_b)`.
#[inline]
pub(crate) fn cz(xa: &[u64], za: &mut [u64], xb: &[u64], zb: &mut [u64], ph: &mut [u64]) {
    for i in 0..ph.len() {
        let (a, b, c, d) = (xa[i], za[i], xb[i], zb[i]);
        ph[i] ^= a & c & (b ^ d);
        za[i] = b ^ c;
        zb[i] = d ^ a;
    }
}

/// `CY`: `z_c ^= x_t ^ z_t`, `x_t ^= x_c`, `z_t ^= x_c`, sign flips where
/// `x_c & (x_t ^ z_t) & !(z_c ^ z_t)`.
#[inline]
pub(crate) fn cy(xc: &[u64], zc: &mut [u64], xt: &mut [u64], zt: &mut [u64], ph: &mut [u64]) {
    for i in 0..ph.len() {
        let (a, b, c, d) = (xc[i], zc[i], xt[i], zt[i]);
        let m = c ^ d;
        ph[i] ^= a & m & !(b ^ d);
        zc[i] = b ^ m;
        xt[i] = c ^ a;
        zt[i] = d ^ a;
    }
}

// ─── Row multiplication ───────────────────────────────────────────────────

/// The Aaronson–Gottesman `g`-rule product of two generators, accumulated over
/// whole words.
///
/// Multiplies `src` into `dst` in place and returns the `ℤ/4` phase the product
/// contributes, *excluding* `src`'s own phase (the caller adds that, matching
/// the old `Row::mul_assign`, which did `add_phase(g)` then `add_phase(rhs.phase)`).
///
/// Only meaningful in [`Orientation::RowMajor`](super::Orientation::RowMajor),
/// where a generator's bits are contiguous.
#[inline]
pub(crate) fn row_multiply(
    dst_x: &mut [u64],
    dst_z: &mut [u64],
    src_x: &[u64],
    src_z: &[u64],
) -> u8 {
    let mut sign_count = 0u32;
    let mut imag_count = 0u32;
    for i in 0..dst_x.len() {
        let (a, b, c, d) = (dst_x[i], dst_z[i], src_x[i], src_z[i]);
        let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
        let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
        sign_count += sign.count_ones();
        imag_count += imag.count_ones();
        dst_x[i] = a ^ c;
        dst_z[i] = b ^ d;
    }
    ((2 * sign_count + imag_count) % 4) as u8
}

/// The index of the lowest set bit of `words`, or `None`.
#[inline]
pub(crate) fn first_set(words: &[u64]) -> Option<usize> {
    words
        .iter()
        .position(|&w| w != 0)
        .map(|i| i * super::BITS_PER_WORD + words[i].trailing_zeros() as usize)
}
