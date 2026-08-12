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

// ─── Column-wise row multiplication ───────────────────────────────────────

/// The running `ℤ/4` phase delta of a column sweep, one entry per generator.
///
/// `sign_parity` is a single XOR plane because `row_multiply`'s `2·sign_count`
/// term is `ℤ/2`-valued; `imag_count` needs the full `mod 4` and so rides the
/// `(imag_lo, imag_hi)` carry-save pair. The delta a generator ends up with is
/// `2·(sign_parity ⊕ imag_hi) + imag_lo`.
#[derive(Debug)]
pub(crate) struct PhaseAccumulator {
    /// Parity of `row_multiply`'s `sign` predicate.
    pub(crate) sign_parity: Vec<u64>,
    /// Low bit of `imag_count mod 4`.
    pub(crate) imag_lo: Vec<u64>,
    /// High bit of `imag_count mod 4`.
    pub(crate) imag_hi: Vec<u64>,
}

impl PhaseAccumulator {
    /// A zeroed accumulator for a frame of the given stride.
    #[inline]
    pub(crate) fn zeroed(stride: usize) -> Self {
        Self {
            sign_parity: vec![0; stride],
            imag_lo: vec![0; stride],
            imag_hi: vec![0; stride],
        }
    }

    /// The `(low, high)` phase-delta planes this accumulator represents.
    #[inline]
    pub(crate) fn delta(&self) -> (&[u64], Vec<u64>) {
        let high = self
            .sign_parity
            .iter()
            .zip(self.imag_hi.iter())
            .map(|(&s, &h)| s ^ h)
            .collect();
        (&self.imag_lo, high)
    }
}

/// The `ℤ/4` phase a set of generators picks up when one fixed Pauli is
/// multiplied into all of them, accumulated one **qubit column** at a time.
///
/// [`row_multiply`] folds one generator against another along the generator's
/// own bits, which is contiguous only in [`Orientation::RowMajor`](super::Orientation::RowMajor).
/// The measurement projection multiplies a *single* pivot into many generators
/// at once, and that shape transposes: hold the qubit fixed and the pivot's two
/// bits `(c, d)` at that qubit are **scalars**, so `row_multiply`'s `sign` and
/// `imag` predicates collapse to plane expressions over the selected generators.
/// Specialising the four `(c, d)` cases of
///
/// ```text
/// sign = (a&b&c&!d) | (a&!b&!c&d) | (!a&b&c&d)
/// imag = (a&!b&d) | (a&!c&d) | (!a&b&c) | (b&c&!d)
/// ```
///
/// gives, with `a` the X plane and `b` the Z plane over generators:
///
/// | pivot at this qubit | `sign` | `imag` |
/// |:--|:--|:--|
/// | `I` (`c=0, d=0`) | `0` | `0` |
/// | `X` (`c=1, d=0`) | `a & b` | `b` |
/// | `Z` (`c=0, d=1`) | `a & !b` | `a` |
/// | `XZ` (`c=1, d=1`) | `!a & b` | `a ^ b` |
///
/// so an identity column contributes nothing and is skipped outright.
///
/// # The two accumulators
///
/// `row_multiply` returns `(2·sign_count + imag_count) mod 4`. The `2·` makes
/// the sign term `ℤ/2`-valued, so only its **parity** survives — one XOR plane.
/// `imag_count` needs the full `mod 4`, so it rides a two-plane carry-save
/// counter. The delta is then `2·(sign_parity ⊕ imag_hi) + imag_lo`, i.e.
/// `imag_lo` is the low phase bit and `sign_parity ⊕ imag_hi` the high one.
///
/// `mask` selects the generators taking part; unselected ones contribute to
/// neither accumulator and are left untouched by [`xor_masked`].
#[inline]
pub(crate) fn accumulate_column_phase(
    x: &[u64],
    z: &[u64],
    mask: &[u64],
    pivot: (bool, bool),
    acc: &mut PhaseAccumulator,
) {
    let (pivot_x, pivot_z) = pivot;
    let PhaseAccumulator {
        sign_parity,
        imag_lo,
        imag_hi,
    } = acc;
    for i in 0..mask.len() {
        let m = mask[i];
        if m == 0 {
            continue;
        }
        let (a, b) = (x[i], z[i]);
        let (sign, imag) = match (pivot_x, pivot_z) {
            (false, false) => continue,
            (true, false) => (a & b, b),
            (false, true) => (a & !b, a),
            (true, true) => (!a & b, a ^ b),
        };
        sign_parity[i] ^= sign & m;
        let bit = imag & m;
        let carry = imag_lo[i] & bit;
        imag_lo[i] ^= bit;
        imag_hi[i] ^= carry;
    }
}

/// `dst ^= mask` where `flag` is set, elementwise. The bit-plane form of
/// "XOR a scalar bit into every selected generator".
#[inline]
pub(crate) fn xor_mask_if(dst: &mut [u64], mask: &[u64], flag: bool) {
    if !flag {
        return;
    }
    for i in 0..dst.len() {
        dst[i] ^= mask[i];
    }
}

/// Add a per-generator `ℤ/4` delta into a `(low, high)` phase plane pair,
/// restricted to `mask`.
///
/// The two-plane ripple `carry = lo & dlo; lo ^= dlo; hi ^= dhi ^ carry` is the
/// `ℤ/4` addition; `mask` is applied to the delta so unselected generators keep
/// their phase.
#[inline]
pub(crate) fn add_phase_planes(
    lo: &mut [u64],
    hi: &mut [u64],
    delta_lo: &[u64],
    delta_hi: &[u64],
    mask: &[u64],
) {
    for i in 0..lo.len() {
        let (dlo, dhi) = (delta_lo[i] & mask[i], delta_hi[i] & mask[i]);
        let carry = lo[i] & dlo;
        lo[i] ^= dlo;
        hi[i] ^= dhi ^ carry;
    }
}

/// Add the *scalar* `ℤ/4` value `delta` to every generator selected by `mask`.
#[inline]
pub(crate) fn add_scalar_phase(lo: &mut [u64], hi: &mut [u64], delta: u8, mask: &[u64]) {
    if delta == 0 {
        return;
    }
    let (dlo, dhi) = (delta & 1 == 1, delta & 2 == 2);
    for i in 0..lo.len() {
        let m = mask[i];
        let d_lo = if dlo { m } else { 0 };
        let d_hi = if dhi { m } else { 0 };
        let carry = lo[i] & d_lo;
        lo[i] ^= d_lo;
        hi[i] ^= d_hi ^ carry;
    }
}

/// The index of the lowest set bit of `words`, or `None`.
#[inline]
pub(crate) fn first_set(words: &[u64]) -> Option<usize> {
    words
        .iter()
        .position(|&w| w != 0)
        .map(|i| i * super::BITS_PER_WORD + words[i].trailing_zeros() as usize)
}
