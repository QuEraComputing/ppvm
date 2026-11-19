// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! One update rule per Clifford gate, in the Heisenberg picture `P ↦ G† P G`.
//!
//! A one-qubit rule takes the pre-gate bits `(x, z)` and returns
//! `([x', z'], flip)` with `G† P G = (−1)^flip · P'`. A two-qubit rule takes
//! `[x_a, z_a, x_b, z_b]` (for controlled gates `a` is the control) and returns
//! the four new bits and the flip.
//!
//! Each rule is named for its gate `G`, so `s` computes `S† P S`. For a forward
//! update `G P G†`, use the rule of `G†` (e.g. `s_dag` for an `S` gate); the
//! self-inverse gates X, Y, Z, H, CX, CZ, CY are the same in both directions.
//!
//! Bits encode the Hermitian Paulis I = (0,0), X = (1,0), Z = (0,1),
//! Y = (1,1); `S† X S = −Y`.
//!
//! Each rule is proved against its gate matrix in
//! `lean/PPVM/Kernels/Rules.lean`; the theorem is named on each function.

use crate::Bits;

/// `X`: sign flip where the Pauli has a Z. (`ruleX_correct`)
#[inline(always)]
pub fn x<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x, z], z)
}

/// `Y`: sign flip where the Pauli has exactly one of X, Z. (`ruleY_correct`)
#[inline(always)]
pub fn y<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x, z], x ^ z)
}

/// `Z`: sign flip where the Pauli has an X. (`ruleZ_correct`)
#[inline(always)]
pub fn z<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x, z], x)
}

/// `H`: swaps X and Z; `Y ↦ −Y`. (`ruleH_correct`)
#[inline(always)]
pub fn h<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([z, x], x & z)
}

/// `S† P S`: `X ↦ −Y`, `Y ↦ X`. (`ruleS_correct`)
#[inline(always)]
pub fn s<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x, z ^ x], x & !z)
}

/// `S P S†`: `X ↦ Y`, `Y ↦ −X`. (`ruleSdag_correct`)
#[inline(always)]
pub fn s_dag<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x, z ^ x], x & z)
}

/// `√X† P √X`: `Y ↦ −Z`, `Z ↦ Y`. (`ruleSqrtX_correct`)
#[inline(always)]
pub fn sqrt_x<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x ^ z, z], x & z)
}

/// `√X P √X†`: `Y ↦ Z`, `Z ↦ −Y`. (`ruleSqrtXdag_correct`)
#[inline(always)]
pub fn sqrt_x_dag<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([x ^ z, z], !x & z)
}

/// `√Y† P √Y`: `X ↦ Z`, `Z ↦ −X`. (`ruleSqrtY_correct`)
#[inline(always)]
pub fn sqrt_y<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([z, x], !x & z)
}

/// `√Y P √Y†`: `X ↦ −Z`, `Z ↦ X`. (`ruleSqrtYdag_correct`)
#[inline(always)]
pub fn sqrt_y_dag<B: Bits>(x: B, z: B) -> ([B; 2], B) {
    ([z, x], x & !z)
}

/// `CX` with control `c` and target `t`: `z_c ^= z_t`, `x_t ^= x_c`.
/// (`ruleCnot_correct`)
#[inline(always)]
pub fn cnot<B: Bits>([xc, zc, xt, zt]: [B; 4]) -> ([B; 4], B) {
    ([xc, zc ^ zt, xt ^ xc, zt], xc & zt & !(xt ^ zc))
}

/// `CZ` on `a`, `b`: `z_a ^= x_b`, `z_b ^= x_a`. (`ruleCz_correct`)
#[inline(always)]
pub fn cz<B: Bits>([xa, za, xb, zb]: [B; 4]) -> ([B; 4], B) {
    ([xa, za ^ xb, xb, zb ^ xa], xa & xb & (za ^ zb))
}

/// `CY` with control `c` and target `t`. (`ruleCy_correct`)
#[inline(always)]
pub fn cy<B: Bits>([xc, zc, xt, zt]: [B; 4]) -> ([B; 4], B) {
    (
        [xc, zc ^ xt ^ zt, xt ^ xc, zt ^ xc],
        xc & (xt ^ zt) & !(zc ^ zt),
    )
}
