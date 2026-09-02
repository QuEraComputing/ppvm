// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Hot-path Pauli product / commutator on raw word chunks.
//!
//! Phase encoding: Pauli product `P·Q = ω · R` where `ω = i^phase` and
//! `phase ∈ {0,1,2,3}` ↔ `ω ∈ {1, i, -1, -i}`. The per-byte XOR/AND
//! formulas are the same ones used by
//! [`ppvm_pauli_word::phase::PhasedPauliWord`]'s `MulAssign`. This module
//! keeps a copy that returns the unpacked `(word, phase)` pair without
//! constructing a phased wrapper.

use crate::word::{W_CHUNKS, Word};
use num::Complex;
use ppvm_traits::PauliWordTrait;

#[inline(always)]
pub(crate) fn phase_factor(phase: u8) -> Complex<f64> {
    match phase & 3 {
        0 => Complex::new(1.0, 0.0),
        1 => Complex::new(0.0, 1.0),
        2 => Complex::new(-1.0, 0.0),
        _ => Complex::new(0.0, -1.0),
    }
}

/// `true` if Pauli words `a` and `b` anti-commute.
///
/// Two Pauli strings anti-commute iff
/// `popcount(a.x & b.z) + popcount(a.z & b.x)` is odd.
#[inline(always)]
pub(crate) fn anti_commutes(a: &Word, b: &Word) -> bool {
    let mut bits: u32 = 0;
    for i in 0..W_CHUNKS {
        bits += (a.xbits.data[i] & b.zbits.data[i]).count_ones();
        bits += (a.zbits.data[i] & b.xbits.data[i]).count_ones();
    }
    bits & 1 == 1
}

/// Commutator product `h · p`: returns `(out, eps)` where `out = h ⊕ p` and
///
/// - `eps =  0` if `h` and `p` commute (caller should skip — `[h,p] = 0`),
/// - `eps = -2.0` if `h·p` has phase `+i` (so `i·[h,p] = -2·out`),
/// - `eps = +2.0` if `h·p` has phase `-i` (so `i·[h,p] = +2·out`).
#[inline(always)]
pub(crate) fn comm_product(h: &Word, p: &Word) -> (Word, f64) {
    let (out, phase) = pauli_mul(h, p);
    let eps = match phase {
        1 => -2.0,
        3 => 2.0,
        _ => 0.0,
    };
    (out, eps)
}

/// Full Pauli product `p · q`: returns `(out, phase)` where the product
/// is `ω · out` with `ω = i^phase`.
#[inline(always)]
pub(crate) fn pauli_mul(p: &Word, q: &Word) -> (Word, u8) {
    let mut out = Word::new(p.n_qubits());
    let mut sign_count: u32 = 0;
    let mut imag_count: u32 = 0;
    for i in 0..W_CHUNKS {
        let a = p.xbits.data[i];
        let b = p.zbits.data[i];
        let c = q.xbits.data[i];
        let d = q.zbits.data[i];
        let sign = (a & b & c & !d) | (a & !b & !c & d) | (!a & b & c & d);
        let imag = (a & !b & d) | (a & !c & d) | (!a & b & c) | (b & c & !d);
        sign_count += sign.count_ones();
        imag_count += imag.count_ones();
        out.xbits.data[i] = a ^ c;
        out.zbits.data[i] = b ^ d;
    }
    out.rehash();
    (out, ((2 * sign_count + imag_count) & 3) as u8)
}
