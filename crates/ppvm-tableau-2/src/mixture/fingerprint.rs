// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::cell::RefCell;

use crate::GeneralizedTableau;

thread_local! {
    static WORD_BYTES: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[inline]
fn mask(index: usize, salt: u64) -> u64 {
    let mut z = (index as u64)
        .wrapping_mul(0x9e37_79b9_7f4a_7c15)
        .wrapping_add(salt);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

#[inline]
pub(crate) fn sign_mask(row: usize) -> u64 {
    mask(row, 0xa1a1_a1a1_a1a1_a1a1)
}

#[inline]
pub(crate) fn loss_mask(qubit: usize) -> u64 {
    mask(qubit, 0xc3c3_c3c3_c3c3_c3c3)
}

/// Digest of the frame's X/Z bits.
///
/// Reads the four quadrants as one contiguous byte range instead of
/// materializing `2n` rows — sound because the arena's padding is held at zero
/// (see [`crate::storage`]), so equal frames have equal bytes. The digest
/// *value* differs from the replaced row-by-row one, which is not observable:
/// a fingerprint only selects a collision bucket, and membership is decided by
/// `structurally_equal`.
pub(crate) fn word_fingerprint<I, H>(tab: &GeneralizedTableau<I, H>) -> u64 {
    let bits = tab.tableau.xz_bytes();
    #[cfg(not(target_arch = "wasm32"))]
    {
        let _ = &WORD_BYTES;
        gxhash::gxhash64(bits, 0)
    }
    #[cfg(target_arch = "wasm32")]
    {
        use std::hash::Hasher;
        let mut hasher = fxhash::FxHasher::default();
        hasher.write(bits);
        hasher.finish()
    }
}

pub(crate) fn phase_loss_hash<I, H>(tab: &GeneralizedTableau<I, H>) -> u64 {
    let mut hash = 0;
    for row in 0..2 * tab.tableau.n_qubits() {
        let phase = tab.tableau.row_phase(row);
        if phase & 1 != 0 {
            hash ^= mask(row, 0xb2b2_b2b2_b2b2_b2b2);
        }
        if phase & 2 != 0 {
            hash ^= sign_mask(row);
        }
    }
    for (qubit, &lost) in tab.is_lost.iter().enumerate() {
        if lost {
            hash ^= loss_mask(qubit);
        }
    }
    hash
}

pub(crate) fn fingerprint<I, H>(tab: &GeneralizedTableau<I, H>) -> u64 {
    word_fingerprint(tab) ^ phase_loss_hash(tab)
}
