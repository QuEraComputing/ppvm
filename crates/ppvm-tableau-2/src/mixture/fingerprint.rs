// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::cell::RefCell;

use crate::{GeneralizedTableau, RowStorage};

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

pub(crate) fn word_fingerprint<A: RowStorage, I, H>(tab: &GeneralizedTableau<A, I, H>) -> u64 {
    #[cfg(not(target_arch = "wasm32"))]
    {
        WORD_BYTES.with(|bytes| {
            let mut bytes = bytes.borrow_mut();
            bytes.clear();
            for (x, z, _) in tab.tableau.rows() {
                bytes.extend_from_slice(bytemuck::bytes_of(&x));
                bytes.extend_from_slice(bytemuck::bytes_of(&z));
            }
            gxhash::gxhash64(&bytes, 0)
        })
    }
    #[cfg(target_arch = "wasm32")]
    {
        use std::hash::Hasher;
        let mut hasher = fxhash::FxHasher::default();
        for (x, z, _) in tab.tableau.rows() {
            x.hash(&mut hasher);
            z.hash(&mut hasher);
        }
        hasher.finish()
    }
}

pub(crate) fn phase_loss_hash<A: RowStorage, I, H>(tab: &GeneralizedTableau<A, I, H>) -> u64 {
    let mut hash = 0;
    for (row, (_, _, phase)) in tab.tableau.rows().enumerate() {
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

pub(crate) fn fingerprint<A: RowStorage, I, H>(tab: &GeneralizedTableau<A, I, H>) -> u64 {
    word_fingerprint(tab) ^ phase_loss_hash(tab)
}
