// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli-word type and the string / `u8`-label codecs.

use crate::Error;
use fxhash::FxBuildHasher;
use ppvm_pauli_word::word::PauliWord;
use ppvm_traits::PauliWordTrait;

/// Pauli-word storage chunk: `u64` on 64-bit targets, `u32` elsewhere
/// (`bitvec` implements `BitStore` for `u64` only on 64-bit targets, so
/// e.g. wasm32 builds use four 32-bit chunks instead of two 64-bit ones).
#[cfg(target_pointer_width = "64")]
pub(crate) type Chunk = u64;
#[cfg(not(target_pointer_width = "64"))]
pub(crate) type Chunk = u32;

/// Chunks per word; words pack up to 128 qubits on every target.
#[cfg(target_pointer_width = "64")]
pub(crate) const W_CHUNKS: usize = 2;
#[cfg(not(target_pointer_width = "64"))]
pub(crate) const W_CHUNKS: usize = 4;

/// Maximum number of qubits supported by [`Word`].
pub const MAX_QUBITS: usize = 128;

/// The Pauli-word storage type used throughout this crate.
///
/// `[Chunk; W_CHUNKS]` covers up to 128 qubits; the `FxBuildHasher`
/// matches the hash used by the `FxHashMap` keys we wrap with;
/// `REHASH=true` means `set()` keeps the cached hash in sync.
pub type Word = PauliWord<[Chunk; W_CHUNKS], FxBuildHasher, true>;

/// Build a [`Word`] from a length-`n_qubits` slice of Pauli labels
/// (`0=I, 1=X, 2=Z, 3=Y` — the [`ppvm_traits::char::Pauli`] discriminants).
/// Sets all bits and rehashes once.
pub fn word_from_codes(codes: &[u8]) -> Result<Word, Error> {
    let n_qubits = codes.len();
    if n_qubits > MAX_QUBITS {
        return Err(Error::TooManyQubits { got: n_qubits });
    }
    let mut w = Word::new(n_qubits);
    for (q, &b) in codes.iter().enumerate() {
        if b > 3 {
            return Err(Error::InvalidPauliCode { code: b });
        }
        if b & 1 != 0 {
            w.xbits.set(q, true);
        }
        if b & 2 != 0 {
            w.zbits.set(q, true);
        }
    }
    w.rehash();
    Ok(w)
}

/// Inverse of [`word_from_codes`]: write `n_qubits` Pauli labels into `out`.
pub fn codes_from_word(w: &Word, out: &mut [u8]) {
    assert_eq!(out.len(), w.n_qubits());
    for (q, slot) in out.iter_mut().enumerate() {
        let xb = w.xbits[q] as u8;
        let zb = w.zbits[q] as u8;
        *slot = xb | (zb << 1);
    }
}

/// Parse a `"IXYZ..."` string into a [`Word`] together with the list of
/// qubits where the Pauli is non-identity (the term's support).
pub fn parse_pauli_string(s: &str, n_qubits: usize) -> Result<(Word, Vec<u32>), Error> {
    if n_qubits > MAX_QUBITS {
        return Err(Error::TooManyQubits { got: n_qubits });
    }
    let chars: Vec<char> = s.chars().filter(|c| *c != '_').collect();
    if chars.len() != n_qubits {
        return Err(Error::WrongLength {
            expected: n_qubits,
            got: chars.len(),
        });
    }
    let mut w = Word::new(n_qubits);
    let mut support = Vec::new();
    for (q, c) in chars.into_iter().enumerate() {
        match c {
            'I' => {}
            'X' => {
                w.xbits.set(q, true);
                support.push(q as u32);
            }
            'Z' => {
                w.zbits.set(q, true);
                support.push(q as u32);
            }
            'Y' => {
                w.xbits.set(q, true);
                w.zbits.set(q, true);
                support.push(q as u32);
            }
            other => return Err(Error::InvalidPauliChar { c: other }),
        }
    }
    w.rehash();
    Ok((w, support))
}

/// Compute the support (non-identity qubits) of `w`.
pub(crate) fn word_support(w: &Word, out: &mut Vec<u32>) {
    out.clear();
    for q in 0..w.n_qubits() {
        if w.xbits[q] || w.zbits[q] {
            out.push(q as u32);
        }
    }
}

/// Compact 64-bit hash of a [`Word`], used as the key in cache-friendly
/// membership tables: an `FxHashMap<u64, ()>` over the basis has a working
/// set ~6× smaller than `FxHashMap<Word, ()>`. The hash mixes the word's
/// cached hash once through `FxHasher` and never touches the 32-byte
/// payload.
#[inline(always)]
pub(crate) fn word_hash(w: &Word) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = fxhash::FxHasher::default();
    w.hash(&mut h);
    h.finish()
}
