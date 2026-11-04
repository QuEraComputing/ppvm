// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Packed Pauli-word type and the string / `u8`-label codecs.

use crate::Error;
use bitvec::array::BitArray;
use fxhash::FxBuildHasher;
use ppvm_traits::traits::HashFinalize;
use std::hash::{Hash, Hasher};

/// Pauli-word storage chunk: `u64` on 64-bit targets, `u32` elsewhere
/// (`bitvec` implements `BitStore` for `u64` only on 64-bit targets, so
/// e.g. wasm32 builds use four 32-bit chunks instead of two 64-bit ones).
#[cfg(target_pointer_width = "64")]
pub(crate) type Chunk = u64;
#[cfg(not(target_pointer_width = "64"))]
pub(crate) type Chunk = u32;

/// Bits per storage chunk.
pub const CHUNK_BITS: usize = Chunk::BITS as usize;

/// Chunks for a 128-qubit word: the default width, and the only one
/// instantiated for registers of at most 128 qubits.
pub const W_CHUNKS: usize = 128 / CHUNK_BITS;

/// Chunk counts of the word widths the Python layer instantiates, narrowest
/// first: 128, 256 and 512 qubits. [`chunks_for`] picks among them.
pub const WIDTHS: [usize; 3] = [128 / CHUNK_BITS, 256 / CHUNK_BITS, 512 / CHUNK_BITS];

/// Largest register any instantiated width supports.
pub const MAX_SUPPORTED_QUBITS: usize = 512;

/// Maximum number of qubits of the default-width [`Word`] (128).
pub const MAX_QUBITS: usize = max_qubits::<W_CHUNKS>();

/// Capacity of a `C`-chunk [`Word`], in qubits.
pub const fn max_qubits<const C: usize>() -> usize {
    C * CHUNK_BITS
}

/// The narrowest entry of [`WIDTHS`] that holds `n_qubits`, or `None` above
/// [`MAX_SUPPORTED_QUBITS`].
pub const fn chunks_for(n_qubits: usize) -> Option<usize> {
    let mut i = 0;
    while i < WIDTHS.len() {
        if n_qubits <= WIDTHS[i] * CHUNK_BITS {
            return Some(WIDTHS[i]);
        }
        i += 1;
    }
    None
}

/// The Pauli-word storage type used throughout this crate, `C` chunks wide:
/// one X and one Z bit plane, 32 bytes at the default 128-qubit width.
///
/// The crate is const-generic in the chunk count so a register of any size
/// up to [`MAX_SUPPORTED_QUBITS`] gets a fixed-size word. The register size
/// belongs to the model and is not stored per word; bits at or above it are
/// always zero. The hash is computed from the planes on demand with the same
/// digest `ppvm_pauli_word::PauliWord` caches, so hash tables keyed by words
/// order them exactly as they would order that type. Ordering compares the X
/// plane, then the Z plane.
#[derive(Clone, Copy, Default)]
pub struct Word<const C: usize = W_CHUNKS> {
    /// X-bit plane, one bit per qubit.
    pub xbits: BitArray<[Chunk; C]>,
    /// Z-bit plane, one bit per qubit.
    pub zbits: BitArray<[Chunk; C]>,
}

impl<const C: usize> Word<C> {
    /// The identity word for a register of `n_qubits` (at most the width's capacity).
    pub fn new(n_qubits: usize) -> Self {
        debug_assert!(n_qubits <= max_qubits::<C>());
        Self::default()
    }

    #[inline(always)]
    fn digest(&self) -> u64 {
        FxBuildHasher::word_hash(
            &[
                bytemuck::bytes_of(&self.xbits.data),
                bytemuck::bytes_of(&self.zbits.data),
            ],
            std::mem::size_of::<[Chunk; C]>(),
        )
    }
}

/// The general `ppvm_pauli_word` word of the same width, which also stores the
/// register size and a cached hash. `ppvm_pauli_sum` symmetry routines take this type.
pub type GeneralWord<const C: usize = W_CHUNKS> =
    ppvm_pauli_word::word::PauliWord<[Chunk; C], FxBuildHasher, true>;

impl<const C: usize> Word<C> {
    /// This word as the general `PauliWord` used by `ppvm_pauli_sum` symmetry groups.
    pub fn to_pauli_word(self, n_qubits: usize) -> GeneralWord<C> {
        use ppvm_traits::PauliWordTrait;
        let mut word = ppvm_pauli_word::word::PauliWord::new(n_qubits);
        word.xbits = self.xbits;
        word.zbits = self.zbits;
        word.rehash();
        word
    }

    /// The planes of a general `PauliWord`.
    pub fn from_pauli_word(word: &GeneralWord<C>) -> Self {
        Self {
            xbits: word.xbits,
            zbits: word.zbits,
        }
    }
}

impl<const C: usize> PartialEq for Word<C> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.xbits.data == other.xbits.data && self.zbits.data == other.zbits.data
    }
}

impl<const C: usize> Eq for Word<C> {}

/// Bit-lexicographic order of one plane from qubit 0 upward, as `BitArray`'s
/// `Ord` defines it, evaluated a chunk at a time.
#[inline(always)]
fn plane_cmp<const C: usize>(a: &[Chunk; C], b: &[Chunk; C]) -> std::cmp::Ordering {
    for (&x, &y) in a.iter().zip(b) {
        let differ = x ^ y;
        if differ != 0 {
            return if (x >> differ.trailing_zeros()) & 1 == 0 {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            };
        }
    }
    std::cmp::Ordering::Equal
}

impl<const C: usize> Ord for Word<C> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        plane_cmp(&self.xbits.data, &other.xbits.data)
            .then_with(|| plane_cmp(&self.zbits.data, &other.zbits.data))
    }
}

impl<const C: usize> PartialOrd for Word<C> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const C: usize> Hash for Word<C> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.digest());
    }
}

impl<const C: usize> std::fmt::Debug for Word<C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Word(x={:x?}, z={:x?})",
            self.xbits.data, self.zbits.data
        )
    }
}

/// `Err(TooManyQubits)` unless a `C`-chunk word holds `n_qubits`.
pub(crate) fn check_width<const C: usize>(n_qubits: usize) -> Result<(), Error> {
    if n_qubits > max_qubits::<C>() {
        return Err(Error::TooManyQubits {
            got: n_qubits,
            max: max_qubits::<C>(),
        });
    }
    Ok(())
}

/// Build a [`Word`] from a length-`n_qubits` slice of Pauli labels
/// (`0=I, 1=X, 2=Z, 3=Y` — the [`ppvm_traits::char::Pauli`] discriminants).
/// Sets all bits and rehashes once.
pub fn word_from_codes<const C: usize>(codes: &[u8]) -> Result<Word<C>, Error> {
    let n_qubits = codes.len();
    check_width::<C>(n_qubits)?;
    let mut w = Word::<C>::new(n_qubits);
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
    Ok(w)
}

/// Inverse of [`word_from_codes`]: write `n_qubits` Pauli labels into `out`.
pub fn codes_from_word<const C: usize>(w: &Word<C>, out: &mut [u8]) {
    debug_assert!(out.len() <= max_qubits::<C>());
    for (q, slot) in out.iter_mut().enumerate() {
        let xb = w.xbits[q] as u8;
        let zb = w.zbits[q] as u8;
        *slot = xb | (zb << 1);
    }
}

/// Parse a `"IXYZ..."` string into a [`Word`] together with the list of
/// qubits where the Pauli is non-identity (the term's support).
pub fn parse_pauli_string<const C: usize>(
    s: &str,
    n_qubits: usize,
) -> Result<(Word<C>, Vec<u32>), Error> {
    check_width::<C>(n_qubits)?;
    let chars: Vec<char> = s.chars().filter(|c| *c != '_').collect();
    if chars.len() != n_qubits {
        return Err(Error::WrongLength {
            expected: n_qubits,
            got: chars.len(),
        });
    }
    let mut w = Word::<C>::new(n_qubits);
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
    Ok((w, support))
}

/// Compute the support (non-identity qubits) of `w`.
pub(crate) fn word_support<const C: usize>(w: &Word<C>, out: &mut Vec<u32>) {
    out.clear();
    for chunk in 0..C {
        let mut bits = w.xbits.data[chunk] | w.zbits.data[chunk];
        while bits != 0 {
            out.push((chunk * CHUNK_BITS) as u32 + bits.trailing_zeros());
            bits &= bits - 1;
        }
    }
}

/// Compact 64-bit hash of a [`Word`], used as the key in cache-friendly
/// membership tables: an `FxHashMap<u64, ()>` over the basis has a working
/// set ~6× smaller than `FxHashMap<Word, ()>`. The hash mixes the word's
/// cached hash once through `FxHasher` and never touches the 32-byte
/// payload.
#[inline(always)]
pub(crate) fn word_hash<const C: usize>(w: &Word<C>) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = fxhash::FxHasher::default();
    w.hash(&mut h);
    h.finish()
}
