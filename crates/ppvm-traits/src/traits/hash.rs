// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Per-hasher finalization of a Pauli word's cached key hash.

/// Post-processing that a `PauliWord`'s hasher
/// applies to its raw 64-bit digest before it is cached as the map-key hash.
///
/// The cached value is what `hashbrown` ultimately splits into a bucket index
/// (low bits) and a control-byte tag (top 7 bits). Whether the raw digest is
/// good enough for that split is a property of the **hasher**, not of the
/// Pauli word: an AES-based hasher such as `gxhash` avalanches well even for an
/// 8-byte key, whereas `FxHasher` — a couple of multiply-rotate rounds — leaves
/// the low bits of a short key highly correlated and needs a fold to fix them.
///
/// Keeping this on the hasher is what makes the abstraction correct. An earlier
/// version folded based on storage width alone, which conflated "narrow
/// storage" with "weak hasher" and so wrongly folded `gxhash` too. Here each
/// hasher declares how (if at all) its output must be adjusted, and
/// `PauliWord::rehash` just defers to it. The default
/// is the identity — the right choice for any hasher that already distributes
/// its low bits well — so a custom hasher opts in with a bare
/// `impl HashFinalize for MyHasher {}`.
pub trait HashFinalize {
    /// Finalize `raw` (the value returned by `Hasher::finish`) for a key whose
    /// backing storage is `storage_bytes` wide per bit-array.
    ///
    /// The width is supplied so a hasher can fold only for the short keys that
    /// actually need it and pass wider keys through. It is a compile-time
    /// constant at every call site (`size_of` of the storage), so any branch
    /// on it is monomorphized away.
    #[inline(always)]
    fn finalize_hash(raw: u64, storage_bytes: usize) -> u64 {
        let _ = storage_bytes;
        raw
    }
}

impl HashFinalize for fxhash::FxBuildHasher {
    /// `FxHasher` consumes one `u64` word at a time and its avalanche is weak
    /// for short inputs: when a word fits in a single `u64` per bit-array
    /// (`[u8; 8]` and narrower) it goes through only a couple of
    /// multiply-rotate rounds, leaving the low bits — the ones `hashbrown`
    /// uses to choose a bucket — correlated. At high fill that clusters
    /// distinct words into a few oversized buckets and the `insert_unique`
    /// probe length explodes (~7x at 64 qubits in `[u8; 8]`; max bucket 2257
    /// vs an ideal of ~6 — see `examples/trotter_storage_cliff.rs`). Folding
    /// the high half into the low half decorrelates them.
    ///
    /// Wider storage gets enough rounds to distribute the low bits already, so
    /// it is passed through: folding there would only mix the top bits (which
    /// `hashbrown` reserves for its tag) back into the bucket index, coupling
    /// two values it wants independent (~4-6% slower in the same benchmark).
    #[inline(always)]
    fn finalize_hash(raw: u64, storage_bytes: usize) -> u64 {
        if storage_bytes <= std::mem::size_of::<u64>() {
            raw ^ (raw >> 32)
        } else {
            raw
        }
    }
}

// `gxhash` is AES-based and avalanches well even for an 8-byte key (measured
// max bucket 6, essentially ideal, vs fxhash's 2257 at 64 qubits), so the
// identity default is exactly right — folding would only pay the tag/bucket
// coupling cost above with no distribution benefit.
#[cfg(feature = "gxhash")]
impl HashFinalize for gxhash::GxBuildHasher {}

#[cfg(test)]
mod tests {
    use super::*;

    // A digest whose high and low halves differ, so a fold is observable.
    const RAW: u64 = 0xDEAD_BEEF_0000_0001;

    #[test]
    fn fxhash_folds_narrow_storage() {
        // `[u8; 8]` and narrower fit in one u64 per bit-array → fold.
        for width in [1, 2, 4, 8] {
            assert_eq!(
                <fxhash::FxBuildHasher as HashFinalize>::finalize_hash(RAW, width),
                RAW ^ (RAW >> 32),
                "fxhash should fold at width {width}"
            );
        }
    }

    #[test]
    fn fxhash_passes_wide_storage_through() {
        // `[u8; 16]` and up are already well distributed → identity.
        for width in [16, 32, 64] {
            assert_eq!(
                <fxhash::FxBuildHasher as HashFinalize>::finalize_hash(RAW, width),
                RAW,
                "fxhash should not fold at width {width}"
            );
        }
    }

    #[cfg(feature = "gxhash")]
    #[test]
    fn gxhash_never_folds() {
        // gxhash already distributes its low bits, so it is the identity at
        // every width — including the narrow ones where fxhash folds.
        for width in [1, 2, 4, 8, 16, 32, 64] {
            assert_eq!(
                <gxhash::GxBuildHasher as HashFinalize>::finalize_hash(RAW, width),
                RAW,
                "gxhash should never fold (width {width})"
            );
        }
    }
}
