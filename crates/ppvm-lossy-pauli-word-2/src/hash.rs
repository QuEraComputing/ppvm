use std::hash::{BuildHasher, Hash, Hasher};
use std::sync::atomic::Ordering;

use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::Indexable;

use crate::LossyPauliWord;

#[inline]
pub(crate) fn structural_hash_lossy<A, H>(x: &A, z: &A, l: &A, _nqubits: usize) -> u64
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    let mut hasher = H::default().build_hasher();
    x.hash(&mut hasher);
    z.hash(&mut hasher);
    l.hash(&mut hasher);
    let structural = H::finalize_hash(hasher.finish(), std::mem::size_of::<A>());
    // Apply the legacy map's index transform before the digest reaches
    // IdentityBuildHasher, which passes it through without further hashing.
    H::index_hash(structural)
}

/// Write the finalized word digest as one `u64` for `IdentityHasher`
/// to pass through unchanged.
impl<A, H> Hash for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn hash<S: Hasher>(&self, state: &mut S) {
        state.write_u64(self.key_hash());
    }
}

/// Cache the X/Z and loss digests independently, computing each on demand.
/// Combine both digests and finalize the word hash on every call.
/// Zero marks an empty cache; a zero digest is recomputed on the next call.
impl<A, H> Indexable for LossyPauliWord<A, H>
where
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
{
    #[inline]
    fn key_hash(&self) -> u64 {
        let xz_hash = match self.xz_hash_cache.load(Ordering::Relaxed) {
            0 => {
                let mut hasher = H::default().build_hasher();
                self.xbits.data.hash(&mut hasher);
                self.zbits.data.hash(&mut hasher);
                let digest = hasher.finish();
                self.xz_hash_cache.store(digest, Ordering::Relaxed);
                digest
            },
            cached => cached
        };
        let loss_hash = match self.loss_hash_cache.load(Ordering::Relaxed) {
            0 => {
                let mut hasher = H::default().build_hasher();
                self.lbits.data.hash(&mut hasher);
                let digest = hasher.finish();
                self.loss_hash_cache.store(digest, Ordering::Relaxed);
                digest
            },
            cached => cached
        };
        let mut hasher = H::default().build_hasher();
        hasher.write_u64(xz_hash);
        hasher.write_u64(loss_hash);

        let digest = H::finalize_hash(hasher.finish(), std::mem::size_of::<[u64;2]>(),);
        H::index_hash(digest)
    }
}
