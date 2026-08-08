// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The packed [`LossyPauliWord`] struct, its inherent constructors/accessors and
//! loss writes, the read-only [`Word`] inspection impl (`Site = LossySite<Pauli>`),
//! the [`PauliBits`] sub-site mutation impl (with direct branch-key builders and
//! `is_lost` overridden), and the
//! structural [`PartialEq`]/[`Eq`]/[`Clone`]/[`Debug`]/[`Display`]/parsing that
//! all agree on the logical identity `(nqubits, X bits, Z bits, loss bits)`.
//!
//! Ported from `ppvm-pauli-word/src/loss/data.rs` (packed layout,
//! `weight`/`loss_weight` popcounts, symbol parsing) to keep the hot paths at
//! parity. Design: `word-data-structures.md` §"Lossy Pauli word" and §"Logical
//! Pauli model"; `traits-2-configuration-and-hashing.md` §"Representation types"
//! and §"Pauli algebra traits" (`PauliBits`).

use std::fmt;
use std::marker::PhantomData;
use std::sync::atomic::{AtomicU64, Ordering};

use bitvec::array::BitArray;
use ppvm_pauli_word_2::{DefaultStorage, PauliStorage};
use ppvm_traits_2::{LossySite, Pauli, PauliBits, Word};

/// A fixed-width **lossy** Pauli word: parallel packed `x`, `z`, and `loss` bit
/// planes.
///
/// Each qubit slot is either a present Pauli — encoded by `(x, z)` exactly as an
/// ordinary word (`I = (0,0)`, `X = (1,0)`, `Z = (0,1)`, `Y = (1,1)`) — or
/// **lost**, with canonical bits `(x, z, lost) = (0, 0, 1)`. Loss is exclusive
/// with the four Paulis (`word-data-structures.md` §"Logical Pauli model"). `A`
/// is the backing-storage blob (e.g. `u64` packs up to 64 qubits) and `H` is the
/// private internal digest algorithm; **both are private representation
/// parameters, neither exposed through `Word` or `Indexable`**.
///
/// # Canonical loss invariant
///
/// A lost site must be identity in its X/Z planes:
/// `lost[q] = 1  ⇒  xbits[q] = 0 ∧ zbits[q] = 0`
/// (`word-data-structures.md` §"Canonical loss invariant"). Every mutator upholds
/// it: [`LossyPauliWord::set_lost`] clears X/Z before setting loss;
/// [`PauliBits::set_x_bit`]/[`set_z_bit`](PauliBits::set_z_bit) clear the loss bit
/// when they set an X/Z bit true; parsing/`set` route through these. This makes
/// the physical encoding of a logical lossy word unique, so the full-blob
/// comparison in [`PartialEq`] is equivalent to comparing only the logical bits.
///
/// The structural identity is `(nqubits, logical X bits, logical Z bits, logical
/// loss bits)`. Equality and hashing exclude unused capacity, the caches, and the
/// `PhantomData` marker.
///
/// One `AtomicU64` lazily caches the finalized structural digest. A warm
/// `key_hash()` (the hashbrown map-lookup hot path) is therefore a single field
/// read. Mutations invalidate the cell through their exclusive `&mut self`;
/// relaxed atomics are enough because racing misses recompute the same pure
/// digest. Keeping only the consumed digest avoids copying and invalidating three
/// atomic cells on every gate-produced key.
///
/// # Examples
///
/// ```
/// use ppvm_lossy_pauli_word_2::LossyPauliWord;
/// use ppvm_traits_2::{LossySite, Pauli, Word};
///
/// let w: LossyPauliWord = "XLZL".into();
/// assert_eq!(w.n_sites(), 4);
/// assert_eq!(w.get(0), LossySite::Present(Pauli::X));
/// assert_eq!(w.get(1), LossySite::Lost);
/// assert_eq!(w.weight(), 4); // X, L, Z, L are all non-identity
/// assert_eq!(w.loss_weight(), 2); // two qubits are lost
/// assert!(w.is_lost(1));
/// ```
pub struct LossyPauliWord<A: PauliStorage = DefaultStorage, H = fxhash::FxBuildHasher> {
    /// X-bit plane (unused high bits are `0`; lost slots are `0`).
    pub(crate) xbits: BitArray<A>,
    /// Z-bit plane (unused high bits are `0`; lost slots are `0`).
    pub(crate) zbits: BitArray<A>,
    /// Loss-bit plane: `1` at index `q` marks qubit `q` lost.
    pub(crate) lbits: BitArray<A>,
    /// Number of qubits (logical width).
    pub(crate) nqubits: usize,
    /// Lazy cache for the finalized structural digest.
    pub(crate) hash_cache: AtomicU64,
    /// The private internal digest algorithm; `fn() -> H` keeps `Send + Sync`.
    pub(crate) _hasher: PhantomData<fn() -> H>,
}

impl<A: PauliStorage, H> LossyPauliWord<A, H> {
    /// Construct the identity word `I…I` on `nqubits` qubits (all planes zero).
    ///
    /// Constructors leave the caches empty (`word-data-structures.md`
    /// §"Invalidation rules").
    #[inline]
    pub fn new(nqubits: usize) -> Self {
        debug_assert!(
            nqubits <= 8 * std::mem::size_of::<A>(),
            "nqubits {nqubits} exceeds the {}-bit backing storage",
            8 * std::mem::size_of::<A>(),
        );
        Self {
            xbits: BitArray::ZERO,
            zbits: BitArray::ZERO,
            lbits: BitArray::ZERO,
            nqubits,
            hash_cache: AtomicU64::new(0),
            _hasher: PhantomData,
        }
    }

    /// Assemble from already-packed planes and a width. Callers must uphold the
    /// canonical-unused-bits and canonical-loss invariants; the caches are left
    /// empty. Crate-internal — used by parsing and the key column, which both
    /// preserve the invariants.
    #[inline]
    pub(crate) fn from_planes(
        xbits: BitArray<A>,
        zbits: BitArray<A>,
        lbits: BitArray<A>,
        nqubits: usize,
    ) -> Self {
        Self {
            xbits,
            zbits,
            lbits,
            nqubits,
            hash_cache: AtomicU64::new(0),
            _hasher: PhantomData,
        }
    }

    /// Clear the structural digest after an X/Z mutation.
    #[inline]
    pub(crate) fn invalidate_xz(&mut self) {
        *self.hash_cache.get_mut() = 0;
    }

    /// Clear the structural digest after a loss mutation.
    #[inline]
    pub(crate) fn invalidate_loss(&mut self) {
        *self.hash_cache.get_mut() = 0;
    }

    /// Whether qubit `q` is lost. Inherent per `word-data-structures.md`
    /// §"Loss-specific behavior"; [`PauliBits::is_lost`] delegates here.
    #[inline]
    pub fn is_lost(&self, q: usize) -> bool {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        self.lbits[q]
    }

    /// Mark qubit `q` lost. Clears the X/Z bits first, then sets the loss bit, so
    /// the canonical loss invariant `lost ⇒ X/Z identity` holds
    /// (`word-data-structures.md` §"Canonical loss invariant"). Inherent, not a
    /// trait method (§"Loss-specific behavior").
    #[inline]
    pub fn set_lost(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.xbits[q] || self.zbits[q] {
            // Marking a *nonidentity* site lost changes the X/Z content.
            self.xbits.set(q, false);
            self.zbits.set(q, false);
            self.invalidate_xz();
        }
        if !self.lbits[q] {
            self.lbits.set(q, true);
            self.invalidate_loss();
        }
    }

    /// Clear loss at qubit `q`, returning it to identity `I`. The X/Z component is
    /// preserved (a lost site already has X/Z identity), only the loss component
    /// is invalidated (`word-data-structures.md` §"Invalidation rules": "Clear
    /// loss to identity | preserve | invalidate"). Inherent (§"Loss-specific
    /// behavior").
    #[inline]
    pub fn clear_loss(&mut self, q: usize) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        if self.lbits[q] {
            self.lbits.set(q, false);
            self.invalidate_loss();
        }
    }

    /// Build the reset-channel branch with `q` marked lost.
    ///
    /// This copies the packed planes directly and starts the derived digest cold;
    /// unlike `clone` followed by `set_lost`, it does not copy an atomic cache
    /// that the mutation immediately invalidates.
    #[inline]
    pub fn with_lost(&self, q: usize) -> Self {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        let mut lbits = self.lbits;
        xbits.set(q, false);
        zbits.set(q, false);
        lbits.set(q, true);
        Self::from_planes(xbits, zbits, lbits, self.nqubits)
    }

    /// Write a whole [`LossySite`] at qubit `q`, upholding the canonical loss
    /// invariant: `Present(p)` clears the loss bit then writes `p`'s X/Z bits;
    /// `Lost` routes through [`set_lost`](Self::set_lost)
    /// (`word-data-structures.md` §"Canonical loss invariant").
    #[inline]
    pub fn set(&mut self, q: usize, site: LossySite<Pauli>) {
        debug_assert!(q < self.nqubits, "qubit {q} out of bounds");
        match site {
            LossySite::Present(p) => {
                if self.lbits[q] {
                    self.lbits.set(q, false);
                    self.invalidate_loss();
                }
                let (x, z) = match p {
                    Pauli::I => (false, false),
                    Pauli::X => (true, false),
                    Pauli::Z => (false, true),
                    Pauli::Y => (true, true),
                };
                self.xbits.set(q, x);
                self.zbits.set(q, z);
                self.invalidate_xz();
            }
            LossySite::Lost => self.set_lost(q),
        }
    }

    /// Number of lost qubits: a popcount over the loss plane. Ported verbatim
    /// from `ppvm-pauli-word`'s `loss_weight` (one `popcnt` per machine word).
    /// Inherent, not a trait method (`word-data-structures.md` §"Loss-specific
    /// behavior").
    #[inline]
    pub fn loss_weight(&self) -> usize {
        let ls: &[u8] = bytemuck::bytes_of(&self.lbits.data);

        let mut total: u32 = 0;
        let (mut i, n) = (0usize, ls.len());

        while i + 8 <= n {
            let l = u64::from_ne_bytes(ls[i..i + 8].try_into().unwrap());
            total += l.count_ones();
            i += 8;
        }
        if i + 4 <= n {
            let l = u32::from_ne_bytes(ls[i..i + 4].try_into().unwrap());
            total += l.count_ones();
            i += 4;
        }
        if i + 2 <= n {
            let l = u16::from_ne_bytes(ls[i..i + 2].try_into().unwrap());
            total += l.count_ones();
            i += 2;
        }
        if i < n {
            total += ls[i].count_ones();
        }

        total as usize
    }
}

impl<A: PauliStorage, H> Word for LossyPauliWord<A, H> {
    type Site = LossySite<Pauli>;

    #[inline]
    fn n_sites(&self) -> usize {
        self.nqubits
    }

    #[inline]
    fn get(&self, index: usize) -> LossySite<Pauli> {
        assert!(index < self.nqubits, "index {index} out of bounds");
        if self.lbits[index] {
            return LossySite::Lost;
        }
        let p = match (self.xbits[index], self.zbits[index]) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (false, true) => Pauli::Z,
            (true, true) => Pauli::Y,
        };
        LossySite::Present(p)
    }

    /// Number of non-identity factors — `X`, `Y`, `Z`, **and** `Lost`: a fused
    /// popcount of `x | z | loss` over the packed planes. Ported verbatim from
    /// `ppvm-pauli-word`'s lossy `weight` to keep it at parity.
    #[inline]
    fn weight(&self) -> usize {
        let xs: &[u8] = bytemuck::bytes_of(&self.xbits.data);
        let zs: &[u8] = bytemuck::bytes_of(&self.zbits.data);
        let ls: &[u8] = bytemuck::bytes_of(&self.lbits.data);
        debug_assert_eq!(xs.len(), zs.len());
        debug_assert_eq!(xs.len(), ls.len());

        let mut total: u32 = 0;
        let (mut i, n) = (0usize, xs.len());

        while i + 8 <= n {
            let x = u64::from_ne_bytes(xs[i..i + 8].try_into().unwrap());
            let z = u64::from_ne_bytes(zs[i..i + 8].try_into().unwrap());
            let l = u64::from_ne_bytes(ls[i..i + 8].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 8;
        }
        if i + 4 <= n {
            let x = u32::from_ne_bytes(xs[i..i + 4].try_into().unwrap());
            let z = u32::from_ne_bytes(zs[i..i + 4].try_into().unwrap());
            let l = u32::from_ne_bytes(ls[i..i + 4].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 4;
        }
        if i + 2 <= n {
            let x = u16::from_ne_bytes(xs[i..i + 2].try_into().unwrap());
            let z = u16::from_ne_bytes(zs[i..i + 2].try_into().unwrap());
            let l = u16::from_ne_bytes(ls[i..i + 2].try_into().unwrap());
            total += (x | z | l).count_ones();
            i += 2;
        }
        if i < n {
            total += (xs[i] | zs[i] | ls[i]).count_ones();
        }

        total as usize
    }

    #[inline]
    fn iter(&self) -> impl Iterator<Item = LossySite<Pauli>> {
        (0..self.nqubits).map(move |i| self.get(i))
    }
}

impl<A: PauliStorage, H> PauliBits for LossyPauliWord<A, H> {
    const PREFER_BORROWED_REKEY: bool = true;

    #[inline]
    fn x_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits[i]
    }

    #[inline]
    fn z_bit(&self, i: usize) -> bool {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits[i]
    }

    /// Set the X bit at `i`. To keep loss exclusive with the Paulis, setting the
    /// bit *true* on a lost site clears the loss bit (the "Replace loss with
    /// Pauli" invalidation row of `word-data-structures.md` §"Invalidation
    /// rules"). Only in-range slots are touched, preserving the
    /// canonical-unused-bits invariant.
    #[inline]
    fn set_x_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits.set(i, v);
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss();
        }
        self.invalidate_xz();
    }

    #[inline]
    fn set_z_bit(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.zbits.set(i, v);
        if v && self.lbits[i] {
            self.lbits.set(i, false);
            self.invalidate_loss();
        }
        self.invalidate_xz();
    }

    #[inline]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        self.xbits.set(i, x);
        self.zbits.set(i, z);
        if (x || z) && self.lbits[i] {
            self.lbits.set(i, false);
        }
        self.invalidate_xz();
    }

    #[inline]
    fn set_xz_bits2(&mut self, i: usize, xi: bool, zi: bool, j: usize, xj: bool, zj: bool) {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        self.xbits.set(i, xi);
        self.zbits.set(i, zi);
        self.xbits.set(j, xj);
        self.zbits.set(j, zj);
        if (xi || zi) && self.lbits[i] {
            self.lbits.set(i, false);
        }
        if (xj || zj) && self.lbits[j] {
            self.lbits.set(j, false);
        }
        self.invalidate_xz();
    }

    /// Build a rotation branch directly from the packed planes. This avoids
    /// cloning an atomic cache only to invalidate it on the first bit write.
    #[inline]
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        if toggle_x {
            let bit = xbits[i];
            xbits.set(i, !bit);
        }
        if toggle_z {
            let bit = zbits[i];
            zbits.set(i, !bit);
        }
        Self::from_planes(xbits, zbits, self.lbits, self.nqubits)
    }

    /// Build a two-site rotation branch with one packed-plane copy.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn toggled_bits2(
        &self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self {
        debug_assert!(i < self.nqubits && j < self.nqubits, "index out of bounds");
        let mut xbits = self.xbits;
        let mut zbits = self.zbits;
        if toggle_x_i {
            let bit = xbits[i];
            xbits.set(i, !bit);
        }
        if toggle_z_i {
            let bit = zbits[i];
            zbits.set(i, !bit);
        }
        if toggle_x_j {
            let bit = xbits[j];
            xbits.set(j, !bit);
        }
        if toggle_z_j {
            let bit = zbits[j];
            zbits.set(j, !bit);
        }
        Self::from_planes(xbits, zbits, self.lbits, self.nqubits)
    }

    /// A lossy word reports genuine loss (overriding the `false` default), so
    /// loss-aware generic propagation and the Clifford loss-guard see it.
    #[inline]
    fn is_lost(&self, i: usize) -> bool {
        LossyPauliWord::is_lost(self, i)
    }

    /// Delegates to the inherent [`LossyPauliWord::set_lost`] — clears the site's
    /// X/Z bits, then sets the loss bit (canonical loss invariant).
    #[inline]
    fn set_lost(&mut self, i: usize) {
        LossyPauliWord::set_lost(self, i);
    }

    /// Delegates to the inherent [`LossyPauliWord::clear_loss`].
    #[inline]
    fn clear_lost(&mut self, i: usize) {
        LossyPauliWord::clear_loss(self, i);
    }

    /// Build a loss-channel recovery branch directly from the packed planes.
    #[inline]
    fn loss_cleared(&self, i: usize) -> Self {
        debug_assert!(i < self.nqubits, "index {i} out of bounds");
        let mut lbits = self.lbits;
        lbits.set(i, false);
        Self::from_planes(self.xbits, self.zbits, lbits, self.nqubits)
    }

    /// Delegates to the inherent [`LossyPauliWord::loss_weight`] — the popcount
    /// over the loss plane the `MaxLossWeight` policy thresholds.
    #[inline]
    fn loss_weight(&self) -> usize {
        LossyPauliWord::loss_weight(self)
    }
}

impl<A: PauliStorage, H> Clone for LossyPauliWord<A, H> {
    /// Cloning copies the cached finalized digest because the clone has identical
    /// structural contents. The remaining relaxed atomic load is intentional:
    /// preserving a warm lookup on an unmodified clone is part of the lazy-cache
    /// performance contract. Branch builders above avoid this cost only when they
    /// immediately change the key. Hand-written so no spurious `H: Clone` bound
    /// is imposed.
    #[inline]
    fn clone(&self) -> Self {
        Self {
            xbits: self.xbits,
            zbits: self.zbits,
            lbits: self.lbits,
            nqubits: self.nqubits,
            hash_cache: AtomicU64::new(self.hash_cache.load(Ordering::Relaxed)),
            _hasher: PhantomData,
        }
    }
}

impl<A: PauliStorage, H> fmt::Debug for LossyPauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LossyPauliWord")
            .field("nqubits", &self.nqubits)
            .field("word", &self.to_string())
            .finish()
    }
}

/// Structural equality over `(nqubits, X bits, Z bits, loss bits)`; the caches and
/// the `PhantomData` marker are excluded. The canonical invariants make the
/// full-blob comparison equivalent to comparing only the logical bits.
impl<A: PauliStorage, H> PartialEq for LossyPauliWord<A, H> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.nqubits == other.nqubits
            && self.xbits.data == other.xbits.data
            && self.zbits.data == other.zbits.data
            && self.lbits.data == other.lbits.data
    }
}

impl<A: PauliStorage, H> Eq for LossyPauliWord<A, H> {}

impl<A: PauliStorage, H> fmt::Display for LossyPauliWord<A, H> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nqubits {
            let c = match self.get(i) {
                LossySite::Present(Pauli::I) => 'I',
                LossySite::Present(Pauli::X) => 'X',
                LossySite::Present(Pauli::Y) => 'Y',
                LossySite::Present(Pauli::Z) => 'Z',
                LossySite::Lost => 'L',
            };
            f.write_str(c.encode_utf8(&mut [0u8; 1]))?;
        }
        Ok(())
    }
}

impl<A: PauliStorage, H> From<&str> for LossyPauliWord<A, H> {
    /// Parse a lossy Pauli string of `I`/`X`/`Y`/`Z`/`L` symbols (underscores are
    /// ignored separators), mirroring `ppvm-pauli-word`'s lossy `From<&str>`.
    /// Panics on any other character or if the width exceeds the backing storage.
    fn from(value: &str) -> Self {
        let mut xbits = BitArray::<A>::ZERO;
        let mut zbits = BitArray::<A>::ZERO;
        let mut lbits = BitArray::<A>::ZERO;
        let mut i = 0usize;
        for ch in value.chars() {
            match ch {
                'I' => {}
                'X' => xbits.set(i, true),
                'Z' => zbits.set(i, true),
                'Y' => {
                    xbits.set(i, true);
                    zbits.set(i, true);
                }
                'L' => lbits.set(i, true),
                '_' => continue,
                other => panic!("invalid lossy Pauli character: {other}"),
            }
            i += 1;
        }
        Self::from_planes(xbits, zbits, lbits, i)
    }
}

impl<A: PauliStorage, H> From<String> for LossyPauliWord<A, H> {
    #[inline]
    fn from(value: String) -> Self {
        Self::from(value.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_get_roundtrip() {
        let w: LossyPauliWord = "XLZL".into();
        assert_eq!(w.n_sites(), 4);
        assert_eq!(w.get(0), LossySite::Present(Pauli::X));
        assert_eq!(w.get(1), LossySite::Lost);
        assert_eq!(w.get(2), LossySite::Present(Pauli::Z));
        assert_eq!(w.get(3), LossySite::Lost);
        assert_eq!(w.to_string(), "XLZL");
    }

    #[test]
    fn all_symbols_roundtrip() {
        let w: LossyPauliWord = "IXYZL".into();
        assert_eq!(w.to_string(), "IXYZL");
        let via_iter: Vec<LossySite<Pauli>> = w.iter().collect();
        let via_get: Vec<LossySite<Pauli>> = (0..w.n_sites()).map(|i| w.get(i)).collect();
        assert_eq!(via_iter, via_get);
    }

    #[test]
    fn weight_and_loss_weight() {
        let w: LossyPauliWord = "IIII".into();
        assert_eq!(w.weight(), 0);
        assert_eq!(w.loss_weight(), 0);

        let w: LossyPauliWord = "XLIL".into();
        assert_eq!(w.weight(), 3); // X, L, L are non-identity
        assert_eq!(w.loss_weight(), 2);

        let w: LossyPauliWord = "LLLL".into();
        assert_eq!(w.weight(), 4);
        assert_eq!(w.loss_weight(), 4);
    }

    #[test]
    fn underscore_is_ignored() {
        let a: LossyPauliWord = "X_L_Z".into();
        let b: LossyPauliWord = "XLZ".into();
        assert_eq!(a, b);
    }

    #[test]
    fn set_lost_clears_xz_and_is_exclusive() {
        let mut w: LossyPauliWord = "XYZ".into();
        w.set_lost(0); // was X (nonidentity)
        assert_eq!(w.get(0), LossySite::Lost);
        assert!(!w.x_bit(0) && !w.z_bit(0), "lost site must be X/Z identity");
        assert_eq!(w.to_string(), "LYZ");
    }

    #[test]
    fn with_lost_builds_canonical_branch_without_mutating_source() {
        let w: LossyPauliWord = "XYZ".into();
        let lost = w.with_lost(1);
        assert_eq!(w.to_string(), "XYZ");
        assert_eq!(lost.to_string(), "XLZ");
        assert!(!lost.x_bit(1) && !lost.z_bit(1));
    }

    #[test]
    fn clear_loss_returns_identity() {
        let mut w: LossyPauliWord = "XLZ".into();
        w.clear_loss(1);
        assert_eq!(w.get(1), LossySite::Present(Pauli::I));
        assert_eq!(w.loss_weight(), 0);
        assert_eq!(w.to_string(), "XIZ");
    }

    #[test]
    fn set_bit_true_clears_loss() {
        // Canonical invariant: a present X/Z bit and loss are exclusive.
        let mut w: LossyPauliWord = "L".into();
        w.set_x_bit(0, true);
        assert!(!w.is_lost(0), "setting an X bit clears loss");
        assert_eq!(w.get(0), LossySite::Present(Pauli::X));
    }

    #[test]
    fn set_writes_sites() {
        let mut w: LossyPauliWord = LossyPauliWord::new(3);
        w.set(0, LossySite::Present(Pauli::Y));
        w.set(1, LossySite::Lost);
        w.set(2, LossySite::Present(Pauli::Z));
        assert_eq!(w.to_string(), "YLZ");
        // Overwrite a lost site with a Pauli and back.
        w.set(1, LossySite::Present(Pauli::X));
        assert_eq!(w.to_string(), "YXZ");
        assert!(!w.is_lost(1));
    }

    #[test]
    fn equality_excludes_width_mismatch() {
        let a: LossyPauliWord = "XL".into();
        let b: LossyPauliWord = "XLI".into();
        assert_ne!(a, b, "different widths are structurally distinct");
    }

    #[test]
    fn send_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<LossyPauliWord>();
    }
}
