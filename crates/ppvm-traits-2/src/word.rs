// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The read-only [`Word`] inspection interface, its site-alphabet leaf types,
//! and the Pauli-specific bit-mutation capability [`PauliBits`].
//!
//! Design: `traits-2-configuration-and-hashing.md` §"Representation types" and
//! §"Pauli algebra traits". `Word` is deliberately *not* the propagation
//! interface — that is the sub-site bit algebra on [`PauliBits`] /
//! [`crate::pauli::SymplecticColumns`].

/// The common **read-only** concept for an indexed algebraic monomial.
///
/// An *inspection* interface — extent, per-site read, weight, and iteration —
/// consumed by display, serialization, tests, and the sparse-sum plumbing. It
/// carries no mutation, so an ordered algebra (a normal-ordered fermionic
/// product) can implement it honestly; structural mutation is relocated to each
/// algebra's own traits.
///
/// Design: §"Representation types".
pub trait Word {
    /// The operator alphabet at one index.
    type Site;

    /// Number of sites (for a dense Pauli word, the qubit width).
    fn n_sites(&self) -> usize;

    /// Read the site at `index`.
    fn get(&self, index: usize) -> Self::Site;

    /// Number of non-identity factors according to the concrete site alphabet.
    ///
    /// A Pauli-motivated read (the `MaxPauliWeight` policy needs it); an ordered
    /// representation that stores no explicit identities may have
    /// `weight() == n_sites()`.
    fn weight(&self) -> usize;

    /// Iterate the sites in index order.
    fn iter(&self) -> impl Iterator<Item = Self::Site>;
}

/// A single-qubit Pauli symbol — the site alphabet of an ordinary packed Pauli
/// word (`Word<Site = Pauli>`).
///
/// Design: §"Representation types".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Pauli {
    /// Identity `I`.
    I,
    /// Pauli `X`.
    X,
    /// Pauli `Y`.
    Y,
    /// Pauli `Z`.
    Z,
}

/// A site that may be lost — the alphabet of a lossy packed word
/// (`Word<Site = LossySite<Pauli>>`).
///
/// Design: §"Representation types".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LossySite<S> {
    /// A present site carrying `S`.
    Present(S),
    /// A lost site — the qubit no longer participates.
    Lost,
}

/// The action a fermionic factor performs on its mode.
///
/// Design: §"Representation types" (the `FermionAction` carried by
/// [`FermionSite`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FermionAction {
    /// A creation operator `a†`.
    Create,
    /// An annihilation operator `a`.
    Annihilate,
}

/// One factor of an ordered fermionic product — the alphabet of a future
/// ordered fermionic word (`Word<Site = FermionSite>`), whose index denotes
/// factor order while the site carries the physical mode.
///
/// Design: §"Representation types".
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FermionSite {
    /// The physical mode this factor acts on.
    pub mode: usize,
    /// Whether this factor creates or annihilates.
    pub action: FermionAction,
}

/// Mutable single-vector X/Z access — a point of `GF(2)^{2n}`.
///
/// Hosts the rotation/branching kernels, which flip individual bits and ship
/// the sign to the coefficient. The narrow bit-level slice of the retired
/// `PauliWordTrait`, separated from key identity ([`crate::hash::Indexable`]),
/// inspection ([`Word`]), and phase. Implemented by `PauliWord` and
/// `LossyPauliWord` (in `ppvm-pauli-word-2`).
///
/// Design: §"Pauli algebra traits" (`PauliBits`). The branch these kernels
/// stage — `c·P → cos·c·P + sin·c·(iGP)` — has a genuinely new key
/// (`lean/PPVM/Instantiations/Rotation.lean` `anticommute_new_key`).
///
/// # Friction: supertrait is `Word`, not `Word<Site = Pauli>`
///
/// The design sketches `PauliBits: Word<Site = Pauli>`, but `LossyPauliWord`
/// implements `PauliBits` while its `Word::Site` is `LossySite<Pauli>` (a lost
/// site is not a bare `Pauli`), so it cannot satisfy a `Word<Site = Pauli>`
/// supertrait. The bound is relaxed to `Word` — this trait's own methods are
/// alphabet-agnostic (they read/write raw X/Z bits and a loss flag), and any
/// code that genuinely needs `Pauli` propagation re-adds `Word<Site = Pauli>`
/// on *its* own methods (the pattern already noted for the propagation layer in
/// `graded.rs`). `PauliWord` still implements `Word<Site = Pauli>`, so nothing
/// on the ordinary path changes.
pub trait PauliBits: Word {
    /// Whether two-site sum re-keying should build from a borrowed source word.
    ///
    /// Copy-like words keep the owned mutation path. Representations whose
    /// `Clone` carries synchronization cost can opt into direct branch builders.
    const PREFER_BORROWED_REKEY: bool = false;

    /// Read the X bit at index `i`.
    fn x_bit(&self, i: usize) -> bool;
    /// Read the Z bit at index `i`.
    fn z_bit(&self, i: usize) -> bool;
    /// Set the X bit at index `i` (invalidates the hash lazily).
    fn set_x_bit(&mut self, i: usize, v: bool);
    /// Set the Z bit at index `i` (invalidates the hash lazily).
    fn set_z_bit(&mut self, i: usize, v: bool);
    /// Set both bit planes at one site.
    ///
    /// The default composes the scalar setters. Packed words whose setters
    /// refresh structural metadata eagerly may override this to refresh once.
    #[inline(always)]
    fn set_xz_bits(&mut self, i: usize, x: bool, z: bool) {
        self.set_x_bit(i, x);
        self.set_z_bit(i, z);
    }
    /// Whether index `i` is lost; `LossyPauliWord` overrides this.
    #[inline(always)]
    fn is_lost(&self, i: usize) -> bool {
        let _ = i;
        false
    }

    /// Packed local Pauli code: `0=I, 1=X, 2=Z, 3=Y`.
    #[inline(always)]
    fn pauli_code(&self, i: usize) -> u8 {
        (self.x_bit(i) as u8) | ((self.z_bit(i) as u8) << 1)
    }

    /// Mark index `i` lost. **No-op by default** — a word with no loss component
    /// cannot represent loss, and old's shared loss kernels reach this arm only
    /// under an `is_lost` guard that is a const `false` there.
    ///
    /// `LossyPauliWord` overrides it (clearing the site's X/Z bits first, per the
    /// canonical loss invariant). Keeping the mutator on `PauliBits` beside the
    /// const-`false` [`is_lost`](PauliBits::is_lost) is what lets **one** loss
    /// kernel serve both word types at zero cost for the non-lossy one
    /// (architecture feature 11): every `if k.is_lost(q) { … }` branch folds away
    /// at monomorphization, taking the `set_lost`/`clear_lost` calls inside it with
    /// it.
    fn set_lost(&mut self, i: usize) {
        let _ = i;
    }

    /// Clear the loss flag at index `i`, returning the site to identity.
    /// **No-op by default**; `LossyPauliWord` overrides it. See
    /// [`set_lost`](PauliBits::set_lost).
    fn clear_lost(&mut self, i: usize) {
        let _ = i;
    }

    /// A copy of this word with the loss flag at `i` cleared.
    ///
    /// The default preserves compatibility for custom word types by cloning and
    /// calling [`clear_lost`](PauliBits::clear_lost). Packed lossy words can
    /// override it to copy their planes directly, avoiding an atomic cache copy
    /// that the loss mutation would immediately invalidate.
    #[inline]
    fn loss_cleared(&self, i: usize) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        out.clear_lost(i);
        out
    }

    /// A copy of this word with the X and/or Z bit at `i` toggled — the
    /// **rotation-branch key builder** (`iGP` from a diagonal `P`).
    ///
    /// Provided as clone-then-flip so every `PauliBits` implementer gets a branch
    /// builder for free and the rotation/branching kernels can be generic over the
    /// word type (the ordinary and the lossy key run the *same* kernel — see
    /// `ppvm-pauli-sum-2`'s rotation and loss modules). `PauliWord` overrides it
    /// with a direct plane copy that leaves the digest cache empty, skipping the
    /// wasted cache load+invalidate a `clone` + `set_*_bit` pair performs.
    fn toggled_bits(&self, i: usize, toggle_x: bool, toggle_z: bool) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        if toggle_x {
            let b = out.x_bit(i);
            out.set_x_bit(i, !b);
        }
        if toggle_z {
            let b = out.z_bit(i);
            out.set_z_bit(i, !b);
        }
        out
    }

    /// A copy of this word with the X and/or Z bits at **two** sites toggled —
    /// the **two-qubit** rotation-branch key builder (`iG_aG_b·P`).
    ///
    /// Chaining [`toggled_bits`](PauliBits::toggled_bits) twice would build two
    /// whole words per produced branch term, and for a packed word each of those is
    /// a full copy of *both* bit planes plus a word rebuild — so the intermediate
    /// is pure waste on the hot path of `rzz`/`rxx`/`ryy`/`rotate_2` (old built
    /// **one** `k.clone()` and then wrote four bits into it,
    /// `ppvm-pauli-sum/src/sum/rot2.rs`). The two-site entry point makes the single
    /// copy the *only* copy, and it scales with the storage tier: at `[u8; 32]` the
    /// chained form moved 64 redundant bytes per branch term.
    ///
    /// The default is clone-then-flip (one clone, up to four bit writes), which is
    /// already old's shape; `PauliWord` overrides it with a direct plane copy that
    /// leaves the digest cache empty, exactly as it does for the single-site form.
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
    ) -> Self
    where
        Self: Sized + Clone,
    {
        let mut out = self.clone();
        if toggle_x_i {
            let b = out.x_bit(i);
            out.set_x_bit(i, !b);
        }
        if toggle_z_i {
            let b = out.z_bit(i);
            out.set_z_bit(i, !b);
        }
        if toggle_x_j {
            let b = out.x_bit(j);
            out.set_x_bit(j, !b);
        }
        if toggle_z_j {
            let b = out.z_bit(j);
            out.set_z_bit(j, !b);
        }
        out
    }

    /// Consume a word and toggle X/Z bits at two sites in place.
    ///
    /// Re-keying kernels already own their key, so this avoids the structural
    /// copy required by [`toggled_bits2`](PauliBits::toggled_bits2). Packed words
    /// may override it to defer derived-hash refresh until all writes complete.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn into_toggled_bits2(
        mut self,
        i: usize,
        toggle_x_i: bool,
        toggle_z_i: bool,
        j: usize,
        toggle_x_j: bool,
        toggle_z_j: bool,
    ) -> Self
    where
        Self: Sized,
    {
        if toggle_x_i {
            let b = self.x_bit(i);
            self.set_x_bit(i, !b);
        }
        if toggle_z_i {
            let b = self.z_bit(i);
            self.set_z_bit(i, !b);
        }
        if toggle_x_j {
            let b = self.x_bit(j);
            self.set_x_bit(j, !b);
        }
        if toggle_z_j {
            let b = self.z_bit(j);
            self.set_z_bit(j, !b);
        }
        self
    }

    /// Number of lost sites. **Zero by default**; `LossyPauliWord` overrides it
    /// with a popcount over its loss plane. This is the `MaxLossWeight` truncation
    /// predicate (old's `PauliWordTrait::loss_weight`), which is why it belongs on
    /// the read side of the same trait rather than on the lossy word alone: a
    /// policy is generic over the key.
    fn loss_weight(&self) -> usize {
        0
    }

    /// Whether this word anticommutes with the single-qubit Pauli
    /// `pauli = (x_bit, z_bit)` at index `i`, i.e. whether the symplectic form
    /// `ω(P, Q) = x_P·z_Q ⊕ z_P·x_Q` is `1` there.
    ///
    /// The old `PauliWordTrait::anticommutes_at`
    /// (`ppvm-traits/src/traits/word_trait.rs`), reproduced verbatim as a
    /// provided method — it is the pivot test the tableau measurement search
    /// runs (`ppvm-tableau/src/data.rs`), and it is derivable from the two bit
    /// reads, so it needs no new required method.
    #[inline]
    fn anticommutes_at(&self, i: usize, pauli: (bool, bool)) -> bool {
        (self.x_bit(i) & pauli.1) ^ (self.z_bit(i) & pauli.0)
    }
}
