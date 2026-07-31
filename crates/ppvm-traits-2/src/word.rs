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
pub trait PauliBits: Word<Site = Pauli> {
    /// Read the X bit at index `i`.
    fn x_bit(&self, i: usize) -> bool;
    /// Read the Z bit at index `i`.
    fn z_bit(&self, i: usize) -> bool;
    /// Set the X bit at index `i` (invalidates the hash lazily).
    fn set_x_bit(&mut self, i: usize, v: bool);
    /// Set the Z bit at index `i` (invalidates the hash lazily).
    fn set_z_bit(&mut self, i: usize, v: bool);
    /// Whether index `i` is lost; `LossyPauliWord` overrides this.
    fn is_lost(&self, i: usize) -> bool {
        let _ = i;
        false
    }
}
