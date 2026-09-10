// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

/// The action a fermionic factor performs on its mode.
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FermionSite {
    /// The physical mode this factor acts on.
    pub mode: usize,
    /// Whether this factor creates or annihilates.
    pub action: FermionAction,
}
