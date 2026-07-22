// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::char::Pauli;
use std::fmt::Display;
use std::hash::Hash;

/// Iterate over a Pauli word slot-by-slot.
pub trait PauliIter {
    /// Yield the [`Pauli`] symbol at each qubit position, in order.
    fn iter(&self) -> impl Iterator<Item = Pauli>;
}

/// Word-level Pauli operations. Types implementing this trait automatically
/// gain [`crate::traits::Clifford`] and [`crate::traits::CliffordExtensions`]
/// behavior via blanket impls in the `crate::traits::clifford` module, using
/// `get_lbit` / `get_xbit` / `set_xbit` / `set_zbit` / `rehash` to transform
/// Pauli words.
///
/// # Word-level vs. state-level gate semantics
///
/// The blanket `Clifford` impl applies gates at the *bit* level of a single
/// Pauli word. X / Y / Z are bit-level no-ops on a Pauli word — they affect
/// phase, not the X/Z bits — so `word.x(i)`, `word.y(i)`, `word.z(i)` are
/// deliberately silent. Phase is tracked separately by
/// `PhasedPauliWord`, which implements `Clifford` manually.
///
/// If you need a word representation whose gate behavior is *not* pure
/// bit manipulation (phase tracking, fused multi-qubit updates, alternative
/// loss semantics), do **not** implement `PauliWordTrait` on it — define a
/// specialized `Clifford` impl instead, the way `PhasedPauliWord` does.
/// Implementing both `PauliWordTrait` and a custom `impl Clifford for ...`
/// for the same type will not compile: the blanket impl overlaps with your
/// custom impl (coherence error), not silently shadows it.
pub trait PauliWordTrait: Clone + Hash + Eq + PauliIter + From<String> + Display {
    /// Construct an identity word over `nqubits` qubits.
    fn new(nqubits: usize) -> Self;

    /// Number of qubits.
    fn n_qubits(&self) -> usize;

    /// X bit at `index`.
    fn get_xbit(&self, index: usize) -> bool;
    /// Z bit at `index`.
    fn get_zbit(&self, index: usize) -> bool;
    /// Loss bit at `index` (always `false` for non-lossy implementations).
    fn get_lbit(&self, index: usize) -> bool;

    /// Set the X bit at `index`.
    fn set_xbit(&mut self, index: usize, value: bool);
    /// Set the Z bit at `index`.
    fn set_zbit(&mut self, index: usize, value: bool);

    /// Number of non-identity slots (counts `X`, `Y`, `Z`, and — for
    /// lossy variants — `L`).
    fn weight(&self) -> usize;

    /// Number of slots marked as lost; always `0` for non-lossy variants.
    fn loss_weight(&self) -> usize;

    /// Recompute the cached hash. Call after batch mutations.
    fn rehash(&mut self);

    /// [`Pauli`] symbol at `index`.
    fn get(&self, index: usize) -> Pauli;

    /// Quick check: is the slot at `index` equal to `pauli`?
    fn is(&self, index: usize, pauli: Pauli) -> bool;

    /// Set the slot at `index` to `pauli`, in place. Returns `&mut self`
    /// for chaining.
    fn set(&mut self, index: usize, pauli: Pauli) -> &mut Self;

    /// Return a clone with the slot at `index` set to `pauli`.
    #[inline(always)]
    fn set_new(&self, index: usize, pauli: Pauli) -> Self {
        if index >= self.n_qubits() {
            panic!("Index out of bounds");
        }
        let mut new = self.clone();
        new.set(index, pauli);
        new
    }

    /// Check if this word anticommutes with a single-qubit Pauli at `addr0`,
    /// where `pauli = (xbit, zbit)`.
    #[inline(always)]
    fn anticommutes_at(&self, addr0: usize, pauli: (bool, bool)) -> bool {
        (self.get_xbit(addr0) & pauli.1) ^ (self.get_zbit(addr0) & pauli.0)
    }
}
