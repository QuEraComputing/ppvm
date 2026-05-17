use crate::char::Pauli;
use crate::traits::PauliStorage;
use std::fmt::Display;
use std::hash::Hash;

pub trait PauliIter {
    fn iter(&self) -> impl Iterator<Item = Pauli>;
}

/// Word-level Pauli operations. Types implementing this trait automatically
/// gain [`crate::traits::Clifford`] and [`crate::traits::CliffordExtensions`]
/// behavior via blanket impls in [`crate::traits::clifford`], using
/// `get_lbit` / `get_xbit` / `set_xbit` / `set_zbit` / `rehash` to transform
/// Pauli words.
///
/// # Word-level vs. state-level gate semantics
///
/// The blanket `Clifford` impl applies gates at the *bit* level of a single
/// Pauli word. X / Y / Z are bit-level no-ops on a Pauli word — they affect
/// phase, not the X/Z bits — so `word.x(i)`, `word.y(i)`, `word.z(i)` are
/// deliberately silent. Phase is tracked separately by
/// [`crate::phase::PhasedPauliWord`], which implements `Clifford` manually.
///
/// If you need a word representation whose gate behavior is *not* pure
/// bit manipulation (phase tracking, fused multi-qubit updates, alternative
/// loss semantics), do **not** implement `PauliWordTrait` on it — define a
/// specialized `Clifford` impl instead, the way `PhasedPauliWord` does.
/// Implementing both `PauliWordTrait` and a custom `impl Clifford for ...`
/// for the same type will not compile: the blanket impl overlaps with your
/// custom impl (coherence error), not silently shadows it.
pub trait PauliWordTrait: Clone + Hash + Eq + PauliIter + From<String> + Display {
    fn new(nqubits: usize) -> Self;

    fn n_qubits(&self) -> usize;

    // getter methods
    fn get_xbit(&self, index: usize) -> bool;
    fn get_zbit(&self, index: usize) -> bool;
    fn get_lbit(&self, index: usize) -> bool;

    // setter methods
    fn set_xbit(&mut self, index: usize, value: bool);
    fn set_zbit(&mut self, index: usize, value: bool);

    fn weight(&self) -> usize;

    fn loss_weight(&self) -> usize;

    fn rehash(&mut self);

    fn get(&self, index: usize) -> Pauli;

    fn get_multiple<const Q: usize>(&self, indices: [usize; Q]) -> Self {
        let mut result = Self::new(Q);
        for (i, &idx) in indices.iter().enumerate() {
            result.set(i, self.get(idx));
        }
        result
    }

    fn set_multiple<const Q: usize, B: PauliStorage>(&mut self, indices: [usize; Q], values: &Self);

    fn get_slice(&self, slice: std::ops::Range<usize>) -> Self;

    fn is(&self, index: usize, pauli: Pauli) -> bool;

    fn set(&mut self, index: usize, pauli: Pauli) -> &mut Self;

    #[inline(always)]
    fn set_new(&self, index: usize, pauli: Pauli) -> Self {
        if index >= self.n_qubits() {
            panic!("Index out of bounds");
        }
        let mut new = self.clone();
        new.set(index, pauli);
        new
    }

    #[inline(always)]
    fn set_new_2(&self, index_0: usize, pauli_0: Pauli, index_1: usize, pauli_1: Pauli) -> Self {
        let mut new = self.clone();
        new.set(index_0, pauli_0);
        new.set(index_1, pauli_1);
        new
    }

    /// Check if this word anticommutes with a single-qubit Pauli at `addr0`,
    /// where `pauli = (xbit, zbit)`.
    #[inline(always)]
    fn anticommutes_at(&self, addr0: usize, pauli: (bool, bool)) -> bool {
        (self.get_xbit(addr0) & pauli.1) ^ (self.get_zbit(addr0) & pauli.0)
    }
}
