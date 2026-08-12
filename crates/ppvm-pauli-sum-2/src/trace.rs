// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Expectation extraction: the [`Trace`] pairing of a [`Sum`] against a
//! **pattern** — the set of Pauli words a measurement contracts onto.
//!
//! This is the port of the old crate's `Trace` surface
//! (`ppvm-pauli-word/src/pattern/`, `ppvm-pauli-sum/src/sum/trace.rs`), which is
//! how every expectation value in the old test-suite is read out: propagate an
//! observable backwards through the circuit, then contract the resulting sum
//! against the state's stabilizer pattern. For the all-zero state that pattern is
//! `Z?*` — every site is `I` or `Z` — because `⟨0…0|P|0…0⟩` is `1` for such a `P`
//! and `0` for any `P` carrying an `X` or `Y`. So
//!
//! ```text
//! ⟨0…0| O |0…0⟩ = Σ_{k ∈ supp(O), k matches Z?*} O[k]
//! ```
//!
//! which is exactly [`Trace::trace`] here. It backs `tests/ghz.rs`,
//! `tests/cnot.rs` and `tests/noise.rs::test_depolarizing_error` in the old
//! crate (e.g. three independent `depolarize1(i, pᵢ)` on `ZZZ` contract to
//! `Πᵢ (1 − 4pᵢ/3)`).
//!
//! Design: §"Compatibility with current names" retains
//! [`Trace<'a, RHS>`](ppvm_traits_2::Trace) in `ppvm-traits-2` precisely so its
//! implementers can land with the pattern port; this is that landing. The
//! same-type sum-against-sum pairing is [`Sum::overlap`], whose bilinearity is
//! machine-checked in `lean/PPVM/Algebra/GradedMap.lean`.
//!
use num::Zero;
use ppvm_traits_2::{Accumulate, Coefficient, Indexable, PauliBits, Trace, Word};

use crate::pattern::{PatternSite, PauliPattern};
use crate::policy::Policy;
use crate::sum::Sum;

/// A pattern traces a *word* to a boolean membership — old's
/// `impl Trace<'a, W> for PauliPattern` (`Output = bool`), which is the shape the
/// sum-side trace below consumes.
impl<'a, W> Trace<'a, W> for PauliPattern
where
    W: Word + PauliBits + 'a,
    W::Site: PatternSite,
{
    type Output = bool;

    #[inline]
    fn trace(&'a self, value: &'a W) -> bool {
        self.matches(value)
    }
}

/// A sum traces against a pattern to the **sum of the coefficients of the
/// matching keys**.
///
/// Old's `impl Trace<'a, Rhs> for PauliSum`, which delegates to the map's
/// `fold(zero, |acc, (k, v)| { value.trace(k).then(|| acc += v); acc })`
/// (`ppvm-traits/src/map/hashmap.rs`). Zero-coefficient terms contribute zero, so
/// the contract is insensitive to whether the support carries them; the
/// *summation order* is the backend's iteration order on both sides, so a
/// floating-point contraction agrees only up to reassociation error.
///
/// This is stated for the concrete [`PauliPattern`] rather than for every
/// `R: Trace<'a, S::Key, Output = bool>`: the engine's read side yields the key
/// with the *scan's* lifetime, not the `'a` of the sum, so it cannot satisfy the
/// `&'a RHS` the generic form would demand. Old could be generic because its map
/// iterator borrows for `'a`. `PauliPattern` is the only implementer old ever
/// had.
///
/// # Only matching coefficients are accumulated
///
/// The fold runs over
/// [`Support::for_each_ref`](ppvm_traits_2::Support::for_each_ref) — the
/// borrowing scan — and calls [`Coefficient::add_assign_ref`] only *after* the
/// pattern accepts its key. Its default is old's
/// `acc += v.clone()` term for term; heap-backed coefficient rings can avoid a
/// whole temporary clone while preserving that value. Reading through
/// [`Support::iter`](ppvm_traits_2::Support::iter) instead would clone every
/// coefficient in the support before the filter looked at it: invisible for
/// `f64`, but on a symbolic coefficient (an owned monomial table per term) that
/// measured 7×–33× slower than old on `sym.random.circuit`, growing with
/// coefficient size — the pattern rejects 255 of 65534 keys there, so the clones
/// were ~99.6% waste.
impl<'a, S, P> Trace<'a, PauliPattern> for Sum<S, P>
where
    S: Accumulate,
    P: Policy<S::Key, S::Coeff>,
    S::Key: Word + PauliBits + Indexable,
    <S::Key as Word>::Site: PatternSite,
{
    type Output = S::Coeff;

    fn trace(&'a self, value: &'a PauliPattern) -> S::Coeff {
        // `zero() + …` in backend order, i.e. old's fold — *not* `Iterator::sum`
        // over a filtered iterator, which would have to own each pair.
        let mut acc = S::Coeff::zero();
        self.for_each_ref(|k, c| {
            if value.matches_bits(k) {
                acc.add_assign_ref(c);
            }
        });
        acc
    }
}
