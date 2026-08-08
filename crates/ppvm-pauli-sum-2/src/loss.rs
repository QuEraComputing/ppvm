// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The neutral-atom **loss** channels on a [`Sum`]: [`LossChannel`],
//! [`CorrelatedLossChannel`] and [`ResetLossChannel`], ported from
//! `ppvm-pauli-sum/src/sum/noise.rs`.
//!
//! # One kernel, two word types, zero cost for the non-lossy one
//!
//! Architecture feature 11. Every loss branch below is guarded by
//! [`PauliBits::is_lost`], which is a **const `false`** for the ordinary
//! `PauliWord` and only genuine for `LossyPauliWord`
//! (`ppvm-lossy-pauli-word-2`). At monomorphization the non-lossy instantiation
//! therefore folds the whole loss ladder away — `loss_channel` collapses to the
//! pure in-place `c *= 1 − p`, `correlated_loss_channel` to `c *= 1 − 2p₁ − p₀`,
//! and the [`PauliBits::clear_lost`] / [`PauliBits::set_lost`] calls inside those
//! branches vanish with them (they are no-ops there anyway). That is exactly old's
//! arrangement, where `get_lbit` is an `#[inline(always)] false` and the same
//! kernels serve both words; it is also why the mutators live on `PauliBits` next
//! to the const-`false` reader rather than on the lossy word alone.
//!
//! [`ResetLossChannel`] is the one exception, and it is old's exception too: it
//! is implemented **only** for a `LossyPauliWord`-keyed sum, because "mark this
//! site lost" has no meaning for a word that cannot represent loss (a no-op
//! `set_lost` would make it emit a duplicate of the input key and silently double
//! coefficients).
//!
//! # Behaviour
//!
//! All three ride the fused in-place paths, so nothing truncates and nothing
//! reduces (behavioural contracts 1 and 2). The sharpest case is
//! `reset_loss_channel` on an already-lost site: old scales the coefficient to
//! **exactly zero and keeps the term**, which old's `tests/loss.rs::test_reset_channel`
//! pins with an exact-map comparison against `state.clone() *= 0.0`. Any implicit
//! zero-drop would fail it.
//!
//! # Not ported: `AsymmetricLossChannel`
//!
//! The trait exists in `ppvm-traits-2` (and existed in `ppvm-traits`), but the old
//! `PauliSum` **never implemented it** — `grep AsymmetricLossChannel
//! crates/ppvm-pauli-sum/src` is empty. Adding one here would be a new surface,
//! not a port, and the design's rule is that a behaviour-preserving port neither
//! widens nor narrows the surface.

use std::hash::BuildHasher;

use ppvm_lossy_pauli_word_2::LossyPauliWord;
use ppvm_pauli_word_2::{HashFinalize, PauliStorage};
use ppvm_traits_2::{
    Accumulate, Coefficient, CorrelatedLossChannel, Indexable, LossChannel, PauliBits,
    ResetLossChannel, Retain, Word,
};

use crate::store::{BranchInPlace, RotateInPlace};
use crate::sum::Sum;

/// A copy of `k` with the loss flag at `q` cleared — the branch key of a loss
/// channel's "the atom comes back as identity" outcome (old's
/// `k.set(addr, Pauli::I)` on a lost site).
#[inline]
fn with_loss_cleared<W: PauliBits + Clone>(k: &W, q: usize) -> W {
    k.loss_cleared(q)
}

/// The single-qubit loss channel: with probability `p` the atom is lost.
///
/// * a site **already** lost branches to its `I` image at `p·c`, and the lost
///   term itself is left **unscaled** (old does not scale `v` in that arm — the
///   `L` term keeps its coefficient and the `I` branch is added beside it, which
///   old's `tests/loss.rs::test_loss_channel` pins as `{L: 1.0, I: 0.2}`);
/// * every other site is a pure in-place `c *= 1 − p`, no branch.
///
/// Ported from `ppvm-pauli-sum/src/sum/noise.rs`'s `impl LossChannel`, which is
/// likewise generic over the word type.
impl<S, P, W, C> LossChannel<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + RotateInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One,
    P: crate::policy::Policy<W, C>,
{
    fn loss_channel(&mut self, qubit: usize, p: C) {
        let survive = C::one() - p.clone();
        self.rotate_in_place(move |k: &W, c: &mut C| {
            if k.is_lost(qubit) {
                let branch = c.clone() * p.clone();
                // NOTE: the survivor is deliberately *not* scaled — old's `L` arm
                // returns the branch and leaves `v` alone.
                Some((with_loss_cleared(k, qubit), branch))
            } else {
                *c *= survive.clone();
                None
            }
        });
    }
}

/// The correlated two-qubit loss channel.
///
/// `p[0]` is the probability of losing **both** atoms when both are in the qubit
/// subspace, `p[1]` of losing **either one** in that case, and `p[2]` of losing
/// one when the other is **already** lost. The four arms are old's verbatim
/// (`ppvm-pauli-sum/src/sum/noise.rs`), including two inconsistencies flagged for
/// the Lean oracle rather than "fixed" here (suspected old bug 4):
///
/// 1. in the one-already-lost arms the survivor is scaled by `1 − p[2]` while the
///    emitted branch is weighted `p[1]`, so the branch weight and the survivor's
///    complement do not pair and the channel is not trace-preserving unless
///    `p[1] == p[2]`;
/// 2. the both-lost arm emits three branches but never scales the survivor at
///    all, unlike every other arm.
///
/// No test in old covers correlated loss with distinct `p[0]`/`p[1]`/`p[2]`, so
/// both are unpinned there; reproducing them keeps the golden master exact and
/// leaves the adjudication (a CPTP / trace-preservation theorem) to Lean.
///
/// This is the crate's only **multi-branch** channel, so it is the one driver of
/// [`Sum::branch_in_place`] and hence of the size-directed merge (architecture
/// feature 3).
impl<S, P, W, C> CorrelatedLossChannel<C> for Sum<S, P>
where
    S: Accumulate<Key = W, Coeff = C> + BranchInPlace<W, C> + Retain<W, C>,
    W: Word + Indexable + PauliBits,
    C: Coefficient + num::One,
    P: crate::policy::Policy<W, C>,
{
    fn correlated_loss_channel(&mut self, qubit0: usize, qubit1: usize, p: [C; 3]) {
        let [p0, p1, p2] = p;
        // `1 − 2·p1 − p0` for the both-present arm, with the doubling as `x + x`
        // (the `Coefficient` ring carries no bare-`f64` scaling).
        let both_present = C::one() - (p1.clone() + p1.clone()) - p0.clone();
        let one_lost_survivor = C::one() - p2.clone();
        self.branch_in_place(move |k: &W, c: &mut C, sink: &mut Vec<(W, C)>| {
            match (k.is_lost(qubit0), k.is_lost(qubit1)) {
                (true, true) => {
                    // Both lost: one branch per single recovery at `p[2]`, one
                    // for the double recovery at `p[0]`. The survivor is not
                    // scaled (old).
                    let single = c.clone() * p2.clone();
                    sink.push((with_loss_cleared(k, qubit0), single.clone()));
                    sink.push((with_loss_cleared(k, qubit1), single));
                    let mut both = k.clone();
                    both.clear_lost(qubit0);
                    both.clear_lost(qubit1);
                    sink.push((both, c.clone() * p0.clone()));
                }
                (false, true) => {
                    sink.push((with_loss_cleared(k, qubit1), c.clone() * p1.clone()));
                    *c *= one_lost_survivor.clone();
                }
                (true, false) => {
                    sink.push((with_loss_cleared(k, qubit0), c.clone() * p1.clone()));
                    *c *= one_lost_survivor.clone();
                }
                (false, false) => *c *= both_present.clone(),
            }
        });
    }
}

/// The reset-loss channel: model a re-cooling event that brings a lost atom back,
/// **only** for a `LossyPauliWord`-keyed sum (old restricts it the same way).
///
/// * `I` and `Z` branch to their `L` image at the *same* coefficient (the term
///   itself stays — old returns `Some` from `map_insert`, which mutates nothing on
///   the survivor);
/// * `X` and `Y` are untouched;
/// * an already-lost site has its coefficient scaled to **exactly zero** and the
///   term **stays in the support** — the sharpest available test of the
///   no-implicit-reduce contract, pinned by old's exact-map
///   `tests/loss.rs::test_reset_channel` against `state.clone() *= 0.0`.
impl<S, P, A, H, C> ResetLossChannel for Sum<S, P>
where
    S: Accumulate<Key = LossyPauliWord<A, H>, Coeff = C>
        + RotateInPlace<LossyPauliWord<A, H>, C>
        + Retain<LossyPauliWord<A, H>, C>,
    A: PauliStorage,
    H: BuildHasher + Default + HashFinalize,
    C: Coefficient,
    P: crate::policy::Policy<LossyPauliWord<A, H>, C>,
{
    fn reset_loss_channel(&mut self, qubit: usize) {
        self.rotate_in_place(move |k: &LossyPauliWord<A, H>, c: &mut C| {
            if k.is_lost(qubit) {
                // Zeroed **in place**: the key stays in the support at 0.0.
                *c *= C::zero();
                return None;
            }
            // I / Z branch to the lost image; X / Y do nothing.
            if k.x_bit(qubit) {
                return None;
            }
            Some((k.with_lost(qubit), c.clone()))
        });
    }
}
