// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`RekeyProducer`] — the bijective re-key term producer a Clifford gate feeds
//! into [`Sum::apply`](crate::Sum::apply).
//!
//! A Clifford conjugation is a pushforward along a Pauli bijection `w ↦ φ(w)`:
//! exactly one produced term per input, and — because `φ` is injective — no
//! collisions. [`Sum::apply`](crate::Sum::apply) runs neither `reduce` nor the
//! policy's truncation (both are caller-driven), so a produced term whose
//! coefficient is exactly zero survives, as in old.
//! Rotation and the diagonal noise channel take dedicated in-place fast paths
//! ([`RotateInPlace`](crate::store::RotateInPlace) /
//! [`ScaleByKey`](crate::store::ScaleByKey)) rather than a batch producer, so this
//! module carries only the bijective re-key; the L4 multiply producer (an outer
//! product) is deferred to a later component.
//!
//! Design: §"Every gate is a producer feeding `accumulate`" (the `RekeyProducer`
//! sketch). The Clifford re-key is the symplectic conjugation bijection: each
//! generator's `Sp(2n, 2)` bit map is proven an involution, hence a bijection, in
//! `lean/PPVM/Pauli/Symplectic.lean` (`hAct_involutive`/`sAct_involutive`/
//! `cnotAct_involutive`/`czAct_involutive`, `*_bijective`), and the corresponding
//! conjugation homs are injective in `lean/PPVM/Pauli/Conjugation.lean`
//! (`conjH_injective`/`conjS_injective`/`conjCNOT_injective`/`conjCZ_injective`) —
//! which is exactly the "no collision" guarantee this producer relies on. (This is
//! the symplectic conjugation, **not** the amplitude-vector relabel
//! `xorRelabel_bijective` of `lean/PPVM/Instantiations/Bitstring.lean`, which
//! validates the GeneralizedTableau non-Clifford branch.)

use ppvm_traits_2::{TermProducer, TermSink};

/// A bijective re-key: one produced term `f(k, c)` per input `(k, c)`.
///
/// The closure `f` is monomorphized and inlined into the accumulate loop; when
/// it captures only `Copy` state (a gate's qubit indices), `RekeyProducer` is
/// effectively a ZST. Never `dyn`.
///
/// Design: §"Every gate is a producer feeding `accumulate`".
pub struct RekeyProducer<F> {
    f: F,
}

impl<F> RekeyProducer<F> {
    /// Wrap a re-key closure `f: Fn(&K, &C) -> (K, C)`.
    #[inline]
    pub fn new(f: F) -> Self {
        Self { f }
    }
}

impl<K, C, F> TermProducer<K, C> for RekeyProducer<F>
where
    // `Send + Sync` on the closure is what makes `RekeyProducer<F>` satisfy
    // `TermProducer`'s `Send + Sync` supertrait — architecture feature 12: a
    // producer is `&self`-only, so a concurrent storage backend can split the
    // produce walk without widening a single signature.
    F: Fn(&K, &C) -> (K, C) + Send + Sync,
{
    #[inline(always)]
    fn produce<S: TermSink<K, C>>(&self, key: &K, coeff: &C, sink: &mut S) {
        let (k, c) = (self.f)(key, coeff);
        sink.push(k, c);
    }
}
