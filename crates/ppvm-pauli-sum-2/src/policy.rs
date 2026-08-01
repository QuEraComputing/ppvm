// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Truncation [`Policy`] — the redesign's renaming of the old `Strategy` — and
//! its concrete implementations [`NoPolicy`], [`MaxPauliWeight`],
//! [`CoefficientThreshold`], and [`CombinedPolicy`].
//!
//! `Policy` retains `Strategy`'s two responsibilities (capacity hint + truncate)
//! but **drops the `Copy` bound** and consumes the non-algebraic
//! [`Retain`](ppvm_traits_2::Retain) capability rather than the map algebra:
//! dropping supported terms breaks module exactness, so truncation lives outside
//! the graded algebra (Design: §"The map is a graded algebra over `C[K]`" and
//! §"`truncate` bounds on `Retain`").
//!
//! Ported from `ppvm-pauli-sum/src/strategy.rs`. Loss-weight truncation
//! (`MaxLossWeight`) is deferred to the lossy component.

use ppvm_traits_2::{Coefficient, Indexable, Retain, Word};

/// A truncation policy: predict an initial map capacity, and drop terms outside
/// the policy's support.
///
/// The key bound is `Word + Indexable`: `MaxPauliWeight` needs
/// [`Word::weight`], and truncation only ever applies to a keyed sum.
/// `truncate` is generic over the storage `M`, bounded on
/// [`Retain`](ppvm_traits_2::Retain), so a policy never names the concrete
/// backend.
///
/// Design: §"`Policy`". The truncation error a policy incurs is bounded in
/// `lean/PPVM/Algebra/Truncation.lean` (`l1_bound`), and the `>=`-vs-`>` keep-rule
/// boundary against the tableau path is `cutoff_mismatch`.
pub trait Policy<W, C>: Default + Clone
where
    W: Word + Indexable,
    C: Coefficient,
{
    /// Predicted initial capacity for a sum of `n_sites`-wide keys.
    fn capacity(&self, n_sites: usize) -> usize;

    /// Drop the terms this policy does not retain.
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>;
}

/// No truncation: keep every term. The default policy.
///
/// Design: §"`Policy`" (`NoStrategy` → `NoPolicy`).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct NoPolicy;

impl<W, C> Policy<W, C> for NoPolicy
where
    W: Word + Indexable,
    C: Coefficient,
{
    #[inline]
    fn capacity(&self, _n_sites: usize) -> usize {
        0
    }

    #[inline]
    fn truncate<M>(&self, _map: &mut M)
    where
        M: Retain<W, C>,
    {
    }
}

/// Drop terms whose Pauli weight (number of non-identity sites) exceeds the
/// bound. `usize::MAX` (the [`Default`]) disables truncation.
///
/// Design: §"`Policy`". Uses [`Word::weight`]; ported from
/// `ppvm-pauli-sum::strategy::MaxPauliWeight`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxPauliWeight(pub usize);

impl MaxPauliWeight {
    /// Maximum Pauli weight retained.
    #[inline]
    pub fn max_weight(&self) -> usize {
        self.0
    }
}

impl Default for MaxPauliWeight {
    #[inline]
    fn default() -> Self {
        Self(usize::MAX)
    }
}

impl<W, C> Policy<W, C> for MaxPauliWeight
where
    W: Word + Indexable,
    C: Coefficient,
{
    #[inline]
    fn capacity(&self, n_sites: usize) -> usize {
        // Conservative: capacity has a direct performance impact (clearing scales
        // with it), so avoid the binomial blow-up. Ported from the old strategy.
        n_sites * 10
    }

    #[inline]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        // `usize::MAX` is the "disabled" sentinel: skip the full support scan.
        if self.0 == usize::MAX {
            return;
        }
        map.retain(|w, _| w.weight() <= self.0);
    }
}

/// Drop terms whose coefficient magnitude falls **below** the threshold.
///
/// The keep-rule is `magnitude() >= threshold`, so a term whose magnitude is
/// *exactly* the threshold is kept. The stabilizer-tableau path keeps on a
/// strict `>`, so the two disagree at `|c| == threshold` — a boundary mismatch
/// machine-checked in `lean/PPVM/Algebra/Truncation.lean` (`cutoff_mismatch`).
///
/// Design: §"`Policy`" (the `CoefficientThreshold` impl, named-field form). The
/// [`Default`] threshold is `0.0`, which keeps every term (magnitude is always
/// `>= 0`).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct CoefficientThreshold {
    /// Terms with `magnitude() < threshold` are dropped.
    pub threshold: f64,
}

impl<W, C> Policy<W, C> for CoefficientThreshold
where
    W: Word + Indexable,
    C: Coefficient,
{
    #[inline]
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites * 10
    }

    #[inline]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        map.retain(|_, coeff| coeff.magnitude() >= self.threshold);
    }
}

/// Two policies run in sequence: `P1` then `P2`. The capacity hint is the
/// smaller of the two.
///
/// Design: §"`Policy`" (`CombinedStrategy` → `CombinedPolicy`). Ported from
/// `ppvm-pauli-sum::strategy::CombinedStrategy`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct CombinedPolicy<P1, P2>(pub P1, pub P2);

impl<W, C, P1, P2> Policy<W, C> for CombinedPolicy<P1, P2>
where
    W: Word + Indexable,
    C: Coefficient,
    P1: Policy<W, C>,
    P2: Policy<W, C>,
{
    #[inline]
    fn capacity(&self, n_sites: usize) -> usize {
        self.0.capacity(n_sites).min(self.1.capacity(n_sites))
    }

    #[inline]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        self.0.truncate(map);
        self.1.truncate(map);
    }
}
