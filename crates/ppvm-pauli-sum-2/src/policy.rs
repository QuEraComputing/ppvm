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
//! Ported from `ppvm-pauli-sum/src/strategy.rs`, including the loss-weight
//! policy [`MaxLossWeight`].

use ppvm_traits_2::{Coefficient, Indexable, PauliBits, Retain, Word};

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

/// The ceiling [`NoPolicy`]'s exponential capacity hint is clamped to.
///
/// `16_384` terms — past the ~10⁴-term support an untruncated deep random
/// circuit reaches, so that workload's map never resizes mid-circuit, while
/// costing well under a megabyte instead of the tens of gigabytes the unclamped
/// `4ⁿ/2` asks for at the widths the engine actually runs at. Old's number is
/// reproduced exactly for `n_sites ≤ 7`.
const NO_POLICY_CAPACITY_CEILING: usize = 1 << 14;

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
    /// Old's `NoStrategy::capacity(n) = 1 << (2n − 1)` — "in exact simulation,
    /// guess `4ⁿ/2` paths" (`ppvm-traits/src/traits/strategy.rs`) — **clamped**
    /// to [`NO_POLICY_CAPACITY_CEILING`].
    ///
    /// Reproducing the formula verbatim is not viable: it is exact for `n ≤ 7`
    /// but at `n = 16` it asks for `2³¹` entries (~80 GB across the primary and
    /// aux maps) and at `n ≥ 32` the shift itself overflows. Old gets away with
    /// it only because every one of its own workloads either runs `NoStrategy`
    /// at a handful of qubits or overrides the hint through the builder's
    /// `.capacity(..)`. The clamp keeps old's number wherever old's number is
    /// allocatable and keeps the *purpose* of the hint (an untruncated
    /// fan-out-driven run must not rehash mid-circuit) everywhere else; a caller
    /// who wants a specific size passes it to
    /// [`Sum::with_capacity`](crate::Sum::with_capacity), the port of old's
    /// builder override. Deliberate, documented divergence.
    #[inline]
    fn capacity(&self, n_sites: usize) -> usize {
        if n_sites == 0 {
            // Old underflows `2*n - 1` here; an empty word admits one term.
            return 1;
        }
        let shift = 2 * n_sites - 1;
        if shift >= usize::BITS as usize {
            return NO_POLICY_CAPACITY_CEILING;
        }
        (1usize << shift).min(NO_POLICY_CAPACITY_CEILING)
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

    #[inline(always)]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        // `usize::MAX` is the "disabled" sentinel: skip the full support scan.
        //
        // Lean: `retain_le_top_eq_self` / `retain_of_all_true` in
        // `lean/PPVM/Algebra/GradedMap.lean` — a cap at the top of the weight
        // order accepts every key, and retain-all is the identity *pointwise*
        // (zero-coefficient terms included), so skipping the walk entirely is
        // observationally exact, not an approximation.
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
/// Design: §"`Policy`" (the `CoefficientThreshold` impl, named-field form).
///
/// The [`Default`] threshold is **`1e-12`**, not `0.0`: that is old's
/// `impl Default for CoefficientThreshold { fn default() -> Self { Self(1e-12) } }`
/// (`ppvm-pauli-sum/src/strategy.rs`), and it is user-facing — a `Sum::new(n)`
/// carrying the default policy drops everything below `1e-12` on `truncate()`.
/// The design sketch sets it to `0.0`; the sketch is what has to change, since a
/// silently loosened default would keep terms old dropped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoefficientThreshold {
    /// Terms with `magnitude() < threshold` are dropped.
    pub threshold: f64,
}

impl Default for CoefficientThreshold {
    /// `1e-12` — old's default (`ppvm-pauli-sum/src/strategy.rs`).
    #[inline]
    fn default() -> Self {
        Self { threshold: 1e-12 }
    }
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

    #[inline(always)]
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
/// `ppvm-pauli-sum::strategy::CombinedStrategy`, including its **two sequential
/// `retain` passes** (old documents this as a structural difference from
/// PauliPropagation.jl's single combined walk).
///
/// Lean: `retain_seq_eq_retain_and` in `lean/PPVM/Algebra/GradedMap.lean` — the
/// sequential form computes the **conjunction** of the two keep-rules — and
/// `retain_comm`, which says the surviving key set is a property of the policy
/// pair, not an artifact of the pass order. Together they license a future
/// backend to fuse or reorder the two walks.
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

    #[inline(always)]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        self.0.truncate(map);
        self.1.truncate(map);
    }
}

/// Drop terms whose **loss weight** (number of lost sites) exceeds the bound.
///
/// Only meaningful for a lossy key: [`PauliBits::loss_weight`] is a const `0` for
/// the ordinary `PauliWord`, so this policy keeps everything there. `usize::MAX`
/// disables the pass entirely (the same sentinel [`MaxPauliWeight`] uses).
///
/// The [`Default`] is **`10`**, not `usize::MAX` — old's
/// `impl Default for MaxLossWeight { fn default() -> Self { Self(10) } }`
/// (`ppvm-pauli-sum/src/strategy.rs`), and it is user-facing, so it is reproduced
/// exactly (behavioural contract 7).
///
/// Design: §"`Policy`". Ported from `ppvm-pauli-sum::strategy::MaxLossWeight`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MaxLossWeight(pub usize);

impl MaxLossWeight {
    /// Maximum loss weight retained.
    #[inline]
    pub fn max_loss_weight(&self) -> usize {
        self.0
    }
}

impl Default for MaxLossWeight {
    /// `10` — old's default.
    #[inline]
    fn default() -> Self {
        Self(10)
    }
}

impl<W, C> Policy<W, C> for MaxLossWeight
where
    W: Word + Indexable + PauliBits,
    C: Coefficient,
{
    /// Old's conservative `n_sites * 10` (clearing a map scales with capacity).
    #[inline]
    fn capacity(&self, n_sites: usize) -> usize {
        n_sites * 10
    }

    #[inline]
    fn truncate<M>(&self, map: &mut M)
    where
        M: Retain<W, C>,
    {
        // The `usize::MAX` disable sentinel, as in `MaxPauliWeight`: skip the
        // whole retain pass rather than walking every bucket for a predicate that
        // is always true (architecture feature 7).
        if self.0 == usize::MAX {
            return;
        }
        map.retain(|w, _| w.loss_weight() <= self.0);
    }
}
