// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! [`TableauLike`](crate::tableau_like::TableauLike) trait: shared interface for stabilizer-tableau backends.
//!
//! Any type implementing [`TableauLike`](crate::tableau_like::TableauLike) gets default implementations of
//! single- and two-qubit Pauli noise channels (Depolarizing, PauliError,
//! Depolarizing2, TwoQubitPauliError). New tableau backends only need to
//! provide [`TableauLike::rng_mut`](crate::tableau_like::TableauLike::rng_mut) and (optionally) override
//! [`TableauLike::is_qubit_lost`](crate::tableau_like::TableauLike::is_qubit_lost)
//! (which reports lost *or* leaked qubits; two-qubit channels marginalize them out).

use num::Zero;
use rand::{Rng, RngExt};

use ppvm_traits::traits::{Clifford, Coefficient};

// `Zero` is used for `Self::Coeff::zero()` inside default method bodies; the
// bound itself is redundant on `Self::Coeff` because `Coefficient: num::Zero`.

/// `true` iff `p` lies in the closed unit interval `[0, 1]`.
///
/// Spelled out with two explicit comparisons rather than
/// `(0.0..=1.0).contains(p)` because the coefficient types here only implement
/// `PartialOrd<f64>`, not the `PartialOrd<Self>` that `RangeInclusive::contains`
/// requires. Keeping the two comparisons in separate statements also avoids the
/// `clippy::manual_range_contains` rewrite, which would not compile here.
#[inline]
fn is_probability<C: PartialOrd<f64>>(p: &C) -> bool {
    if *p < 0.0 {
        return false;
    }
    *p <= 1.0
}

/// Marginal X/Y/Z probabilities on one qubit of a two-qubit Pauli channel.
///
/// `p` uses the [`TableauLike::two_qubit_pauli_error_impl`] layout, where
/// `p[4a + b - 1]` is the probability of `P_a ⊗ P_b` (0 = I, 1 = X, 2 = Y,
/// 3 = Z). Sums out the other qubit's Pauli: `q_s = Σ_o p(s, o)` when keeping
/// the first qubit, `q_s = Σ_o p(o, s)` when keeping the second.
fn marginal_pauli_probs<C: Coefficient>(p: &[C; 15], keep_first: bool) -> [C; 3] {
    core::array::from_fn(|i| {
        let s = i + 1;
        (0..4)
            .map(|o| if keep_first { 4 * s + o } else { 4 * o + s })
            .fold(C::zero(), |acc, k| acc + p[k - 1].clone())
    })
}

/// A stabilizer-tableau-like backend that supports Clifford gates and an RNG.
///
/// Implementing this trait grants default implementations of the Pauli noise
/// channels via [`TableauLike::depolarize_impl`] and friends. The associated
/// `Rng` type lets each backend choose its own RNG; nothing in this trait
/// depends on `SmallRng`.
pub trait TableauLike: Clifford {
    /// Coefficient type used for probabilities.
    type Coeff: Coefficient + PartialOrd<f64>;

    /// RNG type backing the stochastic channels.
    type Rng: Rng + RngExt;

    /// Mutable access to the backend's RNG.
    fn rng_mut(&mut self) -> &mut Self::Rng;

    /// Whether the qubit at `addr` is inactive (lost or leaked), so the
    /// two-qubit Pauli channels act only on the surviving qubit. Default: never inactive.
    #[inline]
    fn is_qubit_lost(&self, _addr: usize) -> bool {
        false
    }

    /// Single-qubit depolarizing channel.
    ///
    /// RNG is consumed unconditionally; the selected Clifford gate is expected
    /// to no-op on lost qubits. This preserves seeded RNG sequences across
    /// loss events.
    #[inline]
    fn depolarize_impl(&mut self, addr0: usize, p: Self::Coeff) {
        debug_assert!(is_probability(&p));
        let r = self.rng_mut().random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p / 3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }

    /// Single-qubit Pauli-error channel (X, Y, Z with given probabilities).
    ///
    /// RNG is consumed unconditionally; the selected Clifford gate is expected
    /// to no-op on lost qubits. This preserves seeded RNG sequences across
    /// loss events.
    #[inline]
    fn pauli_error_impl(&mut self, addr0: usize, p: [Self::Coeff; 3]) {
        debug_assert!(p.iter().all(is_probability));
        let r = self.rng_mut().random::<f64>();
        let mut cumulative = Self::Coeff::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += p_.clone();
            if cumulative > r {
                match i {
                    0 => self.x(addr0),
                    1 => self.y(addr0),
                    _ => self.z(addr0),
                }
                return;
            }
        }
    }

    /// Two-qubit Pauli-error channel (15 non-identity Pauli combinations).
    ///
    /// If exactly one qubit is lost or leaked, the survivor gets the marginal
    /// single-qubit channel (the dead qubit's Pauli is summed out). If both
    /// are, nothing is applied.
    #[inline]
    fn two_qubit_pauli_error_impl(&mut self, addr0: usize, addr1: usize, p: [Self::Coeff; 15]) {
        debug_assert!(p.iter().all(is_probability));

        let is_lost0 = self.is_qubit_lost(addr0);
        let is_lost1 = self.is_qubit_lost(addr1);
        if is_lost0 && is_lost1 {
            return;
        }

        if is_lost0 || is_lost1 {
            let addr = if is_lost0 { addr1 } else { addr0 };
            let marginal_probs = marginal_pauli_probs(&p, !is_lost0);
            self.pauli_error_impl(addr, marginal_probs);
            return;
        }

        let r = self.rng_mut().random::<f64>();
        let sum = Self::Coeff::zero();
        let idx = p
            .iter()
            .scan(sum, |acc, p_| {
                *acc += p_.clone();
                Some(acc.clone())
            })
            .position(|cum_prob| cum_prob > r);

        if let Some(i) = idx {
            // Two-qubit Pauli pairs indexed by `i + 1` (so II at index 0 is
            // skipped). Encoding: 0 = I, 1 = X, 2 = Y, 3 = Z; the first entry
            // acts on addr0, the second on addr1.
            //
            // `rustfmt::skip` keeps the rows grouped by the first Pauli (a
            // readable 4-wide grid); without it rustfmt repacks the tuples to
            // fill the line width and the grouping is lost.
            #[rustfmt::skip]
            const PAULI_PAIRS: [(u8, u8); 16] = [
                (0,0),(0,1),(0,2),(0,3),
                (1,0),(1,1),(1,2),(1,3),
                (2,0),(2,1),(2,2),(2,3),
                (3,0),(3,1),(3,2),(3,3),
            ];
            let cartesian_index = PAULI_PAIRS[i + 1]; // skip II entry

            match cartesian_index.0 {
                0 => {}
                1 => self.x(addr0),
                2 => self.y(addr0),
                _ => self.z(addr0),
            }

            match cartesian_index.1 {
                0 => {}
                1 => self.x(addr1),
                2 => self.y(addr1),
                _ => self.z(addr1),
            }
        }
    }

    /// Two-qubit depolarizing channel: spreads `p` over the 15 non-identity
    /// two-qubit Pauli errors. A lone surviving qubit is depolarized with `4p/5`.
    #[inline]
    fn depolarize2_impl(&mut self, addr0: usize, addr1: usize, p: Self::Coeff) {
        debug_assert!(is_probability(&p));
        let p_arr: [Self::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        self.two_qubit_pauli_error_impl(addr0, addr1, p_arr);
    }
}
