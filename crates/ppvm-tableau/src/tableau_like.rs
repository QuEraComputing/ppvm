// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! [`TableauLike`](crate::tableau_like::TableauLike) trait: shared interface for stabilizer-tableau backends.
//!
//! Any type implementing [`TableauLike`](crate::tableau_like::TableauLike) gets default implementations of
//! single- and two-qubit Pauli noise channels (Depolarizing, PauliError,
//! Depolarizing2, TwoQubitPauliError). New tableau backends only need to
//! provide [`TableauLike::rng_mut`](crate::tableau_like::TableauLike::rng_mut) and (optionally) override
//! [`TableauLike::is_qubit_lost`](crate::tableau_like::TableauLike::is_qubit_lost).

use num::Zero;
use rand::{Rng, RngExt};

use ppvm_runtime::traits::{Clifford, Coefficient};

// `Zero` is used for `Self::Coeff::zero()` inside default method bodies; the
// bound itself is redundant on `Self::Coeff` because `Coefficient: num::Zero`.

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

    /// Whether the qubit at `addr` is lost. Default: never lost.
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
        #[allow(clippy::manual_range_contains)]
        // Can't use RangeInclusive::contains: it requires PartialOrd<Self>,
        // but Self::Coeff only provides PartialOrd<f64>.
        {
            debug_assert!(p >= 0.0 && p <= 1.0);
        }
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
        #[allow(clippy::manual_range_contains)]
        {
            debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
        }
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
    #[inline]
    fn two_qubit_pauli_error_impl(&mut self, addr0: usize, addr1: usize, p: [Self::Coeff; 15]) {
        if self.is_qubit_lost(addr0) || self.is_qubit_lost(addr1) {
            return;
        }
        #[allow(clippy::manual_range_contains)]
        {
            debug_assert!(p.iter().all(|p_| *p_ >= 0.0 && *p_ <= 1.0));
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
    /// two-qubit Pauli errors.
    #[inline]
    fn depolarize2_impl(&mut self, addr0: usize, addr1: usize, p: Self::Coeff) {
        if self.is_qubit_lost(addr0) || self.is_qubit_lost(addr1) {
            return;
        }
        #[allow(clippy::manual_range_contains)]
        {
            debug_assert!(p >= 0.0 && p <= 1.0);
        }
        let p_arr: [Self::Coeff; 15] = core::array::from_fn(|_| p.clone() * (1.0 / 15.0));
        self.two_qubit_pauli_error_impl(addr0, addr1, p_arr);
    }
}
