// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Stochastic channels: the Pauli-error family shared by both tableau types via
//! [`TableauLike`], plus the loss channels that only the generalized tableau has.
//!
//! Ported from `ppvm-tableau/src/{tableau_like,noise}.rs`.
//!
//! # RNG-consumption discipline
//!
//! The comparison boundaries and the *draw* discipline are gate-specific and are
//! reproduced verbatim, because they are what makes a seeded run reproducible:
//!
//! * [`TableauLike::depolarize_impl`] and [`TableauLike::pauli_error_impl`] draw
//!   **unconditionally**, even when the target is lost, and rely on the selected
//!   Clifford no-opping on the lost qubit — deliberately, "to preserve seeded RNG
//!   sequences across loss events".
//! * [`TableauLike::two_qubit_pauli_error_impl`] and
//!   [`TableauLike::depolarize2_impl`] instead return **early** when either qubit
//!   is lost, drawing nothing.
//! * The loss channels fire on `p >= r` (`if p < r { return; }`), the *opposite*
//!   strictness from the depolarizing family's `if p <= r { return; }`.
//!
//! The last two are internally inconsistent with the first (see this crate's
//! `# Deferrals`), but under the behaviour-preserving contract they are
//! reproduced exactly and flagged rather than "fixed": changing either would
//! silently desynchronize every seeded trajectory downstream.

use num::Zero;
use ppvm_traits_2::{
    AsymmetricLossChannel, Clifford, CorrelatedLossChannel, Depolarizing, Depolarizing2,
    LossChannel, Measure, Pauli, PauliError, ResetLossChannel, TwoQubitPauliError,
};
use rand::rngs::SmallRng;
use rand::{Rng, RngExt};

use crate::data::{Bitstring, GeneralizedTableau, RowStorage, Tableau};
use crate::measure::MeasureScratch;

/// `true` iff `p` lies in the closed unit interval `[0, 1]`.
#[inline]
fn is_probability(p: &f64) -> bool {
    if *p < 0.0 {
        return false;
    }
    *p <= 1.0
}

/// A stabilizer-tableau-like backend supporting Clifford gates and an RNG.
///
/// Implementing this grants the four Pauli noise channels through the default
/// bodies below. The associated `Rng` lets each backend choose its own
/// generator; nothing here depends on `SmallRng`.
///
/// Ported from `ppvm-tableau/src/tableau_like.rs`; the `Config`-derived
/// coefficient associated type collapses to `f64` (the probability domain), in
/// line with the crate-wide `Config` removal.
pub trait TableauLike: Clifford {
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
    /// The RNG is consumed **unconditionally**; the selected Clifford is
    /// expected to no-op on a lost qubit. This preserves seeded RNG sequences
    /// across loss events.
    #[inline]
    fn depolarize_impl(&mut self, addr0: usize, p: f64) {
        debug_assert!(is_probability(&p));
        let r = self.rng_mut().random::<f64>();
        if p <= r {
            return;
        }
        if p > r * 3.0 {
            // p / 3 > r >= 0
            self.x(addr0);
        } else if p > r * 1.5 {
            // 2p/3 > r >= p/3
            self.y(addr0);
        } else {
            // p > r >= 2p/3
            self.z(addr0);
        }
    }

    /// Single-qubit Pauli-error channel (`X`, `Y`, `Z` with given
    /// probabilities). The RNG is consumed unconditionally, as above.
    #[inline]
    fn pauli_error_impl(&mut self, addr0: usize, p: [f64; 3]) {
        debug_assert!(p.iter().all(is_probability));
        let r = self.rng_mut().random::<f64>();
        let mut cumulative = f64::zero();
        for (i, p_) in p.iter().enumerate() {
            cumulative += *p_;
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

    /// Two-qubit Pauli-error channel (15 non-identity combinations).
    ///
    /// Returns early — **without drawing** — if either qubit is lost.
    #[inline]
    fn two_qubit_pauli_error_impl(&mut self, addr0: usize, addr1: usize, p: [f64; 15]) {
        if self.is_qubit_lost(addr0) || self.is_qubit_lost(addr1) {
            return;
        }
        debug_assert!(p.iter().all(is_probability));
        let r = self.rng_mut().random::<f64>();
        let sum = f64::zero();
        let idx = p
            .iter()
            .scan(sum, |acc, p_| {
                *acc += *p_;
                Some(*acc)
            })
            .position(|cum_prob| cum_prob > r);

        if let Some(i) = idx {
            // Pairs indexed by `i + 1` (so `II` at index 0 is skipped).
            // Encoding: 0 = I, 1 = X, 2 = Y, 3 = Z; the first entry acts on
            // `addr0`, the second on `addr1`.
            #[rustfmt::skip]
            const PAULI_PAIRS: [(u8, u8); 16] = [
                (0,0),(0,1),(0,2),(0,3),
                (1,0),(1,1),(1,2),(1,3),
                (2,0),(2,1),(2,2),(2,3),
                (3,0),(3,1),(3,2),(3,3),
            ];
            let cartesian_index = PAULI_PAIRS[i + 1];

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
    /// two-qubit Pauli errors. Returns early — without drawing — on loss.
    #[inline]
    fn depolarize2_impl(&mut self, addr0: usize, addr1: usize, p: f64) {
        if self.is_qubit_lost(addr0) || self.is_qubit_lost(addr1) {
            return;
        }
        debug_assert!(is_probability(&p));
        let p_arr: [f64; 15] = core::array::from_fn(|_| p * (1.0 / 15.0));
        self.two_qubit_pauli_error_impl(addr0, addr1, p_arr);
    }
}

impl<A: RowStorage, H> TableauLike for Tableau<A, H> {
    type Rng = SmallRng;

    #[inline]
    fn rng_mut(&mut self) -> &mut Self::Rng {
        &mut self.rng
    }
}

impl<A: RowStorage, I: Bitstring, H> TableauLike for GeneralizedTableau<A, I, H> {
    type Rng = SmallRng;

    #[inline]
    fn rng_mut(&mut self) -> &mut Self::Rng {
        &mut self.tableau.rng
    }

    #[inline]
    fn is_qubit_lost(&self, addr: usize) -> bool {
        self.is_lost[addr]
    }
}

// The orphan rules forbid a blanket `impl<X: TableauLike> Depolarizing<f64> for
// X`, so the four channel traits are expanded per backend.
macro_rules! impl_tableau_noise {
    (generics: [$($gen:tt)*], ty: $ty:ty $(,)?) => {
        impl<A: RowStorage $($gen)*> Depolarizing<f64> for $ty {
            fn depolarize1(&mut self, qubit: usize, p: f64) {
                self.depolarize_impl(qubit, p);
            }
        }

        impl<A: RowStorage $($gen)*> PauliError<f64> for $ty {
            fn pauli_error(&mut self, qubit: usize, probabilities: [f64; 3]) {
                self.pauli_error_impl(qubit, probabilities);
            }
        }

        impl<A: RowStorage $($gen)*> TwoQubitPauliError<f64> for $ty {
            fn two_qubit_pauli_error(&mut self, qubit0: usize, qubit1: usize, p: [f64; 15]) {
                self.two_qubit_pauli_error_impl(qubit0, qubit1, p);
            }
        }

        impl<A: RowStorage $($gen)*> Depolarizing2<f64> for $ty {
            fn depolarize2(&mut self, qubit0: usize, qubit1: usize, p: f64) {
                self.depolarize2_impl(qubit0, qubit1, p);
            }
        }
    };
}

impl_tableau_noise! { generics: [, H], ty: Tableau<A, H> }
impl_tableau_noise! { generics: [, I: Bitstring, H], ty: GeneralizedTableau<A, I, H> }

// ─── Loss channels (generalized tableau only) ─────────────────────────────

impl<A: RowStorage, I: Bitstring, H> GeneralizedTableau<A, I, H> {
    /// Collapse one qubit for a loss event without retaining measurement
    /// scratch that the channel cannot reuse.
    #[inline]
    fn lose_qubit(&mut self, qubit: usize) {
        let outcome = if self.is_lost[qubit] {
            None
        } else {
            let (phase, stab, destab) = self.compute_decomposition(qubit, Pauli::Z);
            self.measure_with_scratch(
                qubit,
                &mut MeasureScratch::new(),
                phase,
                stab,
                destab,
                false,
            )
        };
        if let Some(true) = outcome {
            Clifford::x(self, qubit);
        }
        self.is_lost[qubit] = true;
    }
}

impl<A: RowStorage, I: Bitstring, H> LossChannel<f64> for GeneralizedTableau<A, I, H> {
    /// Lose `qubit` with probability `p`, then collapse it and reset to `|0⟩`.
    ///
    /// Note the comparison: the channel **fires when `p >= r`**, the opposite
    /// strictness from the depolarizing family. A loss event is not a logical
    /// measurement, so the record entry the internal `measure` pushed is popped
    /// — the channel is measurement-record-neutral.
    ///
    /// The convention split is adjudicated in `lean/PPVM/Algebra/Noise.lean`
    /// (§"The Bernoulli firing convention"): for an ideal `Uniform[0,1)` the two
    /// predicates differ only on the null event `{r = p}`
    /// (`fire_conventions_agree_off_diagonal`), so neither is a wrong
    /// `Bernoulli(p)` sampler — but at `p = 0` only the strict convention is a
    /// guaranteed no-op (`fire_strict_zero_noop` vs
    /// `fire_nonstrict_fires_at_zero`), so `loss_channel(q, 0.0)` fires with
    /// probability `2⁻⁵³` here while `depolarize1(q, 0.0)` never does.
    /// Reproduced verbatim from old under the behaviour-preservation directive;
    /// unifying the family is a seeded-stream change needing sign-off.
    fn loss_channel(&mut self, qubit: usize, p: f64) {
        if p < self.tableau.rng.random::<f64>() {
            return;
        }

        // O(n²), but it also potentially removes coefficients, which is nice.
        self.lose_qubit(qubit);
    }
}

impl<A: RowStorage, I: Bitstring, H> AsymmetricLossChannel<f64> for GeneralizedTableau<A, I, H> {
    /// State-dependent single-qubit loss.
    ///
    /// Models a three-level atom whose `|0⟩`/`|1⟩` levels leak into a loss state
    /// at different rates: `p_tot = p0·(1 + ⟨Z⟩)/2 + p1·(1 − ⟨Z⟩)/2`. With
    /// probability `p_tot` the qubit is collapsed, reset to `|0⟩`, and marked
    /// lost.
    ///
    /// # Approximation
    ///
    /// This is the trajectory *approximation*: it reproduces the loss statistics
    /// and is exact in the symmetric limit `p0 == p1`, but it does not apply the
    /// survival back-action `K₀ = √(1−p0)|0⟩⟨0| + √(1−p1)|1⟩⟨1|`, which is
    /// non-Clifford (it would branch the amplitude vector like an `rz`).
    ///
    /// # Known contract inconsistency (reproduced, not fixed)
    ///
    /// Unlike [`LossChannel::loss_channel`] and [`Reset::reset`], this does
    /// **not** pop the record entry the internal `measure` pushed, so an
    /// asymmetric-loss event pollutes the measurement record with a spurious
    /// `Some(bool)` — breaking the crate's own "a loss event is not a logical
    /// measurement" rule and shifting every `rec[-k]` lookback downstream. It is
    /// reproduced verbatim under the behaviour-preserving contract and listed in
    /// this crate's `# Deferrals` for sign-off.
    fn asymmetric_loss_channel(&mut self, qubit: usize, p0: f64, p1: f64) {
        if self.is_lost[qubit] {
            return;
        }
        let z = self.z_expectation(qubit);
        let p_tot = p0 * 0.5 * (1.0 + z) + p1 * 0.5 * (1.0 - z);

        if p_tot < self.tableau.rng.random::<f64>() {
            return;
        }
        if let Some(true) = Measure::measure(self, qubit) {
            Clifford::x(self, qubit);
        }
        self.is_lost[qubit] = true;
    }
}

impl<A: RowStorage, I: Bitstring, H> CorrelatedLossChannel<f64> for GeneralizedTableau<A, I, H> {
    /// Correlated loss on `(qubit0, qubit1)`.
    ///
    /// * `p[0]`: losing both simultaneously when both are in the qubit subspace.
    /// * `p[1]`: losing either one when both are in the qubit subspace.
    /// * `p[2]`: losing one when the other was already lost.
    ///
    /// Draws one `f64` for the `p[0]`/`p[1]` choice and, on the single-loss
    /// branch, one extra `bool` to pick which qubit. Goes through
    /// [`Reset::reset`], so it inherits the record pop.
    fn correlated_loss_channel(&mut self, qubit0: usize, qubit1: usize, p: [f64; 3]) {
        if self.is_lost[qubit0] {
            self.loss_channel(qubit1, p[2]);
            return;
        } else if self.is_lost[qubit1] {
            self.loss_channel(qubit0, p[2]);
            return;
        }

        let r = self.tableau.rng.random::<f64>();
        let mut cumulative = f64::zero();
        for (i, p_i) in p[..2].iter().enumerate() {
            cumulative += *p_i;
            if cumulative > r {
                if i == 0 {
                    // both lost
                    self.lose_qubit(qubit0);
                    self.lose_qubit(qubit1);
                } else {
                    // only a single qubit is lost
                    let choice = self.tableau.rng.random::<bool>();
                    if choice {
                        self.lose_qubit(qubit1);
                    } else {
                        self.lose_qubit(qubit0);
                    }
                }
                return;
            }
        }
    }
}

impl<A: RowStorage, I: Bitstring, H> ResetLossChannel for GeneralizedTableau<A, I, H> {
    /// Clear the loss bit only; the quantum state is untouched.
    fn reset_loss_channel(&mut self, qubit: usize) {
        self.is_lost[qubit] = false;
    }
}
