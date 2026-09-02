// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::fmt::Debug;
use std::ops::Add;

use bitvec::view::BitView;
use num::PrimInt;
use num::complex::{Complex, Complex64, ComplexFloat};
use num::traits::{One, ToPrimitive, Zero};

use crate::prelude::*;
use rand::RngExt;
use rand::rngs::SmallRng;

// === TableauLike impls ===
//
// Implementing TableauLike grants automatic implementations of all
// single- and two-qubit Pauli noise channels via default methods.

impl<T: Config> TableauLike for Tableau<T>
where
    T::Coeff: PartialOrd<f64>,
    // The canonical word-level `Clifford for Tableau<T>` impl (required by the
    // `TableauLike: Clifford` supertrait) operates on raw storage words, so it
    // carries `Store: PrimInt` — the same bound the `GeneralizedTableau` side
    // already requires.
    <T::Storage as BitView>::Store: PrimInt,
{
    type Coeff = T::Coeff;
    type Rng = SmallRng;

    #[inline]
    fn rng_mut(&mut self) -> &mut Self::Rng {
        &mut self.rng
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> TableauLike
    for GeneralizedTableau<T, I, C>
where
    T::Coeff: PartialOrd<f64>,
    Complex<T::Coeff>: From<Complex<f64>>,
    <T::Storage as BitView>::Store: PrimInt,
{
    type Coeff = T::Coeff;
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

/// `true` iff `p` is an admissible correlated-loss parameter triple.
///
/// The channel is completely positive exactly when every event weight is
/// nonnegative, i.e. `p0, p1 >= 0`, `p0 + 2·p1 <= 1` (`p1` is the probability
/// that a *named* one of the pair is lost, so the exactly-one event carries
/// `2·p1`) and `p2 ∈ [0, 1]`. Outside that region the mixture truncates a
/// negative survivor weight and renormalizes, and the trajectory's cumulative
/// scan stops being a categorical sampler — both silently.
///
/// Coefficients are checked in `Coeff` space via `PartialOrd<f64>` (the
/// bound every loss-channel impl already carries). Slack still applies
/// to the `p0 + 2·p1 <= 1` sum so a saturated `[1/3, 1/3, _]` is not
/// rejected by last-bit rounding.
pub fn is_admissible_correlated_loss<C>(p: &[C; 3]) -> bool
where
    C: Clone + PartialOrd<f64> + Add<Output = C>,
{
    // The slack matches `ppvm-tableau-2`'s guard: reject genuinely out-of-region
    // parameters, not last-bit noise in a legitimately saturated triple. An exact
    // `<= 1.0` would spuriously fire on e.g. `[1/3, 1/3, _]`, where
    // `p0 + 2·p1` rounds to `1.0000000000000002`.
    const SLACK: f64 = 1e-9;
    p[0] >= 0.0
        && p[1] >= 0.0
        && p[2] >= 0.0
        && p[2] <= 1.0
        && p[0].clone() + p[1].clone() + p[1].clone() <= 1.0 + SLACK
}

// === Noise trait impls ===
//
// Orphan rules (E0210) forbid `impl<X: TableauLike<Coeff = T::Coeff>> Depolarizing<T> for X`,
// so we expand the four noise traits per backend via a macro. Each backend only
// has to list its generics + where-clause once.

macro_rules! impl_tableau_noise {
    (generics: [$($gen:tt)*], ty: $ty:ty, where: [$($bound:tt)*] $(,)?) => {
        impl<T: Config $($gen)*> Depolarizing<T> for $ty
        where $($bound)*
        {
            fn depolarize1(&mut self, addr0: usize, p: T::Coeff) {
                self.depolarize_impl(addr0, p);
            }
        }

        impl<T: Config $($gen)*> PauliError<T> for $ty
        where $($bound)*
        {
            fn pauli_error(&mut self, addr0: usize, p: [T::Coeff; 3]) {
                self.pauli_error_impl(addr0, p);
            }
        }

        impl<T: Config $($gen)*> TwoQubitPauliError<T> for $ty
        where $($bound)*
        {
            fn two_qubit_pauli_error(&mut self, addr0: usize, addr1: usize, p: [T::Coeff; 15]) {
                self.two_qubit_pauli_error_impl(addr0, addr1, p);
            }
        }

        impl<T: Config $($gen)*> Depolarizing2<T> for $ty
        where $($bound)*
        {
            fn depolarize2(&mut self, addr0: usize, addr1: usize, p: T::Coeff) {
                self.depolarize2_impl(addr0, addr1, p);
            }
        }
    };
}

impl_tableau_noise! {
    generics: [],
    ty: Tableau<T>,
    where: [
        T::Coeff: PartialOrd<f64>,
        <T::Storage as BitView>::Store: PrimInt,
    ],
}

impl_tableau_noise! {
    generics: [, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>],
    ty: GeneralizedTableau<T, I, C>,
    where: [
        T::Coeff: PartialOrd<f64>,
        Complex<T::Coeff>: From<Complex<f64>>,
        <T::Storage as BitView>::Store: PrimInt,
    ],
}

// === GeneralizedTableau-specific loss channels (no Tableau equivalent) ===

impl<T: Config, I: TableauIndex + Send + Sync, C: SparseVector<Complex<T::Coeff>, I>> LossChannel<T>
    for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn loss_channel(&mut self, addr0: usize, p: <T as Config>::Coeff) {
        if p < self.tableau.rng.random::<f64>() {
            return;
        }

        // NOTE: this is O(n^2) but also potentially removes coefficients, which is nice
        let outcome = self.measure(addr0);
        // A loss event is not a logical measurement: keep the measurement
        // record neutral by dropping the entry the internal `measure` pushed.
        self.measurement_record.pop();
        if let Some(true) = outcome {
            // flip back to 0
            self.x(addr0);
        }
        self.is_lost[addr0] = true;
    }
}

impl<T: Config, I: TableauIndex + Send + Sync, C: SparseVector<Complex<T::Coeff>, I>>
    AsymmetricLossChannel<T> for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    /// State-dependent single-qubit loss ("asymmetric loss").
    ///
    /// Models a three-level atom whose `|0⟩` and `|1⟩` levels leak into a loss
    /// state `|L⟩` at different rates: the qubit is lost from `|0⟩` with
    /// probability `p0` and from `|1⟩` with `p1`. The total loss probability is
    /// state-dependent,
    ///
    /// ```text
    ///     p_tot = p0 * (1 + ⟨Z⟩)/2 + p1 * (1 - ⟨Z⟩)/2,
    /// ```
    ///
    /// with `⟨Z⟩` the current Z-expectation of `addr0`. With probability `p_tot`
    /// the qubit is collapsed, reset to `|0⟩`, and marked lost (as in
    /// [`LossChannel::loss_channel`]); otherwise it is left unchanged.
    ///
    /// # Approximation
    ///
    /// This is the trajectory *approximation* of the true loss channel. It
    /// reproduces the loss statistics (which qubits are flagged lost, and how
    /// often) and is exact in the symmetric limit `p0 == p1`, where it reduces
    /// to [`LossChannel::loss_channel`]. It does NOT apply the survival
    /// back-action `K0 = sqrt(1-p0)|0⟩⟨0| + sqrt(1-p1)|1⟩⟨1|`: for `p0 != p1`
    /// the faithful channel reshapes the *surviving* qubit (population tilts
    /// toward the less-leaky level, coherences are damped), and this
    /// implementation skips that reshaping. The back-action is non-Clifford (it
    /// branches the coefficient vector like an `rz`), so it is omitted to keep
    /// the channel cheap enough to apply after every gate. See issue #39.
    fn asymmetric_loss_channel(&mut self, addr0: usize, p0: T::Coeff, p1: T::Coeff) {
        if self.is_lost[addr0] {
            return;
        }
        // State-dependent loss probability from the populations pop0/pop1.
        let z = self.z_expectation(addr0);
        let p_tot = p0.to_f64().unwrap() * 0.5 * (1.0 + z) + p1.to_f64().unwrap() * 0.5 * (1.0 - z);

        if p_tot < self.tableau.rng.random::<f64>() {
            return;
        }
        // Lost: collapse + reset to |0⟩, mirroring loss_channel.
        if let Some(true) = self.measure(addr0) {
            self.x(addr0);
        }
        self.is_lost[addr0] = true;
    }
}

impl<T: Config, I: TableauIndex + Send + Sync, C: SparseVector<Complex<T::Coeff>, I>>
    CorrelatedLossChannel<T> for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd<f64>
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug
        + Send
        + Sync,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    /// Apply a correlated loss channel to qubits at `addr0` and `addr1`.
    ///
    /// The three probabilities are:
    /// * `p[0]`: The probability of losing both qubits simultaneously when
    ///   both of them are in the qubit subspace.
    /// * `p[1]`: The probability of losing a **named** one of the two qubits when
    ///   both of them are in the qubit subspace, so the probability of losing
    ///   *exactly one* is `2·p[1]` and the both-present survivor keeps
    ///   `1 − 2·p[1] − p[0]` (which qubit is lost is 50/50). See
    ///   [`ppvm_traits::traits::CorrelatedLossChannel`] for the normative
    ///   statement.
    /// * `p[2]`: The probability of losing one qubit when the other one has already
    ///   been lost prior to the channel.
    fn correlated_loss_channel(
        &mut self,
        addr0: usize,
        addr1: usize,
        p: [<T as Config>::Coeff; 3],
    ) {
        debug_assert!(
            is_admissible_correlated_loss(&p),
            "correlated loss needs p0, p1 >= 0, p0 + 2*p1 <= 1, p2 in [0, 1]; got {p:?}"
        );
        if self.is_lost[addr0] {
            self.loss_channel(addr1, p[2].clone());
            return;
        } else if self.is_lost[addr1] {
            self.loss_channel(addr0, p[2].clone());
            return;
        }

        let r = self.tableau.rng.random::<f64>();
        let mut cumulative = T::Coeff::zero();
        for (i, p_i) in p[..2].iter().enumerate() {
            // `p[1]` is the probability that a *named* one of the pair is lost,
            // so the exactly-one event carries `2·p[1]` in this categorical
            // scan; the fair coin below then picks which qubit.
            cumulative += if i == 1 {
                p_i.clone() + p_i.clone()
            } else {
                p_i.clone()
            };
            if cumulative > r {
                if i == 0 {
                    // both lost
                    self.reset(addr0);
                    self.reset(addr1);
                    self.is_lost[addr0] = true;
                    self.is_lost[addr1] = true;
                } else {
                    // only losing a single qubit,
                    let choice = self.tableau.rng.random::<bool>();
                    if choice {
                        self.reset(addr1);
                        self.is_lost[addr1] = true;
                    } else {
                        self.reset(addr0);
                        self.is_lost[addr0] = true;
                    }
                }
                return;
            }
        }
    }
}

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> ResetLossChannel<T>
    for GeneralizedTableau<T, I, C>
{
    fn reset_loss_channel(&mut self, addr0: usize) {
        self.is_lost[addr0] = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ppvm_pauli_sum::config::fxhash::ByteF64;

    type TestConfig = ByteF64<1>;
    type TestTab = GeneralizedTableau<TestConfig>;

    fn tab(n: usize) -> TestTab {
        GeneralizedTableau::new(n, 1e-12)
    }

    // === Depolarizing ===

    #[test]
    fn depolarize_p0_no_change() {
        let mut t = tab(1);
        t.depolarize1(0, 0.0);
        assert!(!t.measure(0).unwrap());
    }

    #[test]
    fn depolarize_p1_does_not_mark_lost() {
        // With p=1.0 an error is always applied; verify is_lost is unaffected
        let mut t = tab(1);
        t.depolarize1(0, 1.0);
        assert!(!t.is_lost[0]);
    }

    // === PauliError ===

    #[test]
    fn pauli_error_zero_prob_no_change() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 0.0, 0.0]);
        assert!(!t.measure(0).unwrap());
    }

    #[test]
    fn pauli_error_x_flips_qubit() {
        let mut t = tab(1);
        t.pauli_error(0, [1.0, 0.0, 0.0]); // X|0⟩ = |1⟩
        assert!(t.measure(0).unwrap());
    }

    #[test]
    fn pauli_error_y_flips_qubit() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 1.0, 0.0]); // Y|0⟩ = i|1⟩
        assert!(t.measure(0).unwrap());
    }

    #[test]
    fn pauli_error_z_no_measurement_change() {
        let mut t = tab(1);
        t.pauli_error(0, [0.0, 0.0, 1.0]); // Z|0⟩ = -|0⟩, still measures 0
        assert!(!t.measure(0).unwrap());
    }

    #[test]
    fn pauli_error_x_on_excited_qubit_flips_back() {
        let mut t = tab(1);
        t.x(0); // |1⟩
        t.pauli_error(0, [1.0, 0.0, 0.0]); // X|1⟩ = |0⟩
        assert!(!t.measure(0).unwrap());
    }

    // === TwoQubitPauliError ===

    #[test]
    fn two_qubit_pauli_error_zero_prob_no_change() {
        let mut t = tab(2);
        t.two_qubit_pauli_error(0, 1, [0.0; 15]);
        assert!(!t.measure(0).unwrap());
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn two_qubit_pauli_error_ix_flips_second_only() {
        // p[0] = 1.0 → IX: I on addr0, X on addr1
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[0] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(0).unwrap());
        assert!(t.measure(1).unwrap());
    }

    #[test]
    fn two_qubit_pauli_error_xi_flips_first_only() {
        // p[3] = 1.0 → XI: X on addr0, I on addr1
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[3] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.measure(0).unwrap());
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn two_qubit_pauli_error_xx_flips_both() {
        // p[4] = 1.0 → XX
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[4] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.measure(0).unwrap());
        assert!(t.measure(1).unwrap());
    }

    #[test]
    fn two_qubit_pauli_error_zz_no_measurement_change() {
        // p[14] = 1.0 → ZZ: Z|0⟩ = -|0⟩ on both, still measures 0
        let mut t = tab(2);
        let mut p = [0.0f64; 15];
        p[14] = 1.0;
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(0).unwrap());
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn two_qubit_pauli_error_both_lost_no_change() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.is_lost[1] = true;
        let mut p = [0.0f64; 15];
        p[4] = 1.0; // XX — skipped entirely
        t.two_qubit_pauli_error(0, 1, p);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn two_qubit_pauli_error_first_lost_no_apply() {
        // addr0 lost; p[0] = 1.0 (IX) → marginal p_x for addr1 = 1.0
        let mut t = tab(2);
        t.is_lost[0] = true;
        let mut p = [0.0f64; 15];
        p[0] = 1.0; // IX
        t.two_qubit_pauli_error(0, 1, p);
        assert!(!t.measure(1).unwrap()); // nothing applied to addr1
    }

    // === Depolarizing2 ===

    #[test]
    fn depolarize2_p0_no_change() {
        let mut t = tab(2);
        t.depolarize2(0, 1, 0.0);
        assert!(!t.measure(0).unwrap());
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn depolarize2_both_lost_no_change() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.is_lost[1] = true;
        t.depolarize2(0, 1, 1.0);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn depolarize2_first_lost_p0_second_unchanged() {
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.depolarize2(0, 1, 0.0); // effective p on addr1 = 4/5 * 0 = 0
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn depolarize2_second_lost_p0_first_unchanged() {
        let mut t = tab(2);
        t.is_lost[1] = true;
        t.depolarize2(0, 1, 0.0); // effective p on addr0 = 4/5 * 0 = 0
        assert!(!t.measure(0).unwrap());
    }

    // === LossChannel ===

    #[test]
    fn loss_channel_p0_qubit_not_lost() {
        let mut t = tab(1);
        t.loss_channel(0, 0.0);
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn loss_channel_p1_qubit_marked_lost() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
    }

    #[test]
    fn loss_channel_p1_qubit_reset_to_zero() {
        // Qubit starts in |1⟩; loss_channel should measure, reset to |0⟩, then mark lost
        let mut t = tab(1);
        t.x(0);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
        assert!(t.measure(0).is_none()); // Reset to |0⟩ before marking lost
    }

    #[test]
    fn loss_channel_does_not_pollute_measurement_record() {
        // A loss event is not a logical measurement and must leave the
        // measurement record untouched.
        let mut t = tab(1);
        t.x(0);
        t.loss_channel(0, 1.0);
        assert!(t.current_measurement_record().is_empty());
    }

    #[test]
    fn loss_channel_p1_subsequent_gate_is_noop() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        t.x(0); // No-op: qubit is lost
        assert!(t.measure(0).is_none());
        t.is_lost[0] = false;
        assert!(!t.measure(0).unwrap()); // still 0
    }

    #[test]
    fn loss_channel_p0_second_qubit_unaffected() {
        let mut t = tab(2);
        t.loss_channel(0, 0.0);
        t.loss_channel(1, 0.0);
        assert!(!t.is_lost[0]);
        assert!(!t.is_lost[1]);
    }

    // === ResetLossChannel ===

    #[test]
    fn reset_loss_channel_clears_lost_flag() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
        t.reset_loss_channel(0);
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn reset_loss_channel_qubit_in_ground_state() {
        // loss_channel resets qubit to |0⟩; after reset_loss_channel it should still be |0⟩
        let mut t = tab(1);
        t.x(0); // |1⟩
        t.loss_channel(0, 1.0); // measures, resets to |0⟩, marks lost
        t.reset_loss_channel(0);
        assert!(!t.measure(0).unwrap()); // back in |0⟩
    }

    #[test]
    fn reset_loss_channel_gates_work_again() {
        let mut t = tab(1);
        t.loss_channel(0, 1.0);
        t.reset_loss_channel(0);
        t.x(0); // should no longer be a no-op
        assert!(t.measure(0).unwrap());
    }

    // === Seeded RNG ordering ===
    //
    // `Depolarizing` and `PauliError` must consume RNG unconditionally so
    // that seeded traces are reproducible regardless of loss events. The
    // selected Clifford gate no-ops on lost qubits (see gates/clifford.rs).

    #[test]
    fn depolarize_rng_consumed_on_lost_qubit() {
        let seed = 42u64;
        let mut t_active = tab(1);
        t_active.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
        t_active.depolarize1(0, 0.3);
        let next_active: f64 = t_active.tableau.rng.random();

        let mut t_lost = tab(1);
        t_lost.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
        t_lost.is_lost[0] = true;
        t_lost.depolarize1(0, 0.3);
        let next_lost: f64 = t_lost.tableau.rng.random();

        assert_eq!(next_active, next_lost);
    }

    #[test]
    fn pauli_error_rng_consumed_on_lost_qubit() {
        let seed = 42u64;
        let mut t_active = tab(1);
        t_active.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
        t_active.pauli_error(0, [0.1, 0.1, 0.1]);
        let next_active: f64 = t_active.tableau.rng.random();

        let mut t_lost = tab(1);
        t_lost.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
        t_lost.is_lost[0] = true;
        t_lost.pauli_error(0, [0.1, 0.1, 0.1]);
        let next_lost: f64 = t_lost.tableau.rng.random();

        assert_eq!(next_active, next_lost);
    }

    // === Statistical tests ===

    #[test]
    fn depolarize_statistics() {
        // Starting from |0⟩, P(measure 1) = P(X) + P(Y) = p/3 + p/3 = 2p/3.
        // Z leaves |0⟩ unchanged; I leaves |0⟩ unchanged.
        let p = 0.6_f64;
        let expected = 2.0 * p / 3.0; // 0.4
        let trials = 500;

        let ones = (0..trials)
            .filter(|_| {
                let mut t = tab(1);
                t.depolarize1(0, p);
                t.measure(0).unwrap()
            })
            .count();

        let fraction = ones as f64 / trials as f64;
        // tolerance ~5σ: σ = sqrt(expected*(1-expected)/trials) ≈ 0.022
        assert!(
            (fraction - expected).abs() < 0.1,
            "Expected fraction {expected:.3}, got {fraction:.3}"
        );
    }

    #[test]
    fn depolarize2_statistics() {
        // Starting from |00⟩, errors that flip qubit 0 to |1⟩ are X and Y on that qubit:
        // XI, XX, XY, XZ, YI, YX, YY, YZ — 8 out of 15, so P(q0=1) = 8p/15.
        let p = 0.6_f64;
        let expected = 8.0 * p / 15.0; // 0.32
        let trials = 500;

        let ones = (0..trials)
            .filter(|_| {
                let mut t = tab(2);
                t.depolarize2(0, 1, p);
                t.measure(0).unwrap()
            })
            .count();

        let fraction = ones as f64 / trials as f64;
        // tolerance ~5σ: σ = sqrt(expected*(1-expected)/trials) ≈ 0.021
        assert!(
            (fraction - expected).abs() < 0.1,
            "Expected fraction {expected:.3}, got {fraction:.3}"
        );
    }

    #[test]
    fn test_cnot() {
        let mut t = tab(2);
        t.x(0);
        t.cnot(0, 1);
        t.loss_channel(0, 1.0);
        assert!(t.measure(0).is_none());
        assert!(t.measure(1).unwrap());

        let mut t = tab(2);
        t.loss_channel(0, 1.0);
        t.x(0);
        t.cnot(0, 1);
        assert!(!t.measure(1).unwrap());
        assert!(t.measure(0).is_none());
    }

    #[test]
    fn test_ghz_statistics() {
        let mut t = tab(2);
        t.h(0);
        t.cnot(0, 1);

        let trials = 100u64;
        let mut z_avg = 0.0;
        let p = 0.1;
        for i in 0..trials {
            let mut t_trial = t.fork(Some(i));
            t_trial.loss_channel(0, p);

            let outcome0 = t_trial.measure(0);
            let outcome1 = t_trial.measure(1);
            if outcome0.unwrap_or(false) == outcome1.unwrap_or(false) {
                z_avg += 1.0 / trials as f64;
            } else {
                z_avg += -1.0 / trials as f64;
            }
        }

        println!("{}", z_avg);
        assert!((z_avg - (1.0 - p)).abs() < 10.0 / trials as f64);
    }

    // === CorrelatedLossChannel ===

    #[test]
    fn correlated_loss_p0_no_loss() {
        // All probabilities zero: neither qubit should be lost.
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0]);
        assert!(!t.is_lost[0]);
        assert!(!t.is_lost[1]);
    }

    #[test]
    fn correlated_loss_p0_both_lost() {
        // p[0]=1 → both qubits always lost.
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0]);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn correlated_loss_p1_exactly_one_lost() {
        // p[1]=0.5 → exactly one qubit lost each time, since `p[1]` is the
        // probability that a *named* one of the pair is lost and so the
        // exactly-one event carries 2·p[1] = 1. (`[0.0, 1.0, 0.0]`, which this
        // test used before, is inadmissible under that convention:
        // p0 + 2·p1 = 2 > 1.)
        let trials = 200;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [0.0, 0.5, 0.0]);
            assert!(
                t.is_lost[0] ^ t.is_lost[1],
                "Expected exactly one lost qubit (seed {seed})"
            );
        }
    }

    #[test]
    fn correlated_loss_p1_both_qubits_chosen_equally() {
        // With 2·p[1]=1 the coin flip should lose addr0 and addr1 with equal
        // frequency.
        let trials = 1000u64;
        let mut addr0_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [0.0, 0.5, 0.0]);
            if t.is_lost[0] {
                addr0_lost += 1;
            }
        }
        let fraction = addr0_lost as f64 / trials as f64;
        // Expected 0.5; 5σ tolerance with σ ≈ 0.016
        assert!(
            (fraction - 0.5).abs() < 0.08,
            "Expected ~0.5, got {fraction:.3}"
        );
    }

    #[test]
    fn correlated_loss_both_lost_resets_to_zero() {
        // When both qubits are lost their state should have been reset to |0⟩.
        let mut t = tab(2);
        t.x(0);
        t.x(1);
        t.correlated_loss_channel(0, 1, [1.0, 0.0, 0.0]);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
        // Restore so we can measure.
        t.is_lost[0] = false;
        t.is_lost[1] = false;
        assert!(!t.measure(0).unwrap());
        assert!(!t.measure(1).unwrap());
    }

    #[test]
    fn correlated_loss_single_lost_resets_to_zero() {
        // The lost qubit should be in |0⟩; the surviving qubit keeps its state.
        // Use a seed where addr0 ends up being the lost one.
        // We iterate seeds until we get addr0 lost, then verify.
        for seed in 0..1000u64 {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.x(0); // put addr0 in |1⟩
            // 2·p[1] = 1: exactly one qubit is lost in every trial.
            t.correlated_loss_channel(0, 1, [0.0, 0.5, 0.0]);
            if t.is_lost[0] {
                t.is_lost[0] = false;
                assert!(!t.measure(0).unwrap(), "Lost qubit should be reset to |0⟩");
                return;
            }
        }
        panic!("addr0 was never chosen as the lost qubit in 1000 trials");
    }

    #[test]
    fn correlated_loss_addr0_already_lost_applies_p2_to_addr1() {
        // addr0 already lost → addr1 should be lost with probability p[2]=1.
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 1.0]);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn correlated_loss_addr1_already_lost_applies_p2_to_addr0() {
        // addr1 already lost → addr0 should be lost with probability p[2]=1.
        let mut t = tab(2);
        t.is_lost[1] = true;
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 1.0]);
        assert!(t.is_lost[0]);
        assert!(t.is_lost[1]);
    }

    #[test]
    fn correlated_loss_addr0_already_lost_p2_zero_addr1_survives() {
        // addr0 already lost, p[2]=0 → addr1 stays active.
        let mut t = tab(2);
        t.is_lost[0] = true;
        t.correlated_loss_channel(0, 1, [0.0, 0.0, 0.0]);
        assert!(!t.is_lost[1]);
    }

    #[test]
    fn correlated_loss_statistics_both() {
        // P(both lost) should converge to p[0].
        let p_both = 0.3_f64;
        let trials = 1000u64;
        let mut both_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [p_both, 0.0, 0.0]);
            if t.is_lost[0] && t.is_lost[1] {
                both_lost += 1;
            }
        }
        let fraction = both_lost as f64 / trials as f64;
        // 5σ tolerance: σ = sqrt(0.3*0.7/1000) ≈ 0.014
        assert!(
            (fraction - p_both).abs() < 0.07,
            "Expected ~{p_both:.2}, got {fraction:.3}"
        );
    }

    #[test]
    fn correlated_loss_statistics_single() {
        // `p[1]` is the probability that a *named* one of the pair is lost, so
        // P(exactly one lost) converges to 2·p[1]. (This test previously
        // asserted `p[1]`, i.e. the rejected convention.)
        let p_single = 0.2_f64;
        let expected = 2.0 * p_single;
        let trials = 1000u64;
        let mut one_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [0.0, p_single, 0.0]);
            if t.is_lost[0] ^ t.is_lost[1] {
                one_lost += 1;
            }
        }
        let fraction = one_lost as f64 / trials as f64;
        // 5σ: σ = sqrt(0.4*0.6/1000) ≈ 0.015
        assert!(
            (fraction - expected).abs() < 0.08,
            "Expected ~{expected:.2}, got {fraction:.3}"
        );
    }

    // === G-040 — the correlated-loss `p[1]` convention ===
    //
    // The paper (`ppvm-paper/main.tex:462`, `:523`, `:845`) is the definition of
    // record: `p[1]` is `p_LQ`, the probability that a **named** one of the pair
    // is lost, so `P(exactly one lost) = 2·p[1]` and the both-present survivor
    // scales by `1 − 2·p[1] − p[0]`. `ppvm-pauli-sum` already reads it that way;
    // these tests pin this trajectory sampler to the same number, because the
    // cross-backend disagreement is the actual bug a Python user hits.

    type LossyTestSum = PauliSum<
        ppvm_pauli_sum::config::fxhash::Byte<
            1,
            f64,
            NoStrategy,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    >;

    /// The Heisenberg scale factor `ppvm-pauli-sum` applies to a fully
    /// in-subspace observable, i.e. `1 − p[0] − P(exactly one lost)`.
    fn pauli_sum_survivor(p: [f64; 3]) -> f64 {
        let mut sum = LossyTestSum::builder().n_qubits(2).build();
        sum += ("ZZ", 1.0);
        sum.correlated_loss_channel(0, 1, p);
        let zz: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "ZZ".into();
        *sum.data()
            .get(&zz)
            .expect("the all-present term survives a pure rescale")
    }

    #[test]
    fn correlated_loss_exactly_one_lost_is_two_p1_and_agrees_with_pauli_sum() {
        let p1 = 0.3_f64;
        let p = [0.0, p1, 0.0];
        let trials = 20_000u64;
        let mut one_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, p);
            if t.is_lost[0] ^ t.is_lost[1] {
                one_lost += 1;
            }
        }
        let trajectory = one_lost as f64 / trials as f64;
        // Same number, read off the (already correct) Heisenberg backend.
        let pauli_sum = 1.0 - pauli_sum_survivor(p);
        assert!(
            (pauli_sum - 2.0 * p1).abs() < 1e-12,
            "ppvm-pauli-sum P(exactly one) = {pauli_sum}, want 2*p[1] = {}",
            2.0 * p1
        );
        assert!(
            (trajectory - 2.0 * p1).abs() < 0.02,
            "trajectory P(exactly one lost) = {trajectory:.4}, want 2*p[1] = {} \
             (ppvm-pauli-sum says {pauli_sum})",
            2.0 * p1
        );
    }

    /// The paper's transport prediction: with `[p/3, p/3, p/3]` the per-event
    /// both-present survival is `1 − p[0] − 2·p[1] = 1 − p`.
    #[test]
    fn correlated_loss_paper_transport_survival_is_one_minus_p() {
        let p = 0.15_f64;
        let third = p / 3.0;
        let probabilities = [third, third, third];
        let trials = 20_000u64;
        let mut none_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, probabilities);
            if !t.is_lost[0] && !t.is_lost[1] {
                none_lost += 1;
            }
        }
        let fraction = none_lost as f64 / trials as f64;
        assert!(
            (fraction - (1.0 - p)).abs() < 0.02,
            "trajectory P(no loss) = {fraction:.4}, want 1 - {p} = {}",
            1.0 - p
        );
        let survivor = pauli_sum_survivor(probabilities);
        assert!(
            (survivor - (1.0 - p)).abs() < 1e-12,
            "ppvm-pauli-sum survivor {survivor}, want 1 - {p}"
        );
    }

    // G-043 — the admissible region `p0, p1 >= 0`, `p0 + 2·p1 <= 1`,
    // `p2 ∈ [0, 1]` is what makes the channel completely positive. Outside it
    // the cumulative scan below silently stops being a categorical sampler.

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "p0 + 2*p1 <= 1")]
    fn correlated_loss_rejects_inadmissible_probabilities() {
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [0.6, 0.6, 0.0]);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "p0 + 2*p1 <= 1")]
    fn correlated_loss_rejects_negative_probabilities() {
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [5.0, -3.0, 17.0]);
    }

    /// Coeff stand-in that cannot be converted to `f64`; used to ensure the helper
    /// only relies on `PartialOrd<f64>` comparisons.
    #[derive(Clone, Copy)]
    struct NoF64(f64);

    impl ToPrimitive for NoF64 {
        fn to_i64(&self) -> Option<i64> {
            None
        }
        fn to_u64(&self) -> Option<u64> {
            None
        }
        fn to_f64(&self) -> Option<f64> {
            None
        }
    }

    impl PartialEq<f64> for NoF64 {
        fn eq(&self, other: &f64) -> bool {
            self.0 == *other
        }
    }

    impl PartialOrd<f64> for NoF64 {
        fn partial_cmp(&self, other: &f64) -> Option<std::cmp::Ordering> {
            self.0.partial_cmp(other)
        }
    }

    impl std::ops::Add for NoF64 {
        type Output = Self;
        fn add(self, rhs: Self) -> Self {
            Self(self.0 + rhs.0)
        }
    }

    #[test]
    fn is_admissible_rejects_negative_p1() {
        assert!(!is_admissible_correlated_loss(&[0.0, -0.1, 0.0]));
    }

    #[test]
    fn is_admissible_accepts_f64_saturated_third() {
        assert!(is_admissible_correlated_loss(&[1.0 / 3.0, 1.0 / 3.0, 0.5]));
    }

    #[test]
    fn is_admissible_rejects_when_to_f64_fails() {
        let p = [NoF64(0.6), NoF64(0.6), NoF64(0.0)];
        assert!(
            !is_admissible_correlated_loss(&p),
            "p0 + 2·p1 = 1.8 must be inadmissible even without an f64 conversion"
        );
    }

    #[test]
    fn is_admissible_accepts_saturated_boundary_when_to_f64_fails() {
        let p = [NoF64(1.0 / 3.0), NoF64(1.0 / 3.0), NoF64(0.5)];
        assert!(
            is_admissible_correlated_loss(&p),
            "[1/3, 1/3, _] must stay admissible without an f64 conversion"
        );
    }

    /// The saturated boundary `p0 + 2·p1 == 1` is admissible and must not trip
    /// the guard — including the `[1/3, 1/3, _]` triple whose sum rounds above
    /// one in binary floating point.
    #[test]
    fn correlated_loss_saturated_boundary_is_admissible() {
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [0.2, 0.4, 1.0]);
        let mut t = tab(2);
        t.correlated_loss_channel(0, 1, [1.0 / 3.0, 1.0 / 3.0, 0.5]);
    }

    // === z_expectation ===

    #[test]
    fn z_expectation_ground_state_is_plus_one() {
        let t = tab(1);
        assert!((t.z_expectation(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn z_expectation_excited_state_is_minus_one() {
        let mut t = tab(1);
        t.x(0);
        assert!((t.z_expectation(0) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn z_expectation_superposition_is_zero() {
        let mut t = tab(1);
        t.h(0);
        assert!(t.z_expectation(0).abs() < 1e-12);
    }

    // === AsymmetricLossChannel ===

    #[test]
    fn asymmetric_loss_ground_state_uses_p0() {
        // |0⟩: pop0 = 1, so p_tot = p0.
        let mut t = tab(1);
        t.asymmetric_loss_channel(0, 1.0, 0.0);
        assert!(t.is_lost[0]);

        let mut t = tab(1);
        t.asymmetric_loss_channel(0, 0.0, 1.0); // p_tot = 0
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn asymmetric_loss_excited_state_uses_p1() {
        // |1⟩: pop1 = 1, so p_tot = p1.
        let mut t = tab(1);
        t.x(0);
        t.asymmetric_loss_channel(0, 0.0, 1.0);
        assert!(t.is_lost[0]);

        let mut t = tab(1);
        t.x(0);
        t.asymmetric_loss_channel(0, 1.0, 0.0); // p_tot = 0
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn asymmetric_loss_zero_prob_not_lost() {
        let mut t = tab(1);
        t.asymmetric_loss_channel(0, 0.0, 0.0);
        assert!(!t.is_lost[0]);
    }

    #[test]
    fn asymmetric_loss_already_lost_is_noop() {
        let mut t = tab(1);
        t.is_lost[0] = true;
        t.asymmetric_loss_channel(0, 1.0, 1.0);
        assert!(t.is_lost[0]);
    }

    #[test]
    fn asymmetric_loss_resets_lost_qubit_to_zero() {
        // |1⟩ lost with p1 = 1; after un-marking it should read |0⟩.
        let mut t = tab(1);
        t.x(0);
        t.asymmetric_loss_channel(0, 0.0, 1.0);
        assert!(t.is_lost[0]);
        t.is_lost[0] = false;
        assert!(!t.measure(0).unwrap());
    }

    #[test]
    fn asymmetric_loss_symmetric_matches_loss_channel() {
        // p0 == p1 == p reduces to loss_channel(p); on |+⟩, p_tot = p.
        let p = 0.3;
        let trials = 1000u64;
        let mut lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(1);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.h(0);
            t.asymmetric_loss_channel(0, p, p);
            if t.is_lost[0] {
                lost += 1;
            }
        }
        let frac = lost as f64 / trials as f64;
        assert!((frac - p).abs() < 0.07, "expected ~{p}, got {frac:.3}");
    }

    #[test]
    fn asymmetric_loss_superposition_averages_probs() {
        // |+⟩: ⟨Z⟩ = 0 so p_tot = (p0 + p1) / 2.
        let (p0, p1) = (0.2, 0.6);
        let expected = 0.5 * (p0 + p1); // 0.4
        let trials = 1000u64;
        let mut lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(1);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.h(0);
            t.asymmetric_loss_channel(0, p0, p1);
            if t.is_lost[0] {
                lost += 1;
            }
        }
        let frac = lost as f64 / trials as f64;
        assert!(
            (frac - expected).abs() < 0.07,
            "expected ~{expected}, got {frac:.3}"
        );
    }
}
