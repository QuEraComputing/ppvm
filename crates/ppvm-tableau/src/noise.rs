// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use std::debug_assert;
use std::fmt::Debug;

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
        self.is_lost_or_leaked(addr)
    }
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
    /// * `p[1]`: The probability of losing either one qubit when both of them are
    ///   in the qubit subspace.
    /// * `p[2]`: The probability of losing one qubit when the other one has already
    ///   been lost prior to the channel.
    fn correlated_loss_channel(
        &mut self,
        addr0: usize,
        addr1: usize,
        p: [<T as Config>::Coeff; 3],
    ) {
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
            cumulative += p_i.clone();
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

impl<T: Config, I: TableauIndex, C: SparseVector<Complex<T::Coeff>, I>> LeakageChannel<T>
    for GeneralizedTableau<T, I, C>
where
    <<T as Config>::Storage as BitView>::Store: PrimInt,
    C: std::fmt::Debug,
    T::Coeff: PartialOrd
        + PartialOrd<f64>
        + One
        + Zero
        + Clone
        + num::Num
        + ToPrimitive
        + std::fmt::Debug,
    Complex<T::Coeff>: std::ops::Mul<Output = Complex<T::Coeff>>
        + From<Complex64>
        + std::ops::MulAssign
        + std::ops::AddAssign
        + One
        + ComplexFloat
        + Copy,
    I: Debug,
{
    fn leakage_channel(&mut self, addr0: usize, p0: T::Coeff, p1: T::Coeff) {
        if self.is_lost_or_leaked(addr0) {
            return;
        }

        debug_assert!(T::Coeff::zero() <= p0 && p0 <= T::Coeff::one());
        debug_assert!(T::Coeff::zero() <= p1 && p1 <= T::Coeff::one());
        debug_assert!(
            T::Coeff::zero() <= p0.clone() + p1.clone()
                && p0.clone() + p1.clone() <= T::Coeff::one()
        );

        let p_tot = p0.clone() + p1;
        let r = self.tableau.rng.random::<f64>();

        if p_tot <= r {
            return;
        }

        // Collapse the qubit to a definite basis state. This internal
        // measurement is a mechanism, not a logical measurement, so drop the
        // record entry it pushed (mirrors `loss_channel`).
        let m = self
            .measure(addr0)
            .expect("Loss was checked before, this should be unreachable");
        self.measurement_record.pop();

        // Pin the qubit to |0⟩ (prob p0) or |1⟩ (prob p1). r < p_tot = p0 + p1
        // here, so r < p0 selects |0⟩ and p0 <= r < p_tot selects |1⟩. The pin
        // must be applied before flagging the qubit leaked, otherwise the `x`
        // gate would be skipped by `is_lost_or_leaked`.
        if p0 > r {
            if m {
                self.x(addr0);
            }
        } else if !m {
            self.x(addr0);
        }
        self.is_leaked[addr0] = true;
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
        // p[1]=1 → exactly one qubit lost each time.
        let trials = 200;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]);
            assert!(
                t.is_lost[0] ^ t.is_lost[1],
                "Expected exactly one lost qubit (seed {seed})"
            );
        }
    }

    #[test]
    fn correlated_loss_p1_both_qubits_chosen_equally() {
        // With p[1]=1 the coin flip should lose addr0 and addr1 with equal frequency.
        let trials = 1000u64;
        let mut addr0_lost = 0u64;
        for seed in 0..trials {
            let mut t = tab(2);
            t.tableau.rng = rand::SeedableRng::seed_from_u64(seed);
            t.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]);
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
            t.correlated_loss_channel(0, 1, [0.0, 1.0, 0.0]);
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
        // P(exactly one lost) should converge to p[1].
        let p_single = 0.4_f64;
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
            (fraction - p_single).abs() < 0.08,
            "Expected ~{p_single:.2}, got {fraction:.3}"
        );
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

    // === LeakageChannel ===

    #[test]
    fn leakage_p0_p1_zero_no_leak() {
        // p0 = p1 = 0 → p_tot = 0, never leaks; the qubit stays live.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 0.0);
        assert!(!t.is_leaked[0]);
        assert!(!t.is_lost[0]);
        assert!(!t.measure(0).unwrap());
    }

    #[test]
    fn leakage_to_zero_pins_qubit_to_zero() {
        // Start in |1⟩; leak-to-|0⟩ (p0 = 1) must pin the qubit to |0⟩.
        let mut t = tab(1);
        t.x(0);
        t.leakage_channel(0, 1.0, 0.0);
        assert!(t.is_leaked[0]);
        assert!(!t.is_lost[0]); // leaked, not lost
        assert_eq!(t.measure(0), Some(false));
    }

    #[test]
    fn leakage_to_one_pins_qubit_to_one() {
        // Start in |0⟩; leak-to-|1⟩ (p1 = 1) must pin the qubit to |1⟩.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0);
        assert!(t.is_leaked[0]);
        assert!(!t.is_lost[0]);
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn leaked_qubit_reports_a_bit_unlike_lost() {
        // A leaked qubit measures a definite bit; a lost qubit returns None.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0);
        assert!(t.measure(0).is_some());
    }

    #[test]
    fn leakage_does_not_pollute_measurement_record() {
        // The internal collapse is a mechanism, not a logical measurement;
        // mirrors `loss_channel`.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0);
        assert!(t.current_measurement_record().is_empty());
    }

    #[test]
    fn leakage_collapses_superposition_then_pins() {
        // Leaking a superposed qubit collapses and pins it, so later
        // measurement is deterministic even though |+⟩ alone would be random.
        let mut t = tab(1);
        t.tableau.rng = rand::SeedableRng::seed_from_u64(7);
        t.h(0); // |+⟩
        t.leakage_channel(0, 0.0, 1.0);
        assert!(t.is_leaked[0]);
        assert_eq!(t.measure(0), Some(true));
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn single_qubit_gate_skips_leaked_qubit() {
        // Pinned to |1⟩; a subsequent x must be a no-op.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0);
        t.x(0);
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn two_qubit_gate_skipped_when_control_leaked() {
        // Control leaked to |1⟩; cnot must not flip the (live) target.
        let mut t = tab(2);
        t.leakage_channel(0, 0.0, 1.0);
        t.cnot(0, 1);
        assert_eq!(t.measure(1), Some(false));
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn leaked_qubit_stays_deterministic_after_other_ops() {
        // A leaked qubit is disentangled and pinned: gating/measuring other
        // qubits doesn't disturb its outcome.
        let mut t = tab(2);
        t.leakage_channel(0, 0.0, 1.0);
        t.h(1);
        let _ = t.measure(1);
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn leakage_channel_skips_already_leaked() {
        // A second leakage on a leaked qubit is a no-op (early return), so the
        // pinned value is unchanged.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0); // |1⟩, leaked
        t.leakage_channel(0, 1.0, 0.0); // would pin |0⟩ if it ran
        assert_eq!(t.measure(0), Some(true));
    }

    #[test]
    fn leaked_qubit_can_still_be_lost() {
        // Leaked-then-lost is allowed; loss wins and measurement returns None.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0);
        t.loss_channel(0, 1.0);
        assert!(t.is_lost[0]);
        assert!(t.measure(0).is_none());
    }

    #[test]
    fn reset_skips_leaked_qubit() {
        // reset must not re-zero a leaked qubit.
        let mut t = tab(1);
        t.leakage_channel(0, 0.0, 1.0); // leaked, |1⟩
        t.reset(0);
        assert!(t.is_leaked[0]);
        assert_eq!(t.measure(0), Some(true));
    }
}
