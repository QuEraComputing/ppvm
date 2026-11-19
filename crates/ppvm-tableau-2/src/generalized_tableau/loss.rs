// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Loss and leakage. A lost qubit is collapsed, reset to `|0⟩` and skipped by
//! gates; measuring it gives `None`. A leaked qubit is collapsed and pinned to
//! `|0⟩` or `|1⟩`; gates skip it, and measuring it gives the pinned bit. Loss
//! overwrites leakage. Neither event is a measurement outcome.

use num_complex::Complex64;
use rand::RngExt;

use ppvm_traits_2::{
    AsymmetricLossChannel, Clifford, CorrelatedLossChannel, LeakageChannel, LossChannel,
    ResetLeakageChannel, ResetLossChannel, ResetState,
};

use super::GeneralizedTableau;
use super::status::QubitStatus;
use crate::VirtualState;

impl<V: VirtualState> GeneralizedTableau<V> {
    /// Reset `q`, collapsed onto outcome `m`, to `|0⟩` and mark it lost. The `X`
    /// goes straight to the tableau, as `q` may be leaked.
    fn mark_lost(&mut self, q: usize, m: bool) {
        if m {
            self.tableau.x(q);
        }
        self.status[q] = QubitStatus::Lost;
    }

    /// Collapse `q` by sampling `Z_q`, then reset it and mark it lost.
    fn lose(&mut self, q: usize) {
        let m = self.sample_z(q);
        self.mark_lost(q, m);
    }
}

/// With probability `p`, lose `q`. A lost qubit stays lost.
impl<V: VirtualState> LossChannel<f64> for GeneralizedTableau<V> {
    fn loss_channel(&mut self, q: usize, p: f64) {
        if !self.is_lost(q) && self.rng.random::<f64>() < p {
            self.lose(q);
        }
    }
}

/// Bring a lost qubit back as live; it is already in `|0⟩`.
impl<V: VirtualState> ResetLossChannel for GeneralizedTableau<V> {
    fn reset_loss_channel(&mut self, q: usize) {
        if self.is_lost(q) {
            self.status[q] = QubitStatus::Live;
        }
    }
}

/// With probability `p0` pin `q` to `|0⟩` and mark it leaked, with `p1` to
/// `|1⟩`. Lost and leaked qubits are left alone.
impl<V: VirtualState> LeakageChannel<f64> for GeneralizedTableau<V> {
    fn leakage_channel(&mut self, q: usize, p0: f64, p1: f64) {
        debug_assert!(p0 >= 0.0 && p1 >= 0.0 && p0 + p1 <= 1.0 + 1e-9);
        if self.is_inactive(q) {
            return;
        }
        let r = self.rng.random::<f64>();
        if r >= p0 + p1 {
            return;
        }
        let pin = r >= p0;
        if self.sample_z(q) != pin {
            self.tableau.x(q);
        }
        self.status[q] = QubitStatus::Leaked;
    }
}

/// Bring a leaked qubit back as live `|0⟩`. Live and lost qubits are untouched.
impl<V: VirtualState> ResetLeakageChannel for GeneralizedTableau<V> {
    fn reset_leakage_channel(&mut self, q: usize) {
        if self.is_leaked(q) {
            self.status[q] = QubitStatus::Live;
            self.reset_state(q);
        }
    }
}

/// `p = [p_LL, p_LQ, p_LN]` as the trait documents: if one of the pair is
/// already lost, the other is lost with `p[2]`; otherwise both with `p[0]` and
/// each named one alone with `p[1]`. A leaked partner counts as present.
impl<V: VirtualState> CorrelatedLossChannel<f64> for GeneralizedTableau<V> {
    fn correlated_loss_channel(&mut self, a: usize, b: usize, p: [f64; 3]) {
        debug_assert!(
            p[0] >= 0.0 && p[1] >= 0.0 && p[0] + 2.0 * p[1] <= 1.0 + 1e-9,
            "correlated loss needs p0, p1 >= 0 and p0 + 2·p1 <= 1; got {p:?}"
        );
        debug_assert!((0.0..=1.0).contains(&p[2]), "p2 in [0, 1]; got {p:?}");
        if self.is_lost(a) {
            self.loss_channel(b, p[2]);
        } else if self.is_lost(b) {
            self.loss_channel(a, p[2]);
        } else {
            match self.draw(&[p[0], p[1], p[1]]) {
                Some(0) => {
                    self.lose(a);
                    self.lose(b);
                }
                Some(1) => self.lose(a),
                Some(_) => self.lose(b),
                None => {}
            }
        }
    }
}

/// State-dependent loss, with its back-action (issue #249). With Born weights
/// `P(m)` of `Z_q`, `q` is lost from `|m⟩` with probability `p_m·P(m)`;
/// otherwise the survival Kraus operator
/// `K₀ = √(1−p0)|0⟩⟨0| + √(1−p1)|1⟩⟨1| = a·I + b·Z_q` is applied and the state
/// renormalized. On a qubit in superposition `K₀` is non-Clifford and
/// activates an axis.
impl<V: VirtualState> AsymmetricLossChannel<f64> for GeneralizedTableau<V> {
    fn asymmetric_loss_channel(&mut self, q: usize, p0: f64, p1: f64) {
        debug_assert!((0.0..=1.0).contains(&p0) && (0.0..=1.0).contains(&p1));
        if self.is_lost(q) {
            return;
        }
        let [w0, w1] = self.z_weights(q);
        let pop = [w0 / (w0 + w1), w1 / (w0 + w1)];
        if let Some(k) = self.draw(&[p0 * pop[0], p1 * pop[1]]) {
            // Lost from |m⟩: the outcome comes from the branch, not from `P(m)`.
            let m = k == 1;
            self.collapse(q, m);
            self.mark_lost(q, m);
            return;
        }
        // K₀ ∝ I: the renormalized state is unchanged, and rotating would
        // activate an axis for nothing.
        if p0 == p1 {
            return;
        }
        // Scaled so the norm is unchanged: ‖K₀ψ‖² = Σ (1 − p_m)·P(m).
        let scale = ((1.0 - p0) * pop[0] + (1.0 - p1) * pop[1]).sqrt().recip();
        let (s0, s1) = ((1.0 - p0).sqrt(), (1.0 - p1).sqrt());
        let a = Complex64::from((s0 + s1) / 2.0 * scale);
        let b = Complex64::from((s0 - s1) / 2.0 * scale);
        self.rotate_z(&[q], a, b);
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use ppvm_traits_2::{Measure, Pauli, Reset, RotationOne, RotationTwo, TGate};

    use super::*;
    use crate::{ActiveDormant, DenseState};

    type Sim = GeneralizedTableau<ActiveDormant<DenseState>>;

    /// Whether `count` of `shots` is within 5σ of `rate`.
    fn within(count: usize, shots: usize, rate: f64) -> bool {
        let sigma = (rate * (1.0 - rate) / shots as f64).sqrt();
        ((count as f64 / shots as f64) - rate).abs() <= 5.0 * sigma
    }

    /// A lost qubit reads `None` and no gate reaches it, including a CX from a
    /// live `|1⟩` control. Clearing the flag finds it still in `|0⟩`.
    #[test]
    fn lost_qubit_is_skipped_and_reads_none() {
        let mut s = Sim::new_with_seed(2, 0);
        s.loss_channel(0, 1.0);
        assert!(s.is_lost(0));
        s.x(1);
        s.cnot(1, 0);
        s.x(0);
        s.h(0);
        s.rx(0, 1.0);
        s.t(0);
        assert_eq!(s.measure(0), None);
        assert_eq!(s.measure(1), Some(true));
        s.reset_loss_channel(0);
        assert!(!s.is_lost(0));
        assert_eq!(s.measure(0), Some(false));
    }

    /// `exp(−iπ/2 X⊗X)` with the second qubit lost acts as `exp(−iπ/2 X)` on
    /// the first and flips it; skipping it whole would leave `|0⟩`.
    #[test]
    fn rotation_with_a_lost_partner_acts_on_the_survivor() {
        let mut s = Sim::new_with_seed(2, 0);
        s.loss_channel(1, 1.0);
        s.rotate_2(Pauli::X, Pauli::X, 0, 1, PI);
        assert_eq!(s.measure(0), Some(true));
        assert_eq!(s.measure(1), None);
    }

    /// Loss is an unrecorded measurement: losing half of a Bell pair leaves the
    /// partner in a Z eigenstate, `|1⟩` half of the time.
    #[test]
    fn loss_collapses_the_partner() {
        let shots = 2000;
        let mut ones = 0;
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(2, seed as u64);
            s.h(0);
            s.cnot(0, 1);
            s.loss_channel(0, 1.0);
            let [w0, w1] = s.z_weights(1);
            assert!(w0 == 0.0 || w1 == 0.0, "partner is deterministic");
            ones += usize::from(s.measure(1).unwrap());
        }
        assert!(within(ones, shots, 0.5), "{ones} of {shots}");
    }

    /// `loss_channel(p)` loses a qubit with probability `p`, whatever its state.
    #[test]
    fn loss_rate() {
        let (p, shots) = (0.3, 20_000);
        let mut lost = 0;
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(1, seed as u64);
            s.h(0);
            s.loss_channel(0, p);
            lost += usize::from(s.is_lost(0));
        }
        assert!(within(lost, shots, p), "{lost} of {shots}");
    }

    /// A full reset brings a lost qubit back as live `|0⟩`, so gates reach it.
    #[test]
    fn reset_revives_a_lost_qubit() {
        let mut s = Sim::new_with_seed(1, 0);
        s.loss_channel(0, 1.0);
        s.reset(0);
        assert!(!s.is_lost(0));
        assert_eq!(s.measure(0), Some(false));
        s.x(0);
        assert_eq!(s.measure(0), Some(true));
    }

    /// `leakage_channel(p0, p1)` pins to `|0⟩` with `p0` and `|1⟩` with `p1`,
    /// whatever the state. Gates skip a leaked qubit, so it keeps reading the
    /// pinned bit.
    #[test]
    fn leakage_pins_with_the_branch_probabilities() {
        let (p0, p1, shots) = (0.2, 0.3, 20_000);
        let (mut zeros, mut ones) = (0, 0);
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(1, seed as u64);
            s.h(0);
            s.leakage_channel(0, p0, p1);
            if s.is_leaked(0) {
                let m = s.measure(0).unwrap();
                s.x(0);
                s.h(0);
                s.rx(0, 1.0);
                assert_eq!(s.measure(0), Some(m), "gates skip a leaked qubit");
                if m {
                    ones += 1;
                } else {
                    zeros += 1;
                }
            }
        }
        assert!(within(zeros, shots, p0), "{zeros} of {shots}");
        assert!(within(ones, shots, p1), "{ones} of {shots}");
    }

    /// Leakage reduction returns a leaked qubit as live `|0⟩`; it leaves a
    /// lost qubit lost.
    #[test]
    fn reset_leakage_revives_only_leaked_qubits() {
        let mut s = Sim::new_with_seed(2, 0);
        s.leakage_channel(0, 0.0, 1.0);
        assert_eq!(s.measure(0), Some(true));
        s.reset_leakage_channel(0);
        assert!(!s.is_leaked(0));
        assert_eq!(s.measure(0), Some(false));
        s.x(0);
        assert_eq!(s.measure(0), Some(true));
        s.loss_channel(1, 1.0);
        s.reset_leakage_channel(1);
        assert!(s.is_lost(1));
    }

    /// Loss overwrites leakage: the pinned `|1⟩` is reset to `|0⟩`.
    #[test]
    fn loss_overwrites_leakage() {
        let mut s = Sim::new_with_seed(1, 0);
        s.leakage_channel(0, 0.0, 1.0);
        s.loss_channel(0, 1.0);
        assert!(s.is_lost(0) && !s.is_leaked(0));
        assert_eq!(s.measure(0), None);
        s.reset_loss_channel(0);
        assert_eq!(s.measure(0), Some(false));
    }

    /// A full reset brings a leaked qubit back as live `|0⟩`.
    #[test]
    fn reset_revives_a_leaked_qubit() {
        let mut s = Sim::new_with_seed(1, 0);
        s.leakage_channel(0, 0.0, 1.0);
        s.reset(0);
        assert!(!s.is_leaked(0));
        assert_eq!(s.measure(0), Some(false));
    }

    /// `p = [p_LL, p_LQ, p_LN]`: from two present qubits both are lost with
    /// `p_LL`, each named one alone with `p_LQ`; with one already lost, the
    /// other is lost with `p_LN`.
    #[test]
    fn correlated_loss_event_rates() {
        let (p, shots) = ([0.1, 0.2, 0.5], 20_000);
        let mut counts = [0; 4];
        let mut second = 0;
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(2, seed as u64);
            s.h(0);
            s.cnot(0, 1);
            s.correlated_loss_channel(0, 1, p);
            counts[usize::from(s.is_lost(0)) | usize::from(s.is_lost(1)) << 1] += 1;
            let mut t = Sim::new_with_seed(2, seed as u64);
            t.loss_channel(0, 1.0);
            t.correlated_loss_channel(0, 1, p);
            second += usize::from(t.is_lost(1));
        }
        let rates = [1.0 - p[0] - 2.0 * p[1], p[1], p[1], p[0]];
        for (k, (&count, rate)) in counts.iter().zip(rates).enumerate() {
            assert!(
                within(count, shots, rate),
                "event {k:#b}: {count} of {shots}"
            );
        }
        assert!(within(second, shots, p[2]), "{second} of {shots}");
    }

    /// The Bell-pair example of #249 with `p0 = 0, p1 = 1`: the qubit is lost
    /// exactly from `|1⟩`, so a loss leaves the partner in `|1⟩`; survival
    /// projects onto `|0⟩`, leaving `|00⟩`. Each happens half of the time.
    #[test]
    fn asymmetric_loss_on_a_bell_pair() {
        let shots = 4000;
        let mut lost = 0;
        for seed in 0..shots {
            let mut s = Sim::new_with_seed(2, seed as u64);
            s.h(0);
            s.cnot(0, 1);
            s.asymmetric_loss_channel(0, 0.0, 1.0);
            if s.is_lost(0) {
                lost += 1;
                assert_eq!(s.measure(1), Some(true), "lost from |1⟩");
            } else {
                assert_eq!((s.measure(0), s.measure(1)), (Some(false), Some(false)));
            }
        }
        assert!(within(lost, shots, 0.5), "{lost} of {shots}");
    }

    /// With `p0 = p1` the survival operator is proportional to `I`, so a
    /// surviving qubit in superposition activates no axis.
    #[test]
    fn symmetric_asymmetric_loss_activates_no_axis() {
        for seed in 0..50 {
            let mut s = Sim::new_with_seed(2, seed);
            s.h(0);
            s.cnot(0, 1);
            s.asymmetric_loss_channel(0, 0.3, 0.3);
            assert!(s.state().active().is_empty(), "seed {seed}");
        }
    }
}
