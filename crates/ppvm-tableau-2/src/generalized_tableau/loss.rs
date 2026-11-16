// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Loss and leakage. A lost qubit is collapsed, reset to `|0⟩` and skipped by
//! gates; measuring it gives `None`. A leaked qubit is collapsed and pinned to
//! `|0⟩` or `|1⟩`; gates skip it, and measuring it gives the pinned bit. Loss
//! overwrites leakage. Neither event is a measurement outcome.

use rand::RngExt;

use ppvm_traits_2::{Clifford, LossChannel, ResetLossChannel};

use super::GeneralizedTableau;
use crate::VirtualState;

impl<V: VirtualState> GeneralizedTableau<V> {
    /// Collapse `q` by sampling `Z_q`, reset it to `|0⟩` and mark it lost. The
    /// `X` goes straight to the tableau, as `q` may be leaked.
    fn lose(&mut self, q: usize) {
        if self.sample_z(q) {
            self.tableau.x(q);
        }
        self.status.set_lost(q);
    }
}

/// With probability `p`, lose `q`. A lost qubit stays lost.
impl<V: VirtualState> LossChannel<f64> for GeneralizedTableau<V> {
    fn loss_channel(&mut self, q: usize, p: f64) {
        if !self.status.is_lost(q) && self.rng.random::<f64>() < p {
            self.lose(q);
        }
    }
}

/// Bring a lost qubit back as live; it is already in `|0⟩`.
impl<V: VirtualState> ResetLossChannel for GeneralizedTableau<V> {
    fn reset_loss_channel(&mut self, q: usize) {
        if self.status.is_lost(q) {
            self.status.set_live(q);
        }
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
}
