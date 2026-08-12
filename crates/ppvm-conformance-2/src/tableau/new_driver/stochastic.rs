// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::{Bitstring, GeneralizedTableau};
use ppvm_traits_2::{
    AsymmetricLossChannel, Clifford, CorrelatedLossChannel, Depolarizing, Depolarizing2,
    LossChannel, Measure, PauliError, Reset, ResetLossChannel, TwoQubitPauliError,
};
use rand::{Rng, RngExt};

use super::super::NewDriver;

type Driven<I, H> = NewDriver<GeneralizedTableau<I, H>>;

impl<I: Bitstring, H> Driven<I, H> {
    /// Preserve the old generalized-tableau stream in differential workloads.
    ///
    /// The migrated core intentionally skips deterministic probability draws.
    /// The reference generalized tableau did not, so this adapter samples from a
    /// cloned stream and advances its owned compatibility stream exactly once
    /// for every present measurement.
    pub fn measure(&mut self, qubit: usize) -> Option<bool> {
        if self.tab.is_lost[qubit] {
            return self.tab.measure(qubit, &mut self.rng);
        }
        let mut sampled = self.rng.clone();
        let outcome = self.tab.measure(qubit, &mut sampled);
        let _: f64 = self.rng.random();
        outcome
    }

    pub fn measure_many(&mut self, qubits: &[usize]) -> Vec<Option<bool>> {
        qubits.iter().map(|&q| self.measure(q)).collect()
    }

    pub fn measure_all(&mut self) -> Vec<Option<bool>> {
        (0..self.tab.n_qubits()).map(|q| self.measure(q)).collect()
    }

    pub fn measure_many_with_scratch(
        &mut self,
        qubits: &[usize],
        scratch: &mut ppvm_tableau_2::MeasureScratch<I>,
    ) -> Vec<Option<bool>> {
        self.tab
            .measure_many_with_scratch(qubits, scratch, &mut self.rng)
    }

    pub fn measure_all_with_scratch(
        &mut self,
        scratch: &mut ppvm_tableau_2::MeasureScratch<I>,
    ) -> Vec<Option<bool>> {
        self.tab.measure_all_with_scratch(scratch, &mut self.rng)
    }

    pub fn measure_noisy(&mut self, qubit: usize, p: f64) -> Option<bool> {
        let outcome = self.measure(qubit)?;
        let noisy = if p > 0.0 && self.rng.random::<f64>() < p {
            !outcome
        } else {
            outcome
        };
        self.tab.overwrite_last_measurement_record(Some(noisy));
        Some(noisy)
    }

    pub fn reset(&mut self, qubit: usize) {
        let outcome = self.measure(qubit);
        self.tab.measurement_record.pop();
        if outcome == Some(true) {
            self.tab.x(qubit);
        }
    }

    pub fn reset_z(&mut self, qubit: usize) {
        self.reset(qubit);
    }

    pub fn reset_x(&mut self, qubit: usize) {
        self.reset(qubit);
        self.tab.h(qubit);
    }

    pub fn reset_y(&mut self, qubit: usize) {
        self.reset_x(qubit);
        self.tab.s(qubit);
    }

    pub fn reset_many(&mut self, qubits: &[usize]) {
        for &q in qubits {
            self.reset(q);
        }
    }
}

impl<I: Bitstring, H> Measure for Driven<I, H> {
    fn measure<R: Rng + ?Sized>(&mut self, q: usize, _rng: &mut R) -> Option<bool> {
        Driven::measure(self, q)
    }

    fn measure_many<R: Rng + ?Sized>(
        &mut self,
        qubits: &[usize],
        _rng: &mut R,
    ) -> Vec<Option<bool>> {
        Driven::measure_many(self, qubits)
    }
}

impl<I: Bitstring, H> Reset for Driven<I, H> {
    fn reset<R: Rng + ?Sized>(&mut self, q: usize, _rng: &mut R) {
        Driven::reset(self, q);
    }
}

impl<I: Bitstring, H> PauliError<f64> for Driven<I, H> {
    fn pauli_error<R: Rng + ?Sized>(&mut self, q: usize, p: [f64; 3], _rng: &mut R) {
        self.tab.pauli_error(q, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> TwoQubitPauliError<f64> for Driven<I, H> {
    fn two_qubit_pauli_error<R: Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 15],
        _rng: &mut R,
    ) {
        self.tab.two_qubit_pauli_error(a, b, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> Depolarizing<f64> for Driven<I, H> {
    fn depolarize1<R: Rng + ?Sized>(&mut self, q: usize, p: f64, _rng: &mut R) {
        self.tab.depolarize1(q, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> Depolarizing2<f64> for Driven<I, H> {
    fn depolarize2<R: Rng + ?Sized>(&mut self, a: usize, b: usize, p: f64, _rng: &mut R) {
        self.tab.depolarize2(a, b, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> LossChannel<f64> for Driven<I, H> {
    fn loss_channel<R: Rng + ?Sized>(&mut self, q: usize, p: f64, _rng: &mut R) {
        self.tab.loss_channel(q, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> AsymmetricLossChannel<f64> for Driven<I, H> {
    fn asymmetric_loss_channel<R: Rng + ?Sized>(
        &mut self,
        q: usize,
        p0: f64,
        p1: f64,
        _rng: &mut R,
    ) {
        self.tab.asymmetric_loss_channel(q, p0, p1, &mut self.rng);
    }
}

impl<I: Bitstring, H> CorrelatedLossChannel<f64> for Driven<I, H> {
    fn correlated_loss_channel<R: Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 3],
        _rng: &mut R,
    ) {
        self.tab.correlated_loss_channel(a, b, p, &mut self.rng);
    }
}

impl<I: Bitstring, H> ResetLossChannel for Driven<I, H> {
    fn reset_loss_channel(&mut self, q: usize) {
        self.tab.reset_loss_channel(q);
    }
}
