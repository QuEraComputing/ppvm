// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_tableau_2::prelude::{
    Bitstring, Clifford, CliffordBatch, CliffordExtensions, CliffordExtensionsBatch,
    CorrelatedLossChannel, Depolarizing, Depolarizing2, GeneralizedTableau, LossChannel, Measure,
    PauliError, Reset, RotationOne, TGate, TwoQubitPauliError, U3Gate,
};

use super::StimTableau;
use crate::executor::helpers::{has_repeats, measure_reset_z};
use rand::RngExt;

macro_rules! unary {
    ($name:ident, $trait:ident) => {
        fn $name(&mut self, q: usize) {
            $trait::$name(self, q);
        }
    };
}
macro_rules! binary {
    ($name:ident, $trait:ident) => {
        fn $name(&mut self, a: usize, b: usize) {
            $trait::$name(self, a, b);
        }
    };
}
macro_rules! batch {
    ($name:ident, $trait:ident, $ty:ty) => {
        fn $name(&mut self, q: &[$ty]) {
            $trait::$name(self, q);
        }
    };
}

/// Flip each of the last `outcomes.len()` records with probability `noise`,
/// pushing the recorded bits onto `results`. A lost qubit's `None` stays.
fn record_with_noise<I: Bitstring, H, R: rand::Rng + ?Sized>(
    tab: &mut GeneralizedTableau<I, H>,
    outcomes: &[Option<bool>],
    noise: f64,
    rng: &mut R,
    results: &mut Vec<Option<bool>>,
) {
    let base = tab.measurement_record.len() - outcomes.len();
    for (k, outcome) in outcomes.iter().enumerate() {
        let recorded = outcome.map(|b| GeneralizedTableau::<I, H>::flip_with_prob(b, noise, rng));
        tab.measurement_record[base + k] = recorded;
        results.push(recorded);
    }
}

/// Call `hit(i, rng)` for each index in `0..n` that an independent event of
/// probability `p` lands on, in order. Gaps between hits are geometric, so this
/// draws once per hit rather than once per index — Stim's `RareErrorIterator`.
fn for_each_hit<R: rand::Rng + ?Sized>(
    n: usize,
    p: f64,
    rng: &mut R,
    mut hit: impl FnMut(usize, &mut R),
) {
    if p <= 0.0 {
        return;
    }
    if p >= 1.0 {
        (0..n).for_each(|i| hit(i, rng));
        return;
    }
    let log_miss = (-p).ln_1p();
    let mut i = 0;
    while i < n {
        // `1 - u` is in (0, 1], so the gap is finite and non-negative.
        let gap = (1.0 - rng.random::<f64>()).ln() / log_miss;
        if gap >= (n - i) as f64 {
            return;
        }
        i += gap as usize;
        hit(i, rng);
        i += 1;
    }
}

/// Pauli `1 = X`, `2 = Y`, `3 = Z` on `q`; `0` is the identity. A lost qubit is
/// left alone by the gates themselves.
fn apply_pauli<I: Bitstring, H>(tab: &mut GeneralizedTableau<I, H>, q: usize, pauli: usize) {
    match pauli {
        1 => Clifford::x(tab, q),
        2 => Clifford::y(tab, q),
        3 => Clifford::z(tab, q),
        _ => {}
    }
}

/// `X` on every target whose true outcome was `1`.
fn reset_ones<I: Bitstring, H>(
    tab: &mut GeneralizedTableau<I, H>,
    q: &[usize],
    outcomes: &[Option<bool>],
) {
    for (&q, &outcome) in q.iter().zip(outcomes) {
        if outcome == Some(true) {
            Clifford::x(tab, q);
        }
    }
}

impl<I, H> StimTableau for GeneralizedTableau<I, H>
where
    I: Bitstring,
{
    // The `-2` frame is runtime-sized, so it has no storage-width config to
    // name. The selector's `C` slot is therefore `()` for this backend; see
    // `adapter/mod.rs`.
    type Config = ();
    type Index = I;
    type Coefficients = H;

    fn reset<R: rand::Rng + ?Sized>(&mut self, q: usize, rng: &mut R) {
        Reset::reset(self, q, rng);
    }
    unary!(x, Clifford);
    unary!(y, Clifford);
    unary!(z, Clifford);
    unary!(h, Clifford);
    unary!(s, Clifford);
    unary!(s_dag, CliffordExtensions);
    unary!(sqrt_x, CliffordExtensions);
    unary!(sqrt_x_dag, CliffordExtensions);
    unary!(sqrt_y, CliffordExtensions);
    unary!(sqrt_y_dag, CliffordExtensions);
    binary!(cnot, Clifford);
    binary!(cy, CliffordExtensions);
    binary!(cz, Clifford);
    batch!(x_many, CliffordBatch, usize);
    batch!(y_many, CliffordBatch, usize);
    batch!(z_many, CliffordBatch, usize);
    batch!(h_many, CliffordBatch, usize);
    batch!(s_many, CliffordBatch, usize);
    batch!(s_dag_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_x_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_x_dag_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_y_many, CliffordExtensionsBatch, usize);
    batch!(sqrt_y_dag_many, CliffordExtensionsBatch, usize);
    batch!(cnot_many, CliffordBatch, (usize, usize));
    batch!(cy_many, CliffordExtensionsBatch, (usize, usize));
    batch!(cz_many, CliffordBatch, (usize, usize));
    unary!(t, TGate);
    unary!(t_dag, TGate);

    fn rx(&mut self, q: usize, theta: f64) {
        RotationOne::rx(self, q, theta);
    }
    fn ry(&mut self, q: usize, theta: f64) {
        RotationOne::ry(self, q, theta);
    }
    fn rz(&mut self, q: usize, theta: f64) {
        RotationOne::rz(self, q, theta);
    }
    fn u3(&mut self, q: usize, theta: f64, phi: f64, lambda: f64) {
        U3Gate::u3(self, q, theta, phi, lambda);
    }
    fn depolarize1<R: rand::Rng + ?Sized>(&mut self, q: usize, p: f64, rng: &mut R) {
        Depolarizing::depolarize1(self, q, p, rng);
    }
    fn depolarize2<R: rand::Rng + ?Sized>(&mut self, a: usize, b: usize, p: f64, rng: &mut R) {
        Depolarizing2::depolarize2(self, a, b, p, rng);
    }
    fn pauli_error<R: rand::Rng + ?Sized>(&mut self, q: usize, p: [f64; 3], rng: &mut R) {
        PauliError::pauli_error(self, q, p, rng);
    }
    fn two_qubit_pauli_error<R: rand::Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 15],
        rng: &mut R,
    ) {
        TwoQubitPauliError::two_qubit_pauli_error(self, a, b, p, rng);
    }
    fn loss_channel<R: rand::Rng + ?Sized>(&mut self, q: usize, p: f64, rng: &mut R) {
        LossChannel::loss_channel(self, q, p, rng);
    }
    fn correlated_loss_channel<R: rand::Rng + ?Sized>(
        &mut self,
        a: usize,
        b: usize,
        p: [f64; 3],
        rng: &mut R,
    ) {
        CorrelatedLossChannel::correlated_loss_channel(self, a, b, p, rng);
    }
    fn measure<R: rand::Rng + ?Sized>(&mut self, q: usize, rng: &mut R) -> Option<bool> {
        Measure::measure(self, q, rng)
    }
    fn measure_many<R: rand::Rng + ?Sized>(
        &mut self,
        q: &[usize],
        rng: &mut R,
    ) -> Vec<Option<bool>> {
        GeneralizedTableau::measure_batch(self, q, rng)
    }
    fn measure_noisy<R: rand::Rng + ?Sized>(
        &mut self,
        q: usize,
        p: f64,
        rng: &mut R,
    ) -> Option<bool> {
        GeneralizedTableau::measure_noisy(self, q, p, rng)
    }
    fn flip_with_prob<R: rand::Rng + ?Sized>(&mut self, bit: bool, p: f64, rng: &mut R) -> bool {
        GeneralizedTableau::<I, H>::flip_with_prob(bit, p, rng)
    }
    // Batch the measurements, then apply the `X` resets: `X_q` commutes with
    // `Z_p` for `p != q`, so deferring it is exact. A repeated target can't defer.
    fn reset_many<R: rand::Rng + ?Sized>(&mut self, q: &[usize], rng: &mut R) {
        if has_repeats(q) {
            q.iter().for_each(|&q| Reset::reset(self, q, rng));
            return;
        }
        let outcomes = self.measure_batch(q, rng);
        let kept = self.measurement_record.len() - q.len();
        self.measurement_record.truncate(kept);
        reset_ones(self, q, &outcomes);
    }
    fn measure_reset_many<R: rand::Rng + ?Sized>(
        &mut self,
        q: &[usize],
        noise: f64,
        rng: &mut R,
        results: &mut Vec<Option<bool>>,
    ) {
        if has_repeats(q) {
            q.iter()
                .for_each(|&q| results.push(measure_reset_z(self, q, noise, rng)));
            return;
        }
        let outcomes = self.measure_batch(q, rng);
        record_with_noise(self, &outcomes, noise, rng, results);
        reset_ones(self, q, &outcomes);
    }
    // A repeated target needs no fallback: its second measurement is simply
    // deterministic, and each record's flip is independent.
    fn measure_noisy_many<R: rand::Rng + ?Sized>(
        &mut self,
        q: &[usize],
        noise: f64,
        rng: &mut R,
        results: &mut Vec<Option<bool>>,
    ) {
        let outcomes = self.measure_batch(q, rng);
        record_with_noise(self, &outcomes, noise, rng, results);
    }

    // Noise draws once per error (geometric gaps), not once per target; each
    // error then picks its Pauli. Same distribution as the per-target channels.
    fn depolarize1_many<R: rand::Rng + ?Sized>(&mut self, q: &[usize], p: f64, rng: &mut R) {
        for_each_hit(q.len(), p, rng, |i, rng| {
            apply_pauli(self, q[i], rng.random_range(1..4));
        });
    }
    fn depolarize2_many<R: rand::Rng + ?Sized>(&mut self, q: &[usize], p: f64, rng: &mut R) {
        let pairs = q.len() / 2;
        for_each_hit(pairs, p, rng, |i, rng| {
            let (a, b) = (q[2 * i], q[2 * i + 1]);
            if self.is_lost[a] || self.is_lost[b] {
                return;
            }
            // One of the 15 non-identity pairs `(k / 4, k % 4)`, as `depolarize2`.
            let k = rng.random_range(1..16);
            apply_pauli(self, a, k / 4);
            apply_pauli(self, b, k % 4);
        });
    }
    fn pauli_error_many<R: rand::Rng + ?Sized>(&mut self, q: &[usize], p: [f64; 3], rng: &mut R) {
        let total: f64 = p.iter().sum();
        for_each_hit(q.len(), total, rng, |i, rng| {
            let r = rng.random::<f64>() * total;
            let pauli = if r < p[0] {
                1
            } else if r < p[0] + p[1] {
                2
            } else {
                3
            };
            apply_pauli(self, q[i], pauli);
        });
    }

    fn measurement_record(&self) -> &[Option<bool>] {
        self.current_measurement_record()
    }
    fn append_measurement_record(&mut self, result: Option<bool>) {
        GeneralizedTableau::append_measurement_record(self, result);
    }
    fn overwrite_last_measurement_record(&mut self, result: Option<bool>) {
        GeneralizedTableau::overwrite_last_measurement_record(self, result);
    }
}
