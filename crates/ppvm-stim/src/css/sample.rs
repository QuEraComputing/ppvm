// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

use super::CssProgram;
use super::bits::{Bits, LANES, SHOTS};
use super::compile::{Basis, Channel, Op, PairLoss};

impl CssProgram {
    /// Sample `num_shots` shots and return their measurement records, shot-major, as one byte
    /// per measurement: `0` / `1` for the outcome and `2` for a lost qubit.
    ///
    /// The result depends only on `seed`, not on how shots are split across threads.
    pub fn sample(&self, num_shots: usize, seed: u64) -> Vec<u8> {
        let mut out = vec![0; num_shots * self.num_measurements];
        self.sample_into(&mut out, seed);
        out
    }

    /// Like [`sample`](Self::sample), writing `out.len() / num_measurements()` shots into `out`.
    ///
    /// # Panics
    ///
    /// If `out.len()` is not a multiple of [`num_measurements`](Self::num_measurements).
    pub fn sample_into(&self, out: &mut [u8], seed: u64) {
        let m = self.num_measurements;
        if m == 0 {
            return;
        }
        assert!(
            out.len().is_multiple_of(m),
            "output length {} is not a multiple of {m} measurements",
            out.len()
        );
        let base = SmallRng::seed_from_u64(seed).random::<u64>();
        let run = |batch: &mut Batch, (i, chunk): (usize, &mut [u8])| {
            let mut rng = SmallRng::seed_from_u64(base.wrapping_add(i as u64));
            batch.run(self, &mut rng);
            batch.write(chunk);
        };
        #[cfg(feature = "rayon")]
        {
            use rayon::prelude::*;
            out.par_chunks_mut(SHOTS * m)
                .enumerate()
                .for_each_init(|| Batch::new(self), run);
        }
        #[cfg(not(feature = "rayon"))]
        {
            let mut batch = Batch::new(self);
            for item in out.chunks_mut(SHOTS * m).enumerate() {
                run(&mut batch, item);
            }
        }
    }
}

/// Pauli frames and loss state of `SHOTS` shots. The reference record is all zeros, so a
/// measured value is the frame bit itself.
struct Batch {
    x: Vec<Bits>,
    z: Vec<Bits>,
    /// Lost or leaked: gates skip the qubit, and its frame is unused until a reload rewrites it.
    inactive: Vec<Bits>,
    leaked: Vec<Bits>,
    /// What a leaked qubit reads.
    pin: Vec<Bits>,
    /// Measurement records, indexed `[lane * num_measurements + m]`.
    values: Vec<u64>,
    lost: Vec<u64>,
    num_measurements: usize,
    recorded: usize,
}

impl Batch {
    fn new(program: &CssProgram) -> Self {
        let n = program.num_qubits;
        let records = LANES * program.num_measurements;
        Self {
            x: vec![Bits::ZERO; n],
            z: vec![Bits::ZERO; n],
            inactive: vec![Bits::ZERO; n],
            leaked: vec![Bits::ZERO; n],
            pin: vec![Bits::ZERO; n],
            values: vec![0; records],
            lost: vec![0; records],
            num_measurements: program.num_measurements,
            recorded: 0,
        }
    }

    fn run(&mut self, program: &CssProgram, rng: &mut SmallRng) {
        for q in 0..program.num_qubits {
            self.reload(q, Basis::Z, Bits::ONES, rng);
        }
        self.recorded = 0;
        self.execute(&program.ops, rng);
    }

    fn execute(&mut self, ops: &[Op], rng: &mut SmallRng) {
        for op in ops {
            match *op {
                Op::Reset { q, basis } => self.reload(q, basis, Bits::ONES, rng),
                Op::ResetState(q) => self.prepare(q, Basis::Z, Bits::ONES, rng),
                Op::ResetLoss(q) => {
                    let lost = self.inactive[q] & !self.leaked[q];
                    self.reload(q, Basis::Z, lost, rng);
                }
                Op::ResetLeakage(q) => self.reload(q, Basis::Z, self.leaked[q], rng),
                Op::Measure {
                    q,
                    basis,
                    reset,
                    log1m_flip,
                } => {
                    self.measure(q, basis, log1m_flip, rng);
                    if reset {
                        self.reload(q, basis, Bits::ONES, rng);
                    } else if basis == Basis::Z {
                        self.z[q] ^= Bits::random(rng);
                    } else {
                        self.x[q] ^= Bits::random(rng);
                    }
                }
                Op::Pad { bit, log1m_flip } => {
                    let value = if bit { Bits::ONES } else { Bits::ZERO };
                    self.record(value ^ bernoulli(log1m_flip, rng), Bits::ZERO);
                }
                Op::H(ref qubits) => {
                    for &q in qubits {
                        std::mem::swap(&mut self.x[q], &mut self.z[q]);
                    }
                }
                Op::Pauli { q, x, z } => self.pauli(q, Bits::ONES, x, z),
                Op::Feedback { lookback, q, x, z } => {
                    let m = self.recorded - lookback;
                    let value = Bits(std::array::from_fn(|w| {
                        self.values[w * self.num_measurements + m]
                    }));
                    self.pauli(q, value, x, z);
                }
                Op::Cx(ref pairs) => {
                    for &[c, t] in pairs {
                        let keep = !(self.inactive[c] | self.inactive[t]);
                        let (xc, zt) = (self.x[c], self.z[t]);
                        self.x[t] ^= xc & keep;
                        self.z[c] ^= zt & keep;
                    }
                }
                Op::Cz(ref pairs) => {
                    for &[a, b] in pairs {
                        let keep = !(self.inactive[a] | self.inactive[b]);
                        let (xa, xb) = (self.x[a], self.x[b]);
                        self.z[a] ^= xb & keep;
                        self.z[b] ^= xa & keep;
                    }
                }
                Op::Pauli1(ref c) => self.pauli1(c, rng),
                Op::Pauli2(ref c) => self.pauli2(c, rng),
                Op::Loss(ref c) => self.loss(c, rng),
                Op::Leakage(ref c) => self.leakage(c, rng),
                Op::CorrelatedLoss(ref c) => self.correlated_loss(c, rng),
                Op::Repeat(ref r) => {
                    for _ in 0..r.count {
                        self.execute(&r.body, rng);
                    }
                }
            }
        }
    }

    /// Fresh `|0>` or `|+>` frames on the shots in `mask`.
    fn prepare(&mut self, q: usize, basis: Basis, mask: Bits, rng: &mut SmallRng) {
        let (fixed, random) = match basis {
            Basis::Z => (&mut self.x[q], &mut self.z[q]),
            Basis::X => (&mut self.z[q], &mut self.x[q]),
        };
        *fixed &= !mask;
        *random = (*random & !mask) | (Bits::random(rng) & mask);
    }

    /// Bring a fresh atom into the shots in `mask`.
    fn reload(&mut self, q: usize, basis: Basis, mask: Bits, rng: &mut SmallRng) {
        self.inactive[q] &= !mask;
        self.leaked[q] &= !mask;
        self.prepare(q, basis, mask, rng);
    }

    fn measure(&mut self, q: usize, basis: Basis, log1m_flip: f64, rng: &mut SmallRng) {
        let frame = match basis {
            Basis::Z => self.x[q],
            Basis::X => self.z[q],
        };
        let (inactive, leaked) = (self.inactive[q], self.leaked[q]);
        let value = (frame & !inactive) | (self.pin[q] & leaked);
        self.record(value ^ bernoulli(log1m_flip, rng), inactive & !leaked);
    }

    fn record(&mut self, value: Bits, lost: Bits) {
        let m = self.recorded;
        for w in 0..LANES {
            self.values[w * self.num_measurements + m] = value.0[w] & !lost.0[w];
            self.lost[w * self.num_measurements + m] = lost.0[w];
        }
        self.recorded += 1;
    }

    fn pauli(&mut self, q: usize, mask: Bits, x: bool, z: bool) {
        if x {
            self.x[q] ^= mask;
        }
        if z {
            self.z[q] ^= mask;
        }
    }

    /// Flip Pauli `code` (1 = X, 2 = Y, 3 = Z) on one shot.
    fn pauli_at(&mut self, q: usize, shot: usize, code: usize) {
        if code == 1 || code == 2 {
            self.x[q].flip(shot);
        }
        if code == 2 || code == 3 {
            self.z[q].flip(shot);
        }
    }

    fn is_lost(&self, q: usize, shot: usize) -> bool {
        self.inactive[q].get(shot) && !self.leaked[q].get(shot)
    }

    fn lose(&mut self, q: usize, shot: usize) {
        self.inactive[q].set(shot, true);
        self.leaked[q].set(shot, false);
    }

    fn pauli1(&mut self, c: &Channel<usize>, rng: &mut SmallRng) {
        events(c.log1m_p, c.targets.len(), rng, |rng, i, shot| {
            self.pauli_at(c.targets[i], shot, 1 + c.pick(rng.random()));
        });
    }

    fn pauli2(&mut self, c: &Channel<[usize; 2]>, rng: &mut SmallRng) {
        events(c.log1m_p, c.targets.len(), rng, |rng, i, shot| {
            let code = 1 + c.pick(rng.random());
            let [a, b] = c.targets[i];
            self.pauli_at(a, shot, code >> 2);
            self.pauli_at(b, shot, code & 3);
        });
    }

    fn loss(&mut self, c: &Channel<usize>, rng: &mut SmallRng) {
        events(c.log1m_p, c.targets.len(), rng, |_, i, shot| {
            self.lose(c.targets[i], shot);
        });
    }

    fn leakage(&mut self, c: &Channel<usize>, rng: &mut SmallRng) {
        events(c.log1m_p, c.targets.len(), rng, |rng, i, shot| {
            let q = c.targets[i];
            if !self.inactive[q].get(shot) {
                self.inactive[q].set(shot, true);
                self.leaked[q].set(shot, true);
                self.pin[q].set(shot, c.pick(rng.random()) == 1);
            }
        });
    }

    fn correlated_loss(&mut self, c: &PairLoss, rng: &mut SmallRng) {
        events(c.log1m_p, c.pairs.len(), rng, |rng, i, shot| {
            let [a, b] = c.pairs[i];
            let u: f64 = rng.random();
            match (self.is_lost(a, shot), self.is_lost(b, shot)) {
                (false, false) if u < c.both => {
                    self.lose(a, shot);
                    self.lose(b, shot);
                }
                (false, false) if u < c.first => self.lose(a, shot),
                (false, false) if u < c.second => self.lose(b, shot),
                (true, false) if u < c.other => self.lose(b, shot),
                (false, true) if u < c.other => self.lose(a, shot),
                _ => {}
            }
        });
    }

    /// Write the records of this batch's first `chunk.len() / num_measurements` shots, eight
    /// shots at a time through byte-wide copies of their record bits.
    fn write(&self, chunk: &mut [u8]) {
        let m = self.num_measurements;
        let (mut values, mut lost) = (vec![0u8; m], vec![0u8; m]);
        for (group, rows) in chunk.chunks_mut(8 * m).enumerate() {
            let (lane, shift) = (group / 8, group % 8 * 8);
            for k in 0..m {
                values[k] = (self.values[lane * m + k] >> shift) as u8;
                lost[k] = (self.lost[lane * m + k] >> shift) as u8;
            }
            for (bit, row) in rows.chunks_exact_mut(m).enumerate() {
                for ((out, &v), &l) in row.iter_mut().zip(&values).zip(&lost) {
                    *out = (v >> bit & 1) | (l >> bit & 1) << 1;
                }
            }
        }
    }
}

/// Call `f(rng, target, shot)` for each event of a rate-`p` process on `num_targets * SHOTS`
/// slots, skipping ahead geometrically between events.
fn events(
    log1m_p: f64,
    num_targets: usize,
    rng: &mut SmallRng,
    mut f: impl FnMut(&mut SmallRng, usize, usize),
) {
    if log1m_p >= 0.0 {
        return;
    }
    let len = num_targets * SHOTS;
    let mut pos = skip(log1m_p, rng);
    while pos < len {
        f(rng, pos / SHOTS, pos % SHOTS);
        pos = pos.saturating_add(1 + skip(log1m_p, rng));
    }
}

fn skip(log1m_p: f64, rng: &mut SmallRng) -> usize {
    let u: f64 = rng.random();
    (u.ln() / log1m_p) as usize
}

fn bernoulli(log1m_p: f64, rng: &mut SmallRng) -> Bits {
    let mut bits = Bits::ZERO;
    events(log1m_p, 1, rng, |_, _, shot| bits.flip(shot));
    bits
}
