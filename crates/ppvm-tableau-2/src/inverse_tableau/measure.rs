// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Z measurement and reset of the stabilizer state `C |0…0⟩` (no active
//! qubits).

use ppvm_pauli_kernels::{Bits, heisenberg};
use ppvm_traits_2::Clifford;

use super::bitword::{BitWord, ZERO, set_bits};
use super::transposed::TransposedTableau;
use super::{InverseTableau, QubitMask};

/// Random targets a batch needs before `measure_many` transposes. Measured at
/// 118 qubits: a batch of ~48 small collapses (the d=7 surface code) is faster
/// strided, a batch of 118 is 2.8× faster transposed. The break-even should
/// grow with `n` (transposing costs ~n²/64 words, a strided collapse ~n rows);
/// retune with other sizes.
const TRANSPOSE_THRESHOLD: usize = 64;

/// What measuring needs from a tableau layout: row `iz_q` and the collapse.
trait Collapse {
    /// Columns where row `iz_q` has an X bit, ascending, into `out`.
    fn iz_x_support(&self, q: usize, out: &mut Vec<usize>);
    fn iz_z_bit(&self, q: usize, c: usize) -> bool;
    fn iz_negative(&self, q: usize) -> bool;
    /// `C ← C·V` with `V = CX(p, targets…) · last(p)`, where `last` is
    /// [`collapse_last`] with these flags.
    fn collapse(&mut self, pivot: usize, targets: &[usize], has_y: bool, fix_sign: bool);
}

/// The collapse's final gate on the pivot: `S` if it holds a Y, then `H`, then
/// `X` if the sign calls for it. Generic so both layouts run the same rules.
#[inline(always)]
fn collapse_last<B: Bits>(has_y: bool, fix_sign: bool, x: B, z: B) -> ([B; 2], B) {
    // `x & !x` is the zero flip (`Bits` has no zero constant).
    let ([x, z], f1) = if has_y {
        heisenberg::s(x, z)
    } else {
        ([x, z], x & !x)
    };
    let ([x, z], f2) = heisenberg::h(x, z);
    let ([x, z], f3) = if fix_sign {
        heisenberg::x(x, z)
    } else {
        ([x, z], x & !x)
    };
    ([x, z], f1 ^ f2 ^ f3)
}

/// Measure `Z_q` on `t` and collapse; `targets` is scratch.
///
/// Measuring `Z_q` is measuring `O = C⁻¹ Z_q C` (row `iz_q`) on `|0…0⟩`:
///
/// - **No X in `O`:** `O` is a signed Z string, so the outcome is its sign.
/// - **X on some virtual qubits:** 50/50. Pick the lowest, `p`, as pivot
///   and change basis so that `O` becomes `±Z_p`: `CX(p, v)` clears every
///   other X (trivial on the state, since `p` is in `|0⟩`), leaving
///   `O = s · X_p · Z…` or `s · Y_p · Z…`. The projected state on `p` is
///   `|0⟩ ± |1⟩` (or `|0⟩ ± i|1⟩`), which `C ← C·H` (or `C·S·H`) and an
///   `X` when the sign calls for it map back to `|0⟩`.
fn measure_on<T: Collapse>(
    t: &mut T,
    q: usize,
    coin: impl FnOnce() -> bool,
    targets: &mut Vec<usize>,
) -> bool {
    targets.clear();
    t.iz_x_support(q, targets);
    if targets.is_empty() {
        return t.iz_negative(q);
    }
    let pivot = targets.remove(0);
    // The CXs reduce O to s · X_p or s · Y_p (times Z's); find which, and s,
    // from row iz_q alone before changing the tableau.
    let (mut pivot_bits, mut negative) = ([true, t.iz_z_bit(q, pivot)], t.iz_negative(q));
    for &v in targets.iter() {
        let ([px, pz, _, _], flip) =
            heisenberg::cnot([pivot_bits[0], pivot_bits[1], true, t.iz_z_bit(q, v)]);
        pivot_bits = [px, pz];
        negative ^= flip;
    }
    let (has_y, outcome) = (pivot_bits[1], coin());
    t.collapse(pivot, targets, has_y, outcome != negative);
    outcome
}

impl Collapse for InverseTableau {
    fn iz_x_support(&self, q: usize, out: &mut Vec<usize>) {
        out.extend(set_bits(self.iz(q).x));
    }

    fn iz_z_bit(&self, q: usize, c: usize) -> bool {
        self.iz(q).z_bit(c)
    }

    fn iz_negative(&self, q: usize) -> bool {
        self.iz(q).negative
    }

    fn collapse(&mut self, pivot: usize, targets: &[usize], has_y: bool, fix_sign: bool) {
        self.fan_out(pivot, targets, &[], |x, z| {
            collapse_last(has_y, fix_sign, x, z)
        });
    }
}

impl Collapse for TransposedTableau {
    fn iz_x_support(&self, q: usize, out: &mut Vec<usize>) {
        out.extend((0..self.num_qubits()).filter(|&c| self.bits(1, q, c)[0]));
    }

    fn iz_z_bit(&self, q: usize, c: usize) -> bool {
        self.bits(1, q, c)[1]
    }

    fn iz_negative(&self, q: usize) -> bool {
        self.negative(1, q)
    }

    fn collapse(&mut self, pivot: usize, targets: &[usize], has_y: bool, fix_sign: bool) {
        self.fan_out(pivot, targets, |x: BitWord, z: BitWord| {
            collapse_last(has_y, fix_sign, x, z)
        });
    }
}

impl InverseTableau {
    /// Measure `Z_q` on `C |0…0⟩` and collapse; returns the outcome (`true`
    /// for `|1⟩`). `coin` is called only if the outcome is random. See
    /// [`measure_on`] for the method.
    pub fn measure(&mut self, q: usize, coin: impl FnOnce() -> bool) -> bool {
        measure_on(self, q, coin, &mut Vec::new())
    }

    /// Measure `Z_q` for each of `qubits` in order, collapsing after each, and
    /// write the outcomes to `outcomes`. `coin` is called once per random
    /// outcome. Same result as calling [`measure`](Self::measure) in order.
    ///
    /// With at least `TRANSPOSE_THRESHOLD` random targets, the batch runs on a
    /// column-major copy: two transposes, then each collapse a few words per
    /// 64 rows instead of one bit per row.
    pub fn measure_many(
        &mut self,
        qubits: &[usize],
        coin: impl FnMut() -> bool,
        outcomes: &mut [bool],
    ) {
        // Count random targets only as far as the decision needs.
        let transpose = qubits.len() >= TRANSPOSE_THRESHOLD
            && qubits
                .iter()
                .filter(|&&q| self.iz(q).x.iter().any(|&w| w != ZERO))
                .nth(TRANSPOSE_THRESHOLD - 1)
                .is_some();
        self.measure_many_with(qubits, coin, outcomes, transpose);
    }

    fn measure_many_with(
        &mut self,
        qubits: &[usize],
        mut coin: impl FnMut() -> bool,
        outcomes: &mut [bool],
        transpose: bool,
    ) {
        assert_eq!(qubits.len(), outcomes.len());
        let mut targets = Vec::new();
        if transpose {
            let mut columns = TransposedTableau::from_rows(self);
            for (&q, outcome) in qubits.iter().zip(outcomes) {
                *outcome = measure_on(&mut columns, q, &mut coin, &mut targets);
            }
            columns.write_rows(self);
        } else {
            for (&q, outcome) in qubits.iter().zip(outcomes) {
                *outcome = measure_on(self, q, &mut coin, &mut targets);
            }
        }
    }

    /// Collapse `O = C⁻¹ Z_q C` onto outcome `m` on `|χ⟩ = |φ⟩ ⊗ |0…0⟩`, where
    /// `φ` lives on the `active` qubits and `O` has an X or Y on a dormant one.
    ///
    /// Stim's collapse with the pivot `p` restricted to dormant qubits: `p` is
    /// the lowest of them, `CX(p, v)` clears every other X (active ones too),
    /// and `CZ(p, c)` clears the Z on every active `c`. Both are trivial on
    /// `χ`, since `p` is in `|0⟩`, so `φ` is unchanged. What is left is
    /// `±X_p` or `±Y_p` times Z's on dormant qubits, which act as `+1`, so the
    /// final gates are those of [`measure`](Self::measure).
    pub fn collapse_dormant(&mut self, q: usize, active: &QubitMask, m: bool) {
        let o = self.iz(q);
        let pivot = o.x_outside(active).next().expect("no dormant X in O");
        let cx_targets: Vec<usize> = o.x_support().filter(|&v| v != pivot).collect();
        let cz_targets: Vec<usize> = o.z_inside(active).collect();
        // Pivot bits and sign of O after the CXs and CZs, from row iz_q alone.
        // The CXs leave the other qubits' Z bits and clear their X bits.
        let (mut pivot_bits, mut negative) = ([true, o.z_bit(pivot)], o.negative);
        for &v in &cx_targets {
            let ([px, pz, _, _], flip) =
                heisenberg::cnot([pivot_bits[0], pivot_bits[1], true, o.z_bit(v)]);
            pivot_bits = [px, pz];
            negative ^= flip;
        }
        for _ in &cz_targets {
            let ([px, pz, _, _], flip) =
                heisenberg::cz([pivot_bits[0], pivot_bits[1], false, true]);
            pivot_bits = [px, pz];
            negative ^= flip;
        }
        let has_y = pivot_bits[1];
        self.fan_out(pivot, &cx_targets, &cz_targets, |x, z| {
            collapse_last(has_y, m != negative, x, z)
        });
    }

    /// Reset qubit `q` to `|0⟩`: measure it, then flip it if it was `|1⟩`.
    pub fn reset(&mut self, q: usize, coin: impl FnOnce() -> bool) {
        if self.measure(q, coin) {
            self.x(q);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::test_support::*;
    use super::*;

    /// After a measurement, measuring again is deterministic and repeats it.
    #[test]
    fn measuring_again_repeats_the_outcome() {
        for (i, n) in SIZES.into_iter().enumerate() {
            let mut t = InverseTableau::new(n);
            random_circuit(n, 200, 0x3e5 + i as u64)
                .into_iter()
                .for_each(|g| apply(&mut t, g));
            for q in 0..n {
                let first = t.measure(q, || q % 2 == 0);
                let again = t.measure(q, || panic!("second measurement is random"));
                assert_eq!(first, again, "n = {n}, q = {q}");
            }
        }
    }

    #[test]
    fn measurement_preserves_commutation() {
        for (i, n) in [1, 2, 7, 65].into_iter().enumerate() {
            let mut t = InverseTableau::new(n);
            for (k, gate) in random_circuit(n, 100, 0x5ca1 + i as u64)
                .into_iter()
                .enumerate()
            {
                apply(&mut t, gate);
                t.measure(k % n, || k % 3 == 0);
                assert_canonical_commutation(&t);
            }
        }
    }

    /// `measure_many` equals measuring one qubit at a time with the same coins,
    /// on both the row-major and the column-major path.
    #[test]
    fn measure_many_matches_measure() {
        for (i, n) in SIZES.into_iter().enumerate() {
            let mut start = InverseTableau::new(n);
            random_circuit(n, 150, 0xba7c + i as u64)
                .into_iter()
                .for_each(|g| apply(&mut start, g));
            // Every qubit, then a repeat, so later targets see earlier collapses.
            let qubits: Vec<usize> = (0..n).rev().chain([0]).collect();
            let coins = |seed: u64| {
                let mut s = seed | 1;
                move || {
                    s ^= s << 13;
                    s ^= s >> 7;
                    s ^= s << 17;
                    s & 1 == 1
                }
            };
            let mut single = start.clone();
            let mut coin = coins(7);
            let mut expected = Vec::new();
            // `&mut coin` is needed: `coin` is reused across iterations.
            #[allow(clippy::needless_borrows_for_generic_args)]
            for &q in &qubits {
                expected.push(single.measure(q, &mut coin));
            }
            for transpose in [false, true] {
                let mut batched = start.clone();
                let mut outcomes = vec![false; qubits.len()];
                batched.measure_many_with(&qubits, coins(7), &mut outcomes, transpose);
                assert_eq!(outcomes, expected, "n = {n}, transpose = {transpose}");
                assert_eq!(batched, single, "n = {n}, transpose = {transpose}");
            }
        }
    }

    #[test]
    fn reset_leaves_zero() {
        for (i, n) in [1, 2, 7, 65].into_iter().enumerate() {
            let mut t = InverseTableau::new(n);
            random_circuit(n, 100, 0x7e5e + i as u64)
                .into_iter()
                .for_each(|g| apply(&mut t, g));
            for q in 0..n {
                t.reset(q, || true);
                assert!(!t.measure(q, || panic!("reset qubit is random")));
            }
        }
    }
}
