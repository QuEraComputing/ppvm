// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Z measurement and reset of the stabilizer state `C |0…0⟩` (no active
//! qubits).

use ppvm_pauli_kernels::heisenberg;
use ppvm_traits_2::Clifford;

use super::InverseTableau;
use super::bitword::set_bits;

impl InverseTableau {
    /// Measure `Z_q` on `C |0…0⟩` and collapse; returns the outcome (`true`
    /// for `|1⟩`). `coin` is called only if the outcome is random.
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
    pub fn measure(&mut self, q: usize, coin: impl FnOnce() -> bool) -> bool {
        let mut x_support = set_bits(self.iz(q).x);
        let Some(pivot) = x_support.next() else {
            return self.iz(q).negative;
        };
        let targets: Vec<usize> = x_support.collect();
        // The CXs reduce O to s · X_p or s · Y_p (times Z's); find which, and s,
        // from row iz_q alone before changing the tableau.
        let row = self.iz(q);
        let (mut pivot_bits, mut negative) = ([true, row.z_bit(pivot)], row.negative);
        for &t in &targets {
            let ([px, pz, _, _], flip) =
                heisenberg::cnot([pivot_bits[0], pivot_bits[1], row.x_bit(t), row.z_bit(t)]);
            pivot_bits = [px, pz];
            negative ^= flip;
        }
        let (has_y, outcome) = (pivot_bits[1], coin());
        let fix_sign = outcome != negative;
        self.fan_out(pivot, &targets, &[], |x, z| {
            let ([x, z], f1) = if has_y {
                heisenberg::s(x, z)
            } else {
                ([x, z], false)
            };
            let ([x, z], f2) = heisenberg::h(x, z);
            let ([x, z], f3) = if fix_sign {
                heisenberg::x(x, z)
            } else {
                ([x, z], false)
            };
            ([x, z], f1 ^ f2 ^ f3)
        });
        outcome
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
