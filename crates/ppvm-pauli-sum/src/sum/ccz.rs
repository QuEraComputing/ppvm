// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use crate::sum::PauliSum;
use ppvm_traits::{config::Config, traits::*};

impl<T: Config> CCZGate<T> for PauliSum<T> {
    fn ccz(&mut self, a: usize, b: usize, c: usize) {
        let targets = [a, b, c];
        assert!(a != b && a != c && b != c, "CCZ requires distinct qubits");
        assert!(targets.iter().all(|&q| q < self.n_qubits()));
        self.map_insert_multiple(|word, coeff| {
            if targets.iter().any(|&q| word.get_lbit(q)) {
                return None;
            }
            let x_mask = targets
                .iter()
                .enumerate()
                .fold(0u8, |mask, (i, &q)| mask | ((word.get_xbit(q) as u8) << i));
            if x_mask == 0 {
                return None;
            }

            // CCZ changes sign only at |111>. Conjugating a word with X mask
            // x therefore gives P/2 - Σ_m (-1)^|m| P Z^m/2, where m is nonzero
            // and has even overlap with x. These are exactly three new terms.
            *coeff = coeff.half();
            let mut branches = Vec::with_capacity(3);
            for mask in 1u8..8 {
                if (mask & x_mask).count_ones() % 2 != 0 {
                    continue;
                }
                let mut branch = word.clone();
                let mut phase = 0i8;
                for (i, &q) in targets.iter().enumerate() {
                    if mask & (1 << i) != 0 {
                        // Multiplication on the right: XZ = -iY, YZ = iX.
                        if word.get_xbit(q) {
                            phase += if word.get_zbit(q) { 1 } else { -1 };
                        }
                        branch.set_zbit(q, !word.get_zbit(q));
                    }
                }
                branch.rehash();
                let mut sign = if mask.count_ones() % 2 == 0 { -1 } else { 1 };
                // Even overlap makes the product real: phase is 0 or ±2.
                if phase != 0 {
                    sign = -sign;
                }
                branches.push((branch, coeff.mul_sign(sign)));
            }
            Some(branches)
        });
    }
}
