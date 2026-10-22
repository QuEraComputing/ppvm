// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_sum::prelude::*;
use std::f64::consts::FRAC_PI_4;

type Sum = PauliSum<config::fxhash::ByteF64<1>>;

fn ccz_via_t_gates(sum: &mut Sum, a: usize, b: usize, c: usize) {
    // Reverse the seven-T circuit for Heisenberg propagation. Rz(π/4)
    // and T have the same conjugation action (their global phase cancels).
    sum.cnot(b, c);
    sum.rz(c, -FRAC_PI_4);
    sum.cnot(a, c);
    sum.rz(c, FRAC_PI_4);
    sum.cnot(b, c);
    sum.rz(c, -FRAC_PI_4);
    sum.cnot(a, c);
    sum.cnot(a, b);
    sum.rz(b, -FRAC_PI_4);
    sum.cnot(a, b);
    for q in [a, b, c] {
        sum.rz(q, FRAC_PI_4);
    }
}

#[test]
fn ccz_matches_t_decomposition_on_every_pauli() {
    for targets in [(0, 1, 2), (4, 0, 2)] {
        for p in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
            for q in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
                for r in [Pauli::I, Pauli::X, Pauli::Y, Pauli::Z] {
                    let mut word = PauliWord::new(5);
                    for (addr, pauli) in [(targets.0, p), (targets.1, q), (targets.2, r)] {
                        word.set(addr, pauli);
                    }
                    let mut native = Sum::builder().n_qubits(5).build();
                    native += (word, 1.0);
                    let mut decomposed = native.clone();
                    native.ccz(targets.0, targets.1, targets.2);
                    ccz_via_t_gates(&mut decomposed, targets.0, targets.1, targets.2);

                    // Subtraction also checks that neither sum has extra terms.
                    for (word, coeff) in native.iter() {
                        decomposed += (*word, -coeff);
                    }
                    assert!(decomposed.iter().all(|(_, coeff)| coeff.abs() < 1e-12));
                }
            }
        }
    }
}

#[test]
fn ccz_preserves_loss_sectors() {
    type LossySum = PauliSum<
        config::fxhash::Byte<1, f64, NoStrategy, LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>>,
    >;
    let mut sum = LossySum::builder().n_qubits(3).build();
    for word in ["LXY", "XLZ", "YXL", "LLX", "LYL", "ZLL", "LLL"] {
        sum += (word, 1.0);
    }
    let before = sum.clone();
    sum.ccz(0, 1, 2);
    assert_eq!(sum, before);
}
