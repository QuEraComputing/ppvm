// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;
use ppvm_traits::traits::{
    AmplitudeDamping as OldAmplitudeDamping, Depolarizing as OldDepolarizing,
    Depolarizing2 as OldDepolarizing2, PauliError as OldPauliError,
    TwoQubitPauliError as OldTwoQubit,
};
use ppvm_traits_2::{
    AmplitudeDamping as NewAmplitudeDamping, Depolarizing as NewDepolarizing,
    Depolarizing2 as NewDepolarizing2, PauliError as NewPauliError,
    TwoQubitPauliError as NewTwoQubit,
};

use super::{NewSum, OldSum, bench_mut};

type OneOld = fn(&mut OldSum, usize, f64);
type OneNew = fn(&mut NewSum, usize, f64);
type ManyOld = fn(&mut OldSum, &[usize], f64);
type ManyNew = fn(&mut NewSum, &[usize], f64);

const P: f64 = 0.037;
const P1: [f64; 3] = [0.011, 0.023, 0.031];
const P2: [f64; 15] = [
    0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013,
    0.014, 0.015,
];
const TARGETS: &[usize] = &[0, 2, 4, 6];
const PAIRS: &[(usize, usize)] = &[(0, 1), (2, 3), (4, 5), (6, 7)];

pub fn bench(c: &mut Criterion) {
    pauli_one(c);
    pauli_two(c);
    depolarizing(c);
    bench_mut(
        c,
        "noise/amplitude_damping",
        |s| s.amplitude_damping(3, P),
        |s| s.amplitude_damping(3, P),
    );
}

fn pauli_one(c: &mut Criterion) {
    bench_mut(
        c,
        "noise/pauli_error",
        |s| s.pauli_error(3, P1),
        |s| s.pauli_error(3, P1, &mut ppvm_conformance_2::analytic_rng()),
    );
    let gates: [(&str, OneOld, OneNew); 3] = [
        (
            "x_error",
            |s, q, p| s.x_error(q, p),
            |s, q, p| s.x_error(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
        (
            "y_error",
            |s, q, p| s.y_error(q, p),
            |s, q, p| s.y_error(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
        (
            "z_error",
            |s, q, p| s.z_error(q, p),
            |s, q, p| s.z_error(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
    ];
    for (name, old, new) in gates {
        bench_mut(
            c,
            &format!("noise/{name}"),
            move |s| old(s, 3, P),
            move |s| new(s, 3, P),
        );
    }
    bench_mut(
        c,
        "noise_batch/pauli_error",
        |s| s.pauli_error_many(TARGETS, P1),
        |s| s.pauli_error_many(TARGETS, P1, &mut ppvm_conformance_2::analytic_rng()),
    );
    let batches: [(&str, ManyOld, ManyNew); 3] = [
        (
            "x_error",
            |s, q, p| s.x_error_many(q, p),
            |s, q, p| s.x_error_many(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
        (
            "y_error",
            |s, q, p| s.y_error_many(q, p),
            |s, q, p| s.y_error_many(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
        (
            "z_error",
            |s, q, p| s.z_error_many(q, p),
            |s, q, p| s.z_error_many(q, p, &mut ppvm_conformance_2::analytic_rng()),
        ),
    ];
    for (name, old, new) in batches {
        bench_mut(
            c,
            &format!("noise_batch/{name}"),
            move |s| old(s, TARGETS, P),
            move |s| new(s, TARGETS, P),
        );
    }
}

fn pauli_two(c: &mut Criterion) {
    bench_mut(
        c,
        "noise/two_qubit_pauli_error",
        |s| s.two_qubit_pauli_error(2, 5, P2),
        |s| s.two_qubit_pauli_error(2, 5, P2, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "noise_batch/two_qubit_pauli_error",
        |s| s.two_qubit_pauli_error_many(PAIRS, P2),
        |s| s.two_qubit_pauli_error_many(PAIRS, P2, &mut ppvm_conformance_2::analytic_rng()),
    );
}

fn depolarizing(c: &mut Criterion) {
    bench_mut(
        c,
        "noise/depolarize1",
        |s| s.depolarize1(3, P),
        |s| s.depolarize1(3, P, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "noise_batch/depolarize1",
        |s| s.depolarize1_many(TARGETS, P),
        |s| s.depolarize1_many(TARGETS, P, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "noise/depolarize2",
        |s| s.depolarize2(2, 5, P),
        |s| s.depolarize2(2, 5, P, &mut ppvm_conformance_2::analytic_rng()),
    );
    bench_mut(
        c,
        "noise_batch/depolarize2",
        |s| s.depolarize2_many(PAIRS, P),
        |s| s.depolarize2_many(PAIRS, P, &mut ppvm_conformance_2::analytic_rng()),
    );
}
