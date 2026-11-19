// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! The rules against the matrix-model tables, and integer widths against `bool`.

mod common;
#[allow(dead_code)]
#[path = "generated/tables.rs"]
#[rustfmt::skip]
mod tables;

use common::Rng;
use ppvm_pauli_kernels::{Bits, heisenberg};

type Rule1<B> = fn(B, B) -> ([B; 2], B);
type Rule2<B> = fn([B; 4]) -> ([B; 4], B);
/// A gate's name, its rule at two widths, and its matrix-model table.
type Case1 = (
    &'static str,
    Rule1<bool>,
    Rule1<u64>,
    &'static [[bool; 5]; 4],
);
type Case2 = (
    &'static str,
    Rule2<bool>,
    Rule2<u64>,
    &'static [[bool; 9]; 16],
);

fn one_qubit() -> [Case1; 10] {
    [
        ("X", heisenberg::x, heisenberg::x, &tables::X),
        ("Y", heisenberg::y, heisenberg::y, &tables::Y),
        ("Z", heisenberg::z, heisenberg::z, &tables::Z),
        ("H", heisenberg::h, heisenberg::h, &tables::H),
        ("S", heisenberg::s, heisenberg::s, &tables::S),
        (
            "S_DAG",
            heisenberg::s_dag,
            heisenberg::s_dag,
            &tables::S_DAG,
        ),
        (
            "SQRT_X",
            heisenberg::sqrt_x,
            heisenberg::sqrt_x,
            &tables::SQRT_X,
        ),
        (
            "SQRT_X_DAG",
            heisenberg::sqrt_x_dag,
            heisenberg::sqrt_x_dag,
            &tables::SQRT_X_DAG,
        ),
        (
            "SQRT_Y",
            heisenberg::sqrt_y,
            heisenberg::sqrt_y,
            &tables::SQRT_Y,
        ),
        (
            "SQRT_Y_DAG",
            heisenberg::sqrt_y_dag,
            heisenberg::sqrt_y_dag,
            &tables::SQRT_Y_DAG,
        ),
    ]
}

fn two_qubit() -> [Case2; 3] {
    [
        ("CX", heisenberg::cnot, heisenberg::cnot, &tables::CX),
        ("CZ", heisenberg::cz, heisenberg::cz, &tables::CZ),
        ("CY", heisenberg::cy, heisenberg::cy, &tables::CY),
    ]
}

#[test]
fn one_qubit_rules_match_matrix_model() {
    for (name, rule, _, table) in one_qubit() {
        for row in table {
            let (bits, flip) = rule(row[0], row[1]);
            assert_eq!(
                (bits, flip),
                ([row[2], row[3]], row[4]),
                "{name} on {:?}",
                &row[..2]
            );
        }
    }
}

#[test]
fn two_qubit_rules_match_matrix_model() {
    for (name, rule, _, table) in two_qubit() {
        for row in table {
            let (bits, flip) = rule([row[0], row[1], row[2], row[3]]);
            assert_eq!(
                (bits, flip),
                ([row[4], row[5], row[6], row[7]], row[8]),
                "{name} on {:?}",
                &row[..4]
            );
        }
    }
}

fn bit(w: u64, i: usize) -> bool {
    (w >> i) & 1 == 1
}

#[test]
fn u64_rules_equal_bool_rules_per_bit() {
    let mut rng = Rng::new(7);
    for _ in 0..64 {
        let [x, z] = [rng.next_u64(), rng.next_u64()];
        for (name, rule_b, rule_w, _) in one_qubit() {
            let ([xw, zw], fw) = rule_w(x, z);
            for i in 0..64 {
                let expect = rule_b(bit(x, i), bit(z, i));
                assert_eq!(
                    ([bit(xw, i), bit(zw, i)], bit(fw, i)),
                    expect,
                    "{name} bit {i}"
                );
            }
        }
        let input = [
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
        ];
        for (name, rule_b, rule_w, _) in two_qubit() {
            let (out, fw) = rule_w(input);
            for i in 0..64 {
                let expect = rule_b(input.map(|w| bit(w, i)));
                assert_eq!(
                    (out.map(|w| bit(w, i)), bit(fw, i)),
                    expect,
                    "{name} bit {i}"
                );
            }
        }
    }
}

/// The rules only use `Bits`, so narrower integers behave like `u64`.
fn narrow_matches<B: Bits + Into<u64> + TryFrom<u64>>(rng: &mut Rng)
where
    <B as TryFrom<u64>>::Error: std::fmt::Debug,
{
    let mask = |w: u64| w & (u64::MAX >> (64 - 8 * std::mem::size_of::<B>()));
    let to_b = |w: u64| B::try_from(mask(w)).unwrap();
    let (x, z) = (rng.next_u64(), rng.next_u64());
    let (out, f) = heisenberg::cy([to_b(x), to_b(z), to_b(!x), to_b(z ^ x)]);
    let (out_w, f_w) = heisenberg::cy([x, z, !x, z ^ x]);
    assert_eq!(out.map(Into::into), out_w.map(mask));
    assert_eq!(f.into(), mask(f_w));
}

#[test]
fn narrow_integer_rules_equal_u64_rules() {
    let mut rng = Rng::new(11);
    for _ in 0..64 {
        narrow_matches::<u8>(&mut rng);
        narrow_matches::<u16>(&mut rng);
        narrow_matches::<u32>(&mut rng);
    }
}
