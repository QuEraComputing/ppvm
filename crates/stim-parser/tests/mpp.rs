// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Parsing of multi-qubit Pauli-product measurement targets (`MPP`).

use stim_parser::prelude::{Instruction, MppOp, PauliAxis, PauliFactor, parse};

fn products(src: &str) -> Vec<Vec<PauliFactor>> {
    let prog = parse(src).expect("must parse");
    match &prog.instructions[0] {
        Instruction::Mpp(MppOp { products, .. }) => products.clone(),
        other => panic!("expected Mpp, got {other:?}"),
    }
}

fn f(axis: PauliAxis, qubit: usize) -> PauliFactor {
    PauliFactor { axis, qubit }
}

#[test]
fn single_product_parses_all_factors() {
    let p = products("MPP X0*Y1*Z2");
    assert_eq!(
        p,
        vec![vec![
            f(PauliAxis::X, 0),
            f(PauliAxis::Y, 1),
            f(PauliAxis::Z, 2)
        ]]
    );
}

#[test]
fn space_separated_products_are_distinct() {
    let p = products("MPP Z0 Z1*Z2");
    assert_eq!(
        p,
        vec![
            vec![f(PauliAxis::Z, 0)],
            vec![f(PauliAxis::Z, 1), f(PauliAxis::Z, 2)],
        ]
    );
}

#[test]
fn high_qubit_indices_parse() {
    let p = products("MPP Y0*Y40");
    assert_eq!(p, vec![vec![f(PauliAxis::Y, 0), f(PauliAxis::Y, 40)]]);
}

#[test]
fn invalid_pauli_letter_is_rejected() {
    assert!(parse("MPP Q0").is_err());
    assert!(parse("MPP X0*W1").is_err());
}

#[test]
fn missing_qubit_index_is_rejected() {
    assert!(parse("MPP X").is_err());
    assert!(parse("MPP X0*Y").is_err());
}

#[test]
fn mpp_with_no_targets_is_rejected() {
    assert!(parse("MPP").is_err());
}

#[test]
fn mpp_round_trips_through_display() {
    let src = "MPP X0*Y3*Z7 Z1*Z2";
    let printed = parse(src).unwrap().to_string();
    assert_eq!(printed.trim(), "MPP X0*Y3*Z7 Z1*Z2");
    assert_eq!(parse(printed.trim()).unwrap(), parse(src).unwrap());
}
