// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Circuit tokens parse through the composite's syntax, as they appear in `.sst`.
//!
//! Several tokens share a prefix, so the dialect's declaration order is
//! load-bearing: the longer token must win.

use chumsky::Parser as _;
use ppvm_vihaco::composite::PPVM;
use vihaco::{HasInstructionSet, Parse};

type Syntax = <PPVM as HasInstructionSet>::Syntax;

/// Parse `circuit::circuit::<token>` and return the instruction's debug form.
fn parse(token: &str) -> Option<String> {
    let src = format!("circuit::circuit::{token}");
    Syntax::parser()
        .parse(src.as_str())
        .into_result()
        .ok()
        .map(|instruction| format!("{instruction:?}"))
}

fn assert_parses_to(token: &str, variant: &str) {
    let expected = format!("V_circuit__circuit__{variant}({variant})");
    assert_eq!(parse(token), Some(expected), "token `{token}`");
}

#[test]
fn parses_t_family_without_prefix_collision() {
    // `t` is a prefix of `t_adj`, `trace`, `truncate`, and `twoqubitpaulierror`.
    assert_parses_to("t", "T");
    assert_parses_to("t_adj", "TAdj");
    assert_parses_to("trace", "Trace");
    assert_parses_to("truncate", "Truncate");
    assert_parses_to("twoqubitpaulierror", "TwoQubitPauliError");
}

#[test]
fn parses_s_family_without_prefix_collision() {
    // `s` is a prefix of `s_adj`, `sqrt_x`, `sqrt_y`, etc.
    assert_parses_to("s", "S");
    assert_parses_to("s_adj", "SAdj");
    assert_parses_to("sqrt_x", "SqrtX");
    assert_parses_to("sqrt_x_adj", "SqrtXAdj");
    assert_parses_to("sqrt_y", "SqrtY");
    assert_parses_to("sqrt_y_adj", "SqrtYAdj");
}

#[test]
fn rejects_unknown_token() {
    assert_eq!(parse("nope"), None);
    for token in ["tadj", "sadj", "sqrtx", "sqrtxadj", "sqrty", "sqrtyadj"] {
        assert_eq!(parse(token), None, "legacy compact token `{token}`");
    }
}

#[test]
fn rejects_pascal_case_token() {
    // Tokens are lowercase; the Rust type name must not parse.
    assert_eq!(parse("CNOT"), None);
}

#[test]
fn parses_batched_tokens_without_prefix_collision() {
    // Each base token is a prefix of its batched form, e.g. `x` of `x_many`.
    for (token, variant) in [
        ("x_many", "XMany"),
        ("y_many", "YMany"),
        ("z_many", "ZMany"),
        ("h_many", "HMany"),
        ("s_many", "SMany"),
        ("s_adj_many", "SAdjMany"),
        ("t_many", "TMany"),
        ("t_adj_many", "TAdjMany"),
        ("sqrt_x_many", "SqrtXMany"),
        ("sqrt_x_adj_many", "SqrtXAdjMany"),
        ("sqrt_y_many", "SqrtYMany"),
        ("sqrt_y_adj_many", "SqrtYAdjMany"),
        ("reset_many", "ResetMany"),
        ("measure_many", "MeasureMany"),
        ("cnot_many", "CNOTMany"),
        ("cz_many", "CZMany"),
        ("rx_many", "RXMany"),
        ("ry_many", "RYMany"),
        ("rz_many", "RZMany"),
        ("depolarize_many", "DepolarizeMany"),
        ("loss_many", "LossMany"),
        ("rxx_many", "RXXMany"),
        ("ryy_many", "RYYMany"),
        ("rzz_many", "RZZMany"),
        ("depolarize2_many", "Depolarize2Many"),
        ("u3_many", "U3Many"),
        ("paulierror_many", "PauliErrorMany"),
        ("correlatedloss_many", "CorrelatedLossMany"),
        ("twoqubitpaulierror_many", "TwoQubitPauliErrorMany"),
    ] {
        assert_parses_to(token, variant);
    }
}
