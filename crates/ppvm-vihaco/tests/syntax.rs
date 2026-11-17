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
