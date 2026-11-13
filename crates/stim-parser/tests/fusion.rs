// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Line fusion (#291). Expected outputs are `str(stim.Circuit(src))` from
//! Stim 1.15, except where noted.

use stim_parser::diagnostics::{Collect, FailFast};
use stim_parser::prelude::*;

/// `(source, what Stim prints for it)`.
const STIM_PARITY: &[(&str, &str)] = &[
    // --- fused: the cases from the issue ---
    ("H 0\nH 1", "H 0 1\n"),
    ("CX 0 1\nCX 2 3", "CX 0 1 2 3\n"),
    ("M 0\nM 1", "M 0 1\n"),
    ("X_ERROR(0.1) 0\nX_ERROR(0.1) 1", "X_ERROR(0.1) 0 1\n"),
    ("H[tag] 0\nH[tag] 1", "H[tag] 0 1\n"),
    ("H 0\nH 0", "H 0 0\n"), // targets are not deduplicated
    // --- fused: other families ---
    ("R 0\nR 1", "R 0 1\n"),
    ("MR(0.1) 0\nMR(0.1) 1", "MR(0.1) 0 1\n"),
    ("MPP X0\nMPP Z1", "MPP X0 Z1\n"),
    ("MPAD 0\nMPAD 1", "MPAD 0 1\n"),
    (
        "DEPOLARIZE2(0.1) 0 1\nDEPOLARIZE2(0.1) 2 3",
        "DEPOLARIZE2(0.1) 0 1 2 3\n",
    ),
    (
        "HERALDED_ERASE(0.1) 0\nHERALDED_ERASE(0.1) 1",
        "HERALDED_ERASE(0.1) 0 1\n",
    ),
    (
        "I_ERROR[loss](0.1) 0\nI_ERROR[loss](0.1) 1",
        "I_ERROR[loss](0.1) 0 1\n",
    ),
    ("CX rec[-1] 0\nCX rec[-1] 1", "CX rec[-1] 0 rec[-1] 1\n"),
    // args compare as numbers, not as source text
    ("X_ERROR(0.1) 0\nX_ERROR(0.10) 1", "X_ERROR(0.1) 0 1\n"),
    // comments and blank lines are invisible
    ("H 0\n# c\n\nH 1", "H 0 1\n"),
    // fusion also happens inside REPEAT bodies
    ("REPEAT 2 {\nH 0\nH 1\n}", "REPEAT 2 {\n    H 0 1\n}\n"),
    // runs of more than two lines
    ("H 0\nH 1\nH 2\nX 3\nX 4", "H 0 1 2\nX 3 4\n"),
    // --- not fused ---
    (
        "X_ERROR(0.1) 0\nX_ERROR(0.2) 1",
        "X_ERROR(0.1) 0\nX_ERROR(0.2) 1\n",
    ),
    ("M(0.1) 0\nM 1", "M(0.1) 0\nM 1\n"),
    ("H[a] 0\nH[b] 1", "H[a] 0\nH[b] 1\n"),
    ("H 0\nTICK\nH 1", "H 0\nTICK\nH 1\n"),
    ("H 0\nX 1\nH 2", "H 0\nX 1\nH 2\n"),
    ("TICK\nTICK", "TICK\nTICK\n"),
    (
        "DETECTOR rec[-1]\nDETECTOR rec[-1]",
        "DETECTOR\nDETECTOR\n", // rec targets on annotations are dropped by this parser
    ),
    (
        "QUBIT_COORDS(0.5) 0\nQUBIT_COORDS(0.5) 1",
        "QUBIT_COORDS(0.5) 0\nQUBIT_COORDS(0.5) 1\n",
    ),
    (
        "SHIFT_COORDS(0.5)\nSHIFT_COORDS(0.5)",
        "SHIFT_COORDS(0.5)\nSHIFT_COORDS(0.5)\n",
    ),
    (
        "CORRELATED_ERROR(0.1) 0\nCORRELATED_ERROR(0.1) 1",
        "CORRELATED_ERROR(0.1) 0\nCORRELATED_ERROR(0.1) 1\n",
    ),
    (
        "REPEAT 2 {\nH 0\n}\nREPEAT 2 {\nH 0\n}",
        "REPEAT 2 {\n    H 0\n}\nREPEAT 2 {\n    H 0\n}\n",
    ),
    ("REPEAT 2 {\nH 0\n}\nH 1", "REPEAT 2 {\n    H 0\n}\nH 1\n"),
    ("H 0\nREPEAT 2 {\nH 1\n}", "H 0\nREPEAT 2 {\n    H 1\n}\n"),
    // spelling aliases are deliberately not fused (Stim would fuse these)
    ("CX 0 1\nCNOT 2 3", "CX 0 1\nCNOT 2 3\n"),
    ("M 0\nMZ 1", "M 0\nMZ 1\n"),
];

/// Parse each source with `parse_fn` and compare the printed program.
fn assert_prints<P: std::fmt::Display>(
    cases: &[(&str, &str)],
    parse_fn: impl Fn(&str) -> Result<P, Diagnostics>,
) {
    for (src, expected) in cases {
        let printed = parse_fn(src)
            .unwrap_or_else(|e| panic!("{src:?}: {e}"))
            .to_string();
        assert_eq!(&printed, expected, "source: {src:?}");
    }
}

#[test]
fn fusion_matches_stim() {
    assert_prints(STIM_PARITY, parse);
}

#[test]
fn per_line_arity_is_checked_before_fusing() {
    // Each line alone has an odd target count; Stim rejects this too.
    let err = parse("CX 0\nCX 1").unwrap_err();
    assert_eq!(err.iter().next().unwrap().code, Some("target-count"));
}

#[test]
fn fused_span_covers_the_whole_run() {
    // Spans cover the name token: first line's name to last line's name.
    let src = "H 0\n# c\nH 1\nX 2";
    let p = parse(src).unwrap();
    let Instruction::Gate(op) = &p.instructions[0] else {
        panic!("{:?}", p.instructions[0]);
    };
    assert_eq!(&src[op.span.start..op.span.end], "H 0\n# c\nH");
    assert_eq!(op.span.line(&p.line_map), 1);
}

// ---- extended AST: fusion runs after lowering ----------------------------

/// `(source, what parse_extended prints for it)`.
const EXTENDED: &[(&str, &str)] = &[
    ("S[T] 0\nS[T] 1", "S[T] 0 1\n"),
    // The extended AST describes operations, not spelling: both lines are `T`.
    ("S[T] 0\nT 1", "S[T] 0 1\n"),
    ("T[a] 0\nT[b] 1", "S[T] 0 1\n"),
    ("S_DAG[T] 0\nT 1", "S_DAG[T] 0\nS[T] 1\n"),
    (
        "I[R_X(theta=0.5*pi)] 0\nI[R_X(theta=0.5*pi)] 1",
        "I[R_X(theta=0.5*pi)] 0 1\n",
    ),
    (
        "I[R_X(theta=0.5*pi)] 0\nI[R_X(theta=0.25*pi)] 1",
        "I[R_X(theta=0.5*pi)] 0\nI[R_X(theta=0.25*pi)] 1\n",
    ),
    (
        "I[R_X(theta=0.5*pi)] 0\nI[R_Z(theta=0.5*pi)] 1",
        "I[R_X(theta=0.5*pi)] 0\nI[R_Z(theta=0.5*pi)] 1\n",
    ),
    (
        "I[U3(theta=0.5*pi, phi=1.0*pi, lambda=1.5*pi)] 0\nI[U3(theta=0.5*pi, phi=1.0*pi, lambda=1.5*pi)] 1",
        "I[U3(theta=0.5*pi, phi=1.0*pi, lambda=1.5*pi)] 0 1\n",
    ),
    (
        "I_ERROR[loss](0.1) 0\nI_ERROR[loss](0.1) 1",
        "I_ERROR[loss](0.1) 0 1\n",
    ),
    (
        "I_ERROR[loss](0.1) 0\nI_ERROR[loss](0.2) 1",
        "I_ERROR[loss](0.1) 0\nI_ERROR[loss](0.2) 1\n",
    ),
    (
        "I_ERROR[leakage](0.1, 0.2) 0\nI_ERROR[leakage](0.1, 0.2) 1",
        "I_ERROR[leakage](0.1, 0.2) 0 1\n",
    ),
    (
        "I_ERROR[correlated_loss](0.1) 0 1\nI_ERROR[correlated_loss](0.1) 2 3",
        "I_ERROR[correlated_loss](0.1, 0.0, 0.0) 0 1 2 3\n",
    ),
    (
        "I_ERROR[reset_loss] 0\nI_ERROR[reset_loss] 1\nI_ERROR[reset_leakage] 2",
        "I_ERROR[reset_loss] 0 1\nI_ERROR[reset_leakage] 2\n",
    ),
    ("H 0\nH 1\nM 0\nM 1", "H 0 1\nM 0 1\n"),
    (
        "REPEAT 2 {\nS[T] 0\nT 1\n}",
        "REPEAT 2 {\n    S[T] 0 1\n}\n",
    ),
];

#[test]
fn extended_fusion() {
    assert_prints(EXTENDED, parse_extended);
}

#[test]
fn vanilla_keeps_spellings_apart() {
    // Plain `parse` follows Stim: different names never fuse.
    assert_eq!(parse("S[T] 0\nT 1").unwrap().instructions.len(), 2);
}

#[test]
fn lowering_checks_each_line_before_fusing() {
    // Fusing first would have turned two odd target counts into a valid `0 1`.
    let err = parse_extended("I_ERROR[correlated_loss](0.1) 0\nI_ERROR[correlated_loss](0.1) 1")
        .unwrap_err();
    assert_eq!(err.iter().next().unwrap().code, Some("invalid-tag"));
    assert!(err.to_string().starts_with("error at line 1,"), "{err}");
}

#[test]
fn lowering_errors_point_at_their_own_line() {
    // Not reported as part of a fused `T 0 rec[-1]` starting on line 1.
    let err = parse_extended("T 0\nT rec[-1]").unwrap_err();
    assert_eq!(
        err.iter().next().unwrap().code,
        Some("record-target-not-allowed")
    );
    assert!(err.to_string().starts_with("error at line 2,"), "{err}");
}

#[test]
fn unfused_pipeline_keeps_one_instruction_per_line() {
    let mut sink = FailFast::new();
    let p = Pipeline::new("H 0\nH 1")
        .parse(&mut sink)
        .and_then(|p| p.validate(&mut sink))
        .map(|p| p.finish())
        .unwrap();
    assert_eq!(p.instructions.len(), 2);

    let mut sink = FailFast::new();
    let p = Pipeline::new("H 0\nH 1")
        .parse(&mut sink)
        .and_then(|p| p.validate(&mut sink))
        .and_then(|p| p.lower(&mut sink))
        .map(|p| p.finish())
        .unwrap();
    assert_eq!(p.instructions.len(), 2);
}

// ---- error recovery: never fuse across a skipped line ---------------------

#[test]
fn does_not_fuse_across_a_line_skipped_by_validate() {
    // `TIKC` may have been meant as a barrier.
    let mut sink = Collect::new();
    let p = Pipeline::new("H 0\nTIKC\nH 1")
        .parse(&mut sink)
        .and_then(|p| p.validate(&mut sink))
        .map(|p| p.fuse().finish())
        .expect("Collect never aborts");
    assert_eq!(sink.into_items().len(), 1);
    assert_eq!(p.to_string(), "H 0\nH 1\n");
}

#[test]
fn does_not_fuse_across_a_line_skipped_by_lower() {
    let mut sink = Collect::new();
    let p = Pipeline::new("S[T] 0\nI_ERROR[bogus] 0\nS[T] 1\nREPEAT 2 {\nT 0\nI_ERROR 0\nT 1\n}")
        .parse(&mut sink)
        .and_then(|p| p.validate(&mut sink))
        .and_then(|p| p.lower(&mut sink))
        .map(|p| p.fuse().finish())
        .expect("Collect never aborts");
    assert_eq!(sink.into_items().len(), 2);
    assert_eq!(
        p.to_string(),
        "S[T] 0\nS[T] 1\nREPEAT 2 {\n    S[T] 0\n    S[T] 1\n}\n"
    );
}
