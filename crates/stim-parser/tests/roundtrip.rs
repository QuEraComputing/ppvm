// SPDX-FileCopyrightText: 2026 QuEra Computing Inc.
// SPDX-License-Identifier: Apache-2.0

//! Printer-fixpoint roundtrip tests for both AST layers.
//!
//! The grammar accepts a much wider input surface than the printer
//! produces (comments, blank lines, varying whitespace, `rec[-1]`
//! annotation targets, `pi`-expressions, …). The printer normalizes all
//! of that away. After one normalization round the output is a fixpoint:
//! parsing the printed form and reprinting must produce byte-identical
//! output.
//!
//! Property checked, for both `parse` and `parse_extended`:
//!
//! ```text
//! src --parse--> ast1 --print--> s1 --parse--> _ast2 --print--> s2
//! assert s1 == s2
//! ```
//!
//! Note: we do not compare ASTs directly because `line` numbers are
//! source-derived and shift when leading comments / blank lines get
//! normalized away.

use stim_parser::extended::parse_extended;
use stim_parser::prelude::parse;

const VANILLA_CORPUS: &[(&str, &str)] = &[
    ("bell_pair", "H 0\nCX 0 1\nM 0 1\n"),
    (
        "ghz4",
        "# 4-qubit GHZ\nH 0\nCX 0 1\nCX 1 2\nCX 2 3\nM 0 1 2 3\n",
    ),
    ("repeat_block", "REPEAT 5 {\n    X 0\n}\nM 0\n"),
    ("depolarize_smoke", "H 0\nDEPOLARIZE1(0.05) 0\nM 0\n"),
    (
        "mx_unsupported",
        // parses fine — only execution rejects it
        "MX 0 1\n",
    ),
    ("swap_unsupported", "SWAP 0 1\nM 0 1\n"),
    (
        "annotation_drops_rec",
        // Non-numeric annotation targets get dropped during parse; the
        // fixpoint property still holds because both rounds drop them.
        "R 0 1\nMR 0\nDETECTOR rec[-1]\nOBSERVABLE_INCLUDE(0) rec[-1]\nTICK\n",
    ),
    (
        "rep_code_d3_r3",
        "R 0 1 2 3 4\nREPEAT 3 {\n    CX 0 3\n    CX 1 3\n    CX 1 4\n    CX 2 4\n    MR 3 4\n    DETECTOR rec[-1]\n    DETECTOR rec[-2]\n    TICK\n}\nM 0 1 2\nOBSERVABLE_INCLUDE(0) rec[-1]\n",
    ),
    (
        "nested_repeat",
        "REPEAT 2 {\n    REPEAT 3 {\n        H 0\n        M 0\n    }\n}\n",
    ),
    (
        "tags_and_args",
        "DEPOLARIZE1(0.1) 0 1 2\nPAULI_CHANNEL_1(0.01, 0.02, 0.03) 0\nM(0.001) 0 1\nMPAD(0.05) 0 1 0 1\n",
    ),
    (
        "noisy_measurements",
        "MZ(0.25) 0\nMR(0.001) 1\nMPAD 0 1 0\n",
    ),
    (
        "comments_and_blank_lines",
        "# leading comment\n\nH 0  # trailing\n\n# mid\nM 0\n",
    ),
    (
        "pi_expression_args",
        // `pi` and `0.5*pi` parse to plain f64, then print as decimal.
        // The fixpoint check ignores the original spelling.
        "PAULI_CHANNEL_1(0.5, 0.5, 0) 0\n",
    ),
];

const EXTENDED_CORPUS: &[(&str, &str)] = &[
    ("vanilla_h", "H 0\n"),
    ("t_sugar", "S[T] 0 1\nS_DAG[T] 2\n"),
    (
        "rotations",
        "I[R_X(theta=0.5)] 0\nI[R_Y(theta=1.25)] 1\nI[R_Z(theta=-0.5)] 2\n",
    ),
    ("u3", "I[U3(theta=0.5, phi=1.0, lambda=1.5)] 0\n"),
    ("loss", "I_ERROR[loss](0.01) 0 1\n"),
    (
        "correlated_loss",
        "I_ERROR[correlated_loss](0.1, 0.05, 0.05) 0 1 2 3\n",
    ),
    ("mpad_bits", "MPAD 0 1 0\nMPAD(0.01) 1 1 0 0\n"),
    (
        "extended_in_repeat",
        "REPEAT 3 {\n    S[T] 0\n    I[R_Z(theta=0.25)] 0\n    M 0\n}\n",
    ),
    (
        "mixed_extended_and_vanilla",
        "H 0\nS[T] 1\nCX 0 1\nI_ERROR[loss](0.01) 0\nM 0 1\n",
    ),
];

fn assert_raw_fixpoint(name: &str, src: &str) {
    let ast1 = parse(src).unwrap_or_else(|e| panic!("[{name}] parse failed: {e}"));
    let s1 = format!("{ast1}");
    let ast2 = parse(&s1).unwrap_or_else(|e| panic!("[{name}] reparse failed: {e}\n--s1--\n{s1}"));
    let s2 = format!("{ast2}");
    assert_eq!(s1, s2, "[{name}] print/reparse not a fixpoint");
    let _ = ast2;
}

fn assert_extended_fixpoint(name: &str, src: &str) {
    let ast1 = parse_extended(src).unwrap_or_else(|e| panic!("[{name}] parse failed: {e}"));
    let s1 = format!("{ast1}");
    let ast2 = parse_extended(&s1)
        .unwrap_or_else(|e| panic!("[{name}] reparse failed: {e}\n--s1--\n{s1}"));
    let s2 = format!("{ast2}");
    assert_eq!(s1, s2, "[{name}] print/reparse not a fixpoint");
    let _ = ast2;
}

#[test]
fn vanilla_corpus_fixpoint() {
    for (name, src) in VANILLA_CORPUS {
        assert_raw_fixpoint(name, src);
    }
}

#[test]
fn extended_corpus_fixpoint() {
    for (name, src) in EXTENDED_CORPUS {
        assert_extended_fixpoint(name, src);
    }
}

#[test]
fn empty_program_fixpoint() {
    assert_raw_fixpoint("empty", "");
    assert_extended_fixpoint("empty", "");
}

#[test]
fn whitespace_only_fixpoint() {
    assert_raw_fixpoint("whitespace", "  \n# only comments\n\n");
}

#[test]
fn vanilla_printed_form_is_canonical_shape() {
    // Spot-check the printer output for one representative input so
    // that the chosen canonical shape is documented.
    let src = "H 0  # trail\nCX  0   1\nDEPOLARIZE1(0.05) 0 1\nREPEAT 2 { X 0 }\n";
    let ast = parse(src).unwrap();
    let printed = format!("{ast}");
    let expected = "\
H 0
CX 0 1
DEPOLARIZE1(0.05) 0 1
REPEAT 2 {
    X 0
}
";
    assert_eq!(printed, expected);
}

#[test]
fn extended_printed_form_lowers_sugar_into_canonical_stim() {
    let src = "S[T] 0\nI[R_X(theta=0.25)] 1\nI_ERROR[loss](0.01) 2\n";
    let ast = parse_extended(src).unwrap();
    let printed = format!("{ast}");
    let expected = "\
S[T] 0
I[R_X(theta=0.25)] 1
I_ERROR[loss](0.01) 2
";
    assert_eq!(printed, expected);
}
