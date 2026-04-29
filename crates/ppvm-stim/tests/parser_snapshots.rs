//! Snapshot tests for the Stim parser.
//!
//! Each test parses a small representative program and snapshots the resulting
//! `Program` debug representation. Re-running flags any AST drift cheaply.
//!
//! On first run insta writes pending snapshots into `tests/snapshots/`.
//! Review with `cargo insta review` (install with `cargo install cargo-insta`).

use ppvm_stim::parse;

fn snap(name: &str, src: &str) {
    let prog = parse(src).unwrap_or_else(|e| panic!("parse failed for {name}: {e}"));
    insta::assert_debug_snapshot!(name, prog);
}

#[test]
fn snapshot_bare_gate() {
    snap("bare_gate", "H 0\n");
}

#[test]
fn snapshot_tagged_gate() {
    snap("tagged_gate", "S[T] 0\n");
}

#[test]
fn snapshot_tagged_gate_with_named_params() {
    snap(
        "tagged_gate_with_named_params",
        "I[R_X(theta=0.5*pi)] 0\n",
    );
}

#[test]
fn snapshot_multi_tag() {
    snap("multi_tag", "S[T,debug] 0\n");
}

#[test]
fn snapshot_args_only_noise() {
    snap("args_only_noise", "DEPOLARIZE1(0.5) 0\n");
}

#[test]
fn snapshot_tag_plus_args() {
    snap(
        "tag_plus_args",
        "I_ERROR[correlated_loss](0.1, 0.2, 0.3) 0 1\n",
    );
}

#[test]
fn snapshot_single_target_measurement() {
    snap("single_target_measurement", "M 0\n");
}

#[test]
fn snapshot_multi_target_measurement_with_noise() {
    snap("multi_target_measurement_with_noise", "M(0.001) 0 1 2\n");
}

#[test]
fn snapshot_annotation_with_rec_target() {
    snap("annotation_with_rec_target", "DETECTOR rec[-1]\n");
}

#[test]
fn snapshot_annotation_with_args() {
    snap("annotation_with_args", "OBSERVABLE_INCLUDE(0)\n");
}

#[test]
fn snapshot_empty_repeat_body() {
    snap("empty_repeat_body", "REPEAT 5 { }\n");
}

#[test]
fn snapshot_multi_instruction_repeat() {
    snap(
        "multi_instruction_repeat",
        "REPEAT 3 {\n    H 0\n    CX 0 1\n    M 0 1\n}\n",
    );
}

#[test]
fn snapshot_nested_repeat() {
    snap(
        "nested_repeat",
        "REPEAT 2 {\n    REPEAT 3 {\n        X 0\n    }\n}\n",
    );
}

#[test]
fn snapshot_comment_heavy_program() {
    snap(
        "comment_heavy_program",
        "# Top comment.\n\n# Pre-gate comment.\nH 0\n# Mid comment.\nM 0\n# Trailing-line comment.\n",
    );
}

#[test]
fn snapshot_whitespace_stress_program() {
    snap(
        "whitespace_stress_program",
        "  \t  H   0  \t \nCX  0   1\n   M  0  1  \n",
    );
}
