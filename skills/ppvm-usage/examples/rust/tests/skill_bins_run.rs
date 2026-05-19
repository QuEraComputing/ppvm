//! Integration test: every Rust skill example actually runs.
//!
//! `cargo build` catches signature changes; this catches runtime
//! regressions in the public API surface that the skill teaches.

use std::process::Command;

fn run_bin(name: &str) {
    let status = Command::new(env!("CARGO"))
        .args(["run", "--quiet", "-p", "ppvm-skill-examples", "--bin", name])
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn cargo for `{name}`: {e}"));
    assert!(
        status.success(),
        "skill example `{name}` exited with {status:?}"
    );
}

#[test]
fn paulisum_skill_example_runs() {
    run_bin("paulisum");
}

#[test]
fn tableau_skill_example_runs() {
    run_bin("tableau");
}

#[test]
fn stim_sample_skill_example_runs() {
    run_bin("stim_sample");
}

#[test]
fn sym_skill_example_runs() {
    run_bin("sym");
}
