//! Skill code block: ppvm-stim multi-shot sampling.
//!
//! Mirrors the §6 "Running Stim programs (Rust)" example in
//! `skills/ppvm-usage/SKILL.md`. The pattern matters: parse once,
//! reuse the parsed program across shots via a factory closure.

use ppvm_stim::{parse_extended, sample};
use ppvm_tableau::prelude::*;

type Tab = GeneralizedTableau<ppvm_runtime::config::indexmap::ByteFxHashF64<1>, usize>;

fn main() -> Result<(), ppvm_stim::Error> {
    let stim_src = r#"
        H 0
        CX 0 1
        M 0 1
    "#;
    let n_qubits = 2;

    let prog = parse_extended(stim_src)?;
    // The factory receives the shot index; derive a per-shot seed from it
    // when you need deterministic, order-independent results.
    let shots = sample(&prog, 32, |_| Tab::new(n_qubits, 1e-10))?;

    // GHZ: every shot has the two qubits in agreement.
    for (i, shot) in shots.iter().enumerate() {
        let pair: Vec<_> = shot.iter().collect();
        assert_eq!(pair.len(), 2, "shot {i} has {} outcomes", pair.len());
        match (pair[0], pair[1]) {
            (Some(a), Some(b)) => assert_eq!(a, b, "shot {i} uncorrelated: {a} vs {b}"),
            other => panic!("shot {i} had a loss outcome: {other:?}"),
        }
    }
    println!("shots={} correlated=true", shots.len());
    Ok(())
}
