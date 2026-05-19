//! Skill code block: Pauli propagation with truncation.
//!
//! Mirrors the §4 "ppvm-runtime — Pauli propagation" / "Pauli propagation"
//! example in `skills/ppvm-usage/SKILL.md`. Compiled by
//! `cargo build --examples` in CI, so a method-rename or strategy-API
//! change here will break the build instead of silently leaving the
//! skill stale.

use ppvm_runtime::{prelude::*, strategy::CoefficientThreshold};

type State = PauliSum<ppvm_runtime::config::indexmap::ByteFxHashF64<4, CoefficientThreshold>>;

fn main() {
    let mut state: State = PauliSum::builder()
        .n_qubits(20)
        .strategy(CoefficientThreshold(1e-6))
        .capacity(400)
        .build();

    state += ("ZZIIIIIIIIIIIIIIIIII", 1.0);

    // Textbook H(0); CNOT(0, 1) → reversed for Heisenberg propagation.
    state.cnot(0, 1);
    state.h(0);
    state.truncate();

    let zero_state: PauliPattern = "Z?*".into();
    println!("trace = {}", state.trace(&zero_state));
}
