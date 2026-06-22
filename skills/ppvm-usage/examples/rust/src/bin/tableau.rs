//! Skill code block: generalized stabilizer tableau in Rust.
//!
//! Mirrors the §5 "ppvm-tableau" / "Generalized stabilizer tableau"
//! example in `skills/ppvm-usage/SKILL.md`. CI compiles + runs this,
//! so a signature change to `GeneralizedTableau::new` or to the gate
//! trait surface will break the build instead of silently leaving the
//! skill stale.

use ppvm_pauli_sum::prelude::*;
use ppvm_tableau::prelude::*;

fn main() {
    // GeneralizedTableau::new takes (n_qubits, coefficient_threshold).
    // The third generic parameter (sparse-vector backing) defaults to
    // Vec<(Complex64, IndexType)>, so we leave it implicit.
    let mut tab: GeneralizedTableau<ppvm_pauli_sum::config::indexmap::ByteFxHashF64<1>, usize> =
        GeneralizedTableau::new(2, 1e-10);

    tab.h(0);
    tab.cnot(0, 1);

    let r0 = tab.measure(0);
    let r1 = tab.measure(1);

    // GHZ: outcomes are perfectly correlated on each shot.
    assert_eq!(r0, r1, "GHZ correlation broken: {r0:?} vs {r1:?}");
    println!("q0={r0:?} q1={r1:?} correlated=true");
}
