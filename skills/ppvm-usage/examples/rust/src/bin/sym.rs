//! Skill code block: symbolic Pauli propagation.
//!
//! Mirrors the §7 "ppvm-sym — symbolic propagation" example in
//! `skills/ppvm-usage/SKILL.md`. We propagate a Pauli observable
//! through a small parametric circuit with `Term`-valued coefficients,
//! then evaluate the resulting expression at a concrete pair of
//! angles to check that the API surface still matches what the skill
//! tells agents to write.

use ppvm_pauli_sum::prelude::*;
use ppvm_sym::Term;

fn main() {
    let mut sum = PauliSum::<ppvm_pauli_sum::config::fxhash::Byte<2, Term>>::builder()
        .n_qubits(2)
        .build();
    sum += ("ZZ", Term::from(1.0));

    sum.rz(0, Term::var(0));
    sum.ry(0, Term::var(1));
    sum.cnot(0, 1);
    sum.rx(1, Term::var(1));

    let pat: PauliPattern = "Z?*".into();
    let trace = sum.trace(&pat);

    // Substitute concrete angles. The skill claims this finishing step
    // works — assert it does and that the value is finite.
    let value = trace.eval(&[0.3_f64, 0.5_f64]).expect("trace evaluation");
    assert!(value.is_finite(), "non-finite trace value: {value}");
    println!("symbolic_trace_finite=true value={value:.6}");
}
