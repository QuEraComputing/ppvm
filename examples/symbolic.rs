use ppvm_runtime::prelude::*;
use ppvm_sym::*;

fn main() {
    let pat: PauliPattern = "Z?*".into();
    let mut sum = PauliSum::<config::fxhash::Byte<2, Term>>::builder()
        .n_qubits(2)
        .build();
    sum += ("ZZ", Term::from(1.0));

    sum.rz(0, Term::var(0));
    sum.ry(0, Term::var(1));
    sum.rz(0, Term::var(0));

    sum.rz(1, Term::var(0));
    sum.ry(1, Term::var(1));
    sum.rz(1, Term::var(0));

    sum.cnot([0, 1]);

    sum.rx(0, Term::var(1));
    sum.rx(1, Term::var(1));

    let trace = sum.trace(&pat);
    println!("Trace expression: {}", trace);
    let value = trace.eval(&[1.1, 2.1]).unwrap();
    println!("Trace: {}", value);
}
