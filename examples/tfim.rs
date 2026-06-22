use ppvm_pauli_sum::prelude::*;
use ppvm_sym::*;

fn main() {
    let mut o_t = PauliSum::<config::fxhash::Byte<2, Term>>::builder()
        .n_qubits(4)
        .build();
    o_t += ("ZZII", Term::from(1.0));
    for i in 0..=3 {
        o_t.rx(i, Term::var(0));
    }

    for i in 0..3 {
        o_t.rzz([i, i + 1], Term::var(0));
    }

    let mut b_t = PauliSum::<config::fxhash::Byte<2, Term>>::builder()
        .n_qubits(4)
        .build();
    b_t += ("ZIII", Term::from(1.0));
    for i in 0..=3 {
        b_t.rx(i, Term::var(0));
    }
    for i in 0..3 {
        b_t.rzz([i, i + 1], Term::var(0));
    }
}
