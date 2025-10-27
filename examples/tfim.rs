use ppvm_runtime::prelude::*;
use ppvm_sym::*;

fn main() {
    let mut sum = PauliSum::<config::fxhash::Byte<2, Term>>::builder()
        .n_qubits(4)
        .build();
    sum += ("ZZII", Term::from(1.0));
    for i in 0..=3 {
        sum.rx(i, Term::var(0));
    }

    for i in 0..3 {
        sum.rzz(i, i + 1, Term::var(0));
    }

    println!("{:?}", sum);
}
