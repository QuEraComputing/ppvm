use ppvm_runtime::{config::dashmap::ByteFxHashF64, prelude::*};

#[test]
fn test_tableau() {
    // let conf =
    let mut tableau: Tableau<2, ByteFxHashF64<1>> = Tableau::new();

    tableau.h(0);
    tableau.cnot(0, 1);

    println!("{}", tableau);
}
