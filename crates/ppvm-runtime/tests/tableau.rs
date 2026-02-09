use itertools::Itertools;
use num::complex::Complex64;
use ppvm_runtime::{config::dashmap::ByteFxHashF64, prelude::*};

#[test]
fn test_tableau() {
    // let conf =
    let mut tableau: Tableau<2, ByteFxHashF64<1>> = Tableau::new();

    tableau.h(0);
    tableau.cnot(0, 1);

    println!("{}", tableau);
}

#[test]
fn generalized_tableau() {
    let mut tableau: GeneralizedTableau<2, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);
    tableau.cnot(0, 1);
    tableau.t(0);

    assert_eq!(tableau.coefficients.len(), 2);
    let idx: Vec<_> = tableau
        .coefficients
        .clone()
        .into_iter()
        .map(|(_, i)| i)
        .sorted()
        .collect();
    assert_eq!(idx, vec![0, 1]);

    tableau.t_adj(0);

    println!("{}", tableau);

    assert_eq!(tableau.coefficients.len(), 1);
}
