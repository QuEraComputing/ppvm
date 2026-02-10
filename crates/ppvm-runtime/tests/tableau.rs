use itertools::Itertools;
use num::complex::{Complex, Complex64};
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

    assert_eq!(tableau.coefficients.len(), 1);

    tableau.t(0);
    tableau.t(1);

    // NOTE: since IZ|psi> = (IZ) * ZZ |psi> = ZI|psi>, we don't branch again
    assert_eq!(tableau.coefficients.len(), 2);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    const PI: f64 = std::f64::consts::PI;
    let cos_pi_8: f64 = (PI / 8.0).cos();
    let sin_pi_8: f64 = (PI / 8.0).sin();
    let expected_coefficients = vec![
        Complex {
            re: (PI / 4.0).cos() * (cos_pi_8 * cos_pi_8 - sin_pi_8 * sin_pi_8),
            im: (PI / 4.0).sin() * (cos_pi_8 * cos_pi_8 - sin_pi_8 * sin_pi_8),
        },
        Complex {
            re: (PI / 4.0).cos() * 2.0 * sin_pi_8 * cos_pi_8,
            im: (PI / 4.0).sin() * -2.0 * sin_pi_8 * cos_pi_8,
        },
    ];

    for ((val1, idx1), (idx2, val2)) in sorted_coefficients
        .iter()
        .zip(expected_coefficients.iter().enumerate())
    {
        assert_eq!(idx1, &idx2);
        assert!((val1.re - val2.re).abs() < 1e-11);
        assert!((val1.im - val2.im).abs() < 1e-11);
    }

    println!("{}", tableau);
}

#[test]
fn test_generalized_tableau_phase() {
    let mut tableau: GeneralizedTableau<1, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);
    tableau.t(0);
    tableau.t(0);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    let expected_coefficients = vec![Complex { re: 0.5, im: 0.5 }, Complex { re: 0.5, im: -0.5 }];

    for ((val1, idx1), (idx2, val2)) in sorted_coefficients
        .iter()
        .zip(expected_coefficients.iter().enumerate())
    {
        assert_eq!(idx1, &idx2);
        assert!((val1.re - val2.re).abs() < 1e-11);
        assert!((val1.im - val2.im).abs() < 1e-11);
    }

    let mut tableau: GeneralizedTableau<1, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);
    tableau.t(0);
    tableau.t(0);

    println!("{}", tableau);
}
