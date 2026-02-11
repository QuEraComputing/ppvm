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

    println!("{}", tableau);

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
    tableau.x(0);
    tableau.t(0);

    let mut sorted_coefficients = tableau.coefficients.clone();
    sorted_coefficients.sort_by(|entry1, entry2| entry1.1.cmp(&entry2.1));

    let expected_coefficients = vec![
        Complex {
            re: 0.8535533905932737,
            im: 0.3535533905932738,
        },
        Complex {
            re: -0.14644660940672624,
            im: 0.3535533905932738,
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
fn test_generalized_tableau_multiple_ts() {
    let mut tableau: GeneralizedTableau<1, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);

    tableau.t(0);
    tableau.t(0);
    tableau.t(0);
    tableau.t(0);

    // four T gates should be equivalent to a Z
    assert_eq!(tableau.coefficients.len(), 1);
}

#[test]
fn test_generalized_tableau_multiple_ts2() {
    let mut tableau: GeneralizedTableau<2, ByteFxHashF64<1>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);
    tableau.h(1);

    tableau.t(0);
    tableau.t(0);
    tableau.t(0);
    tableau.t(0);

    tableau.t(1);
    tableau.t(1);
    tableau.t(1);
    tableau.t(1);

    // four T gates should be equivalent to a Z
    assert_eq!(tableau.coefficients.len(), 1);
}

#[test]
fn test_generalized_tableau_multiqubit_branching() {
    const N: usize = 18;
    let mut tableau: GeneralizedTableau<N, ByteFxHashF64<3>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    for i in 0..N {
        tableau.h(i);
    }

    // make sure to branch, but watch out since we have 2 ^ t scaling
    let mut tgate_counter: u32 = 0;
    for i in (0..10).step_by(2) {
        tableau.t(i);
        tgate_counter += 1;
    }

    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter));

    // test random measurement
    let outcome = tableau.measure(0);

    // should remove a branch
    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter - 1));

    // let's move it back
    if outcome {
        tableau.x(0);
    }

    tableau.h(0);
    tableau.t(0);
    assert_eq!(tableau.coefficients.len(), 2_usize.pow(tgate_counter));
}

#[test]
fn test_multiqubit_ghz_state() {
    const N: usize = 18;
    let mut tableau: GeneralizedTableau<N, ByteFxHashF64<3>, Vec<(Complex64, usize)>> =
        GeneralizedTableau::new(1e-12);

    tableau.h(0);
    tableau.t(0);
    // Let's generate a GHZ state
    for i in 0..N - 1 {
        tableau.cnot(i, i + 1);
    }

    assert_eq!(tableau.coefficients.len(), 2);

    let outcome = tableau.measure(0);
    println!("{}", tableau);
    println!("{}", tableau.coefficients.len());

    for i in 0..N {
        let outcome_i = tableau.measure(i);
        assert_eq!(outcome, outcome_i)
    }
}
