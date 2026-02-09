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

    const COS_PI_OVER_8: f64 = 0.9238795325112867; // cos(pi/8)
    const SIN_PI_OVER_8: f64 = 0.3826834323650898; // sin(pi/8)
    let expected_coefficients = vec![
        Complex {
            re: COS_PI_OVER_8 * COS_PI_OVER_8 - SIN_PI_OVER_8 * SIN_PI_OVER_8,
            im: 0.0,
        },
        Complex {
            re: 0.0,
            im: -2.0 * SIN_PI_OVER_8 * COS_PI_OVER_8,
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
