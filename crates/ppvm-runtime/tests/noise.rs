use ppvm_runtime::{prelude::*, strategy::CoefficientThreshold};
use std::f64::consts::PI;

#[test]
fn test_two_qubit_pauli_error() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("IZ", 1.0);

    let mut state2 = state.clone();

    let mut p = [0.0; 15];
    p[0] = 1.0;

    state.two_qubit_pauli_error(0, 1, p);
    state2.x(1);
    assert_eq!(state, state2);

    p[0] = 0.0;
    p[1] = 1.0;
    state.two_qubit_pauli_error(0, 1, p);
    state2.y(1);
    assert_eq!(state, state2);

    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("IX", 1.0);

    let mut state2 = state.clone();

    p[1] = 0.0;
    p[2] = 1.0;
    state.two_qubit_pauli_error(0, 1, p);
    state2.z(1);
    assert_eq!(state, state2);

    p[2] = 0.0;

    // {IX, IY, IZ, XI, XX, XY, XZ, YI, YX, YY, YZ, ZI, ZX, ZY, ZZ}
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("YI", 1.0);
    let mut state2 = state.clone();

    p[4] = 1.0;
    state.two_qubit_pauli_error(0, 1, p);
    state2.rxx(0, 1, PI);
    state2.truncate();

    assert_eq!(state, state2);

    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("YZ", 1.0);
    let state2 = state.clone();

    p[4] = 1.0; // XX
    state.two_qubit_pauli_error(0, 1, p);

    assert_eq!(state, state2);

    p[4] = 0.0;
    p[9] = 1.0; // YY

    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("XI", 1.0);
    let mut state2 = state.clone();
    state2 *= -1.0;

    p[9] = 0.0;
    p[14] = 1.0; // ZZ
    state.two_qubit_pauli_error(0, 1, p);

    assert_eq!(state, state2);

    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("XI", 1.0);
    let mut state2 = state.clone();
    state2 *= -1.0;

    p[4] = 1.0; // XX
    state.two_qubit_pauli_error(0, 1, p);

    assert_eq!(state, state2);
}

#[test]
fn test_depolarizing_error() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(3).build();

    state += ("ZZZ", 1.0);

    let ps = [0.1, 0.2, 0.3];
    state.depolarize(0, ps[0]);
    state.depolarize(1, ps[1]);
    state.depolarize(2, ps[2]);

    println!("State: {}", state);
    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);
    println!("Overlap: {}", overlap);

    let result: f64 = ps.map(|p| 1.0 - 4.0 * p / 3.0).iter().product();

    assert!((overlap - result).abs() < 1e-10);
}
