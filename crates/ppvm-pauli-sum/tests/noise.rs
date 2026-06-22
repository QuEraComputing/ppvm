// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use ppvm_pauli_sum::{prelude::*, strategy::CoefficientThreshold};
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

#[test]
fn test_depolarize2() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();

    state += ("ZZ", 1.0);

    let p = 0.1_f64;
    state.depolarize2(0, 1, p);

    let pattern: PauliPattern = "Z0Z1".into();
    let overlap = state.trace(&pattern);

    // 8 of the 15 non-trivial two-qubit Paulis anticommute with ZZ,
    // so the coefficient scales by 1 - 2 * 8 * (p/15) = 1 - 16p/15.
    let expected = 1.0 - 16.0 * p / 15.0;
    assert!((overlap - expected).abs() < 1e-10);
}

#[test]
fn test_pauli_error() {
    let px = 0.05_f64;
    let py = 0.10_f64;
    let pz = 0.15_f64;
    let p = [px, py, pz];

    // I at addr0: coefficient should be unchanged
    let mut state_i: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(1).build();
    state_i += ("I", 1.0);
    let state_i_before = state_i.clone();
    state_i.pauli_error(0, p);
    assert_eq!(state_i, state_i_before);

    // X at addr0: scales by 1 - 2*py - 2*pz (Y and Z errors anticommute with X)
    let mut state_x: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(1).build();
    state_x += ("X", 1.0);
    state_x.pauli_error(0, p);
    let x_coeff = state_x.trace(&PauliPattern::from("X0"));
    let expected_x = 1.0 - 2.0 * py - 2.0 * pz;
    assert!((x_coeff - expected_x).abs() < 1e-10);

    // Y at addr0: scales by 1 - 2*px - 2*pz (X and Z errors anticommute with Y)
    let mut state_y: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(1).build();
    state_y += ("Y", 1.0);
    state_y.pauli_error(0, p);
    let y_coeff = state_y.trace(&PauliPattern::from("Y0"));
    let expected_y = 1.0 - 2.0 * px - 2.0 * pz;
    assert!((y_coeff - expected_y).abs() < 1e-10);

    // Z at addr0: scales by 1 - 2*px - 2*py (X and Y errors anticommute with Z)
    let mut state_z: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(1).build();
    state_z += ("Z", 1.0);
    state_z.pauli_error(0, p);
    let z_coeff = state_z.trace(&PauliPattern::from("Z0"));
    let expected_z = 1.0 - 2.0 * px - 2.0 * py;
    assert!((z_coeff - expected_z).abs() < 1e-10);

    // L at addr0: coefficient should be unchanged (lost qubit is skipped)
    type LossyPauliSum = PauliSum<
        config::indexmap::ByteFxHashF64<
            1,
            NoStrategy,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    >;
    let mut state_l: LossyPauliSum = LossyPauliSum::builder().n_qubits(1).build();
    state_l += ("L", 1.0);
    let state_l_before = state_l.clone();
    state_l.pauli_error(0, p);
    assert_eq!(state_l, state_l_before);
}

#[test]
fn test_amplitude_damping() {
    let gamma = 0.3_f64;

    // In the backwards (Heisenberg) propagation picture:
    //   E†[Z] = (1-γ)Z + γI   (longitudinal / T₁ decay)
    //   E†[X] = √(1-γ) X      (transverse / T₂ decay)
    //   E†[I] = I              (trace-preserving)

    // Longitudinal: Z decays and leaks into I
    let mut state_z: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();
    state_z += ("ZI", 1.0);
    state_z.amplitude_damping(0, gamma);
    let z_coeff = state_z.trace(&PauliPattern::from("Z0"));
    assert!((z_coeff - (1.0 - gamma)).abs() < 1e-10);

    // Transverse: X decays by √(1-γ)
    let mut state_x: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();
    state_x += ("XI", 1.0);
    state_x.amplitude_damping(0, gamma);
    let x_coeff = state_x.trace(&PauliPattern::from("X0"));
    assert!((x_coeff - (1.0 - gamma).sqrt()).abs() < 1e-10);

    // Physical cross-check: T₂ = 2T₁.
    // For pure amplitude damping, the transverse decay rate is half the
    // longitudinal, so the transverse factor squared equals the longitudinal
    // factor: (√(1-γ))² = (1-γ).
    // These are computed from separate arms of the channel, so this catches
    // any mismatch between the Z and X/Y scaling.
    assert!((x_coeff * x_coeff - z_coeff).abs() < 1e-10);

    // Trace-preserving: E†[I] = I, so the state should be entirely unchanged.
    let mut state_i: PauliSum<config::indexmap::ByteFxHashF64<1>> =
        PauliSum::builder().n_qubits(2).build();
    state_i += ("II", 1.0);
    let state_i_before = state_i.clone();
    state_i.amplitude_damping(0, gamma);
    assert_eq!(state_i, state_i_before);
}
