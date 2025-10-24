use ppvm_runtime::prelude::*;

#[test]
fn test_h() {
    // test for H * Y * H -> -Y since tests on PauliWords don't track the phase
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(1).build();

    state += ("Y", 1.0);
    state.h(0);

    let mut len = 0;
    for (k, v) in state.data().iter() {
        len += 1;
        assert_eq!(v, &-1.0);
        assert_eq!(k.to_string(), "Y");
    }
    assert_eq!(len, 1);
}

#[test]
fn test_ghz_forward() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();

    // build initial state as explicit 00 state
    state += ("ZI", 1.0);
    state += ("IZ", 1.0);
    state += ("ZZ", 1.0);
    state += ("II", 1.0);

    // GHZ prep
    state.h(0);
    state.cnot(0, 1);

    // explicitly construct GHZ state
    let mut ghz_state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();

    ghz_state += ("XX", 1.0);
    ghz_state += ("YY", -1.0);
    ghz_state += ("ZZ", 1.0);
    ghz_state += ("II", 1.0);

    println!("{}", state);
    println!("{}", ghz_state);
    assert_eq!(state, ghz_state);
}

#[test]
fn test_ghz_backward() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();

    // prepare "state" as the final expectation value we want
    state += ("ZZ", 1.0);

    // propagate through the circuit in backwards order
    state.cnot(0, 1);
    state.h(0);

    // zero state
    let zero_state: PauliPattern = "Z?*".into();
    let result = state.trace(&zero_state);
    state.data().trace(&zero_state);
    state.data().iter();

    println!("{:?}", state.data());

    assert_eq!(result, 1.0);
}

#[test]
fn test_pauli_sum_trace() {
    let mut state1: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();

    // prepare "state" as the final expectation value we want
    state1 += ("ZZ", 0.5);
    let state2 = state1.clone();
    let result = state1.trace(&state2);
    assert_eq!(result, 0.25);
}
