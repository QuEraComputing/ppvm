use ppvm_runtime::{prelude::*, strategy::CoefficientThreshold};

#[test]
fn test_coefficient_threshold() {
    let strat = CoefficientThreshold(1e-4);
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>> =
        PauliSum::builder().n_qubits(2).strategy(strat).build();

    state += ("ZZ", 1.0);
    state += ("XX", 1e-3);
    state += ("YY", 1e-5);
    state += ("XZ", -0.5);
    state.truncate();

    let mut target_state: PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>> =
        PauliSum::builder().strategy(strat).n_qubits(2).build();
    target_state += ("ZZ", 1.0);
    target_state += ("XX", 1e-3);
    target_state += ("XZ", -0.5);

    assert_eq!(state, target_state);
}
