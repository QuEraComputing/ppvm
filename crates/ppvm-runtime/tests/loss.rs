use ppvm_runtime::{prelude::*, strategy::MaxLossWeight};

type LossyPauliSum = PauliSum<
    config::indexmap::ByteFxHashF64<1, NoStrategy, LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>>,
>;

#[test]
fn test_lossy_sum() {
    let mut state: LossyPauliSum = LossyPauliSum::builder().n_qubits(2).build();

    state += ("ZI", 1.0);
    state += ("IZ", 1.0);
    state += ("LL", 1.0);

    println!("{}", state);

    assert_eq!(state.data().len(), 3);
}

#[test]
fn test_reset_channel() {
    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("X", 1.0);
    let state2 = state.clone();
    state.reset_loss_channel(0);
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("Y", 1.0);
    let state2 = state.clone();
    state.reset_loss_channel(0);
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("I", 1.0);
    let mut state2 = state.clone();
    state.reset_loss_channel(0);
    state2 += ("L", 1.0);
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("Z", 1.0);
    let mut state2 = state.clone();
    state.reset_loss_channel(0);
    state2 += ("L", 1.0);
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("L", 1.0);
    let mut state2 = state.clone();
    state.reset_loss_channel(0);
    state2 *= 0.0;
    assert_eq!(state, state2);
}

#[test]
fn test_loss_channel() {
    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("X", 1.0);
    let mut state2 = state.clone();
    state.loss_channel(0, 0.2);
    state2 *= 0.8;
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("Y", 1.0);
    let mut state2 = state.clone();
    state.loss_channel(0, 0.2);
    state2 *= 0.8;
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("I", 1.0);
    let mut state2 = state.clone();
    state.loss_channel(0, 0.2);
    state2 *= 0.8;
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("Z", 1.0);
    let mut state2 = state.clone();
    state.loss_channel(0, 0.2);
    state2 *= 0.8;
    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("L", 1.0);
    let mut state2 = state.clone();
    state.loss_channel(0, 0.2);
    state2 += ("I", 0.2);
    assert_eq!(state, state2);
}

#[test]
fn test_single_qubit_loss() {
    let mut state = LossyPauliSum::builder().n_qubits(1).build();

    state += ("Z", 1.0);

    state.reset_loss_channel(0);

    let intermediate = state.clone();

    // apply identity
    state.x(0);
    state.x(0);

    assert_eq!(state, intermediate);

    // overall circuit is X + loss
    state.loss_channel(0, 0.1);
    state.x(0);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);
    assert!((overlap + 0.8).abs() < 1e-10);
}

#[test]
fn test_ghz_final_loss() {
    // GHZ state circuit, with loss channels at the end, causing uncorrelated ZZ
    // expectation values some of the time.
    let mut state = LossyPauliSum::builder().n_qubits(2).build();

    let p_l = 0.1;

    state += ("ZZ", 1.0);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    // Applying some identity gates shouldn't affect loss
    state.x(0);
    state.x(1);
    state.x(0);
    state.x(1);

    // just lose at the end, before this we should have a perfect GHZ state
    state.loss_channel(0, p_l);
    state.loss_channel(1, p_l);

    state.cnot(0, 1);
    state.h(0);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);

    // Compute the expected <ZZ>: half the time it's |00> anyway, |11> state is affected by loss, but reset to 0
    let prob = 0.5 + 0.5 * ((1.0 - p_l) * (1.0 - p_l) - 2.0 * p_l * (1.0 - p_l) + p_l * p_l);

    assert!((overlap - prob).abs() < 1e-10);
}

#[test]
fn test_ghz() {
    let mut state = LossyPauliSum::builder().n_qubits(2).build();

    let p_l = 0.1;

    state += ("ZZ", 1.0);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    state.loss_channel(0, p_l);
    state.loss_channel(1, p_l);

    state.cnot(0, 1);

    state.loss_channel(0, 2.0 * p_l);
    state.h(0);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);

    // in 2p_l cases, the first qubit is lost after hadamard
    let prob = 2.0 * p_l
        + (1.0 - 2.0 * p_l)
            * (0.5 + 0.5 * ((1.0 - p_l) * (1.0 - p_l) - 2.0 * p_l * (1.0 - p_l) + p_l * p_l));
    assert!((overlap - prob).abs() < 1e-10);
}

#[test]
fn test_loss_truncation() {
    let strat = ppvm_runtime::strategy::MaxLossWeight(2);
    let mut state: PauliSum<
        config::indexmap::ByteFxHashF64<
            1,
            MaxLossWeight,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    > = PauliSum::<
        config::indexmap::ByteFxHashF64<
            1,
            MaxLossWeight,
            LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>,
        >,
    >::builder()
    .n_qubits(3)
    .strategy(strat)
    .build();

    state += ("ZZZ", 1.0);
    state.reset_loss_channel(0);
    state.reset_loss_channel(1);
    state.reset_loss_channel(2);

    state.loss_channel(0, 0.1);
    state.loss_channel(1, 0.1);
    state.loss_channel(2, 0.1);

    let original_len = state.data().len();

    state.truncate();
    assert_eq!(state.data().len(), original_len - 1);
}

#[test]
fn test_ghz_final_correlated_loss() {
    // GHZ state circuit, with loss channels at the end, causing uncorrelated ZZ
    // expectation values some of the time.
    let mut state = LossyPauliSum::builder().n_qubits(2).build();

    let p_l = 0.1;

    state += ("ZZ", 1.0);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    // Applying some identity gates shouldn't affect loss
    state.x(0);
    state.x(1);
    state.x(0);
    state.x(1);

    // just lose at the end, before this we should have a perfect GHZ state
    state.correlated_loss_channel(0, 1, p_l);

    state.cnot(0, 1);
    state.h(0);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);

    // correlated loss should leave ZZ invariant
    assert!((overlap - 1.0).abs() < 1e-10);

    // same thing again, but this time we flip one bit leading to |10> + |01>
    let mut state = LossyPauliSum::builder().n_qubits(2).build();

    let p_l = 0.1;

    state += ("ZZ", 1.0);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    state.x(0); // flip first qubit

    // just lose at the end, before this we should have a perfect GHZ state
    state.correlated_loss_channel(0, 1, p_l);

    state.cnot(0, 1);
    state.h(0);

    println!("{}", state);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);

    // in 1 - p cases we have the state invarian (ZZ = -1), in p cases we have |00>
    let expected_overlap = -1.0 * (1.0 - p_l) + 1.0 * p_l;

    // correlated loss should leave ZZ invariant
    assert!((overlap - expected_overlap).abs() < 1e-10);
}
