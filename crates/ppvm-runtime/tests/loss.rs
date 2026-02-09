use ppvm_runtime::{prelude::*, strategy::MaxLossWeight};
use std::panic::{AssertUnwindSafe, catch_unwind};

type LossyPauliSum = PauliSum<
    config::indexmap::ByteFxHashF64<1, NoStrategy, LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>>,
>;
type LossyPauliSumHashMap = PauliSum<
    config::fxhash::Byte<1, f64, NoStrategy, LossyPauliWord<[u8; 1], fxhash::FxBuildHasher>>,
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
fn test_reset_loss_channel_accumulates_duplicate_target_indexmap() {
    let mut state = LossyPauliSum::builder().n_qubits(1).build();
    state += ("I", 2.0);
    state += ("Z", 3.0);

    state.reset_loss_channel(0);

    let i: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "I".into();
    let z: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "Z".into();
    let l: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "L".into();
    assert!(state.contains(&i, &2.0));
    assert!(state.contains(&z, &3.0));
    assert!(state.contains(&l, &5.0));
}

#[test]
fn test_reset_loss_channel_accumulates_duplicate_target_hashmap() {
    let mut state = LossyPauliSumHashMap::builder().n_qubits(1).build();
    state += ("I", 2.0);
    state += ("Z", 3.0);

    state.reset_loss_channel(0);

    let i: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "I".into();
    let z: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "Z".into();
    let l: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "L".into();
    assert!(state.contains(&i, &2.0));
    assert!(state.contains(&z, &3.0));
    assert!(state.contains(&l, &5.0));
}

#[test]
fn test_rx_on_lost_qubit_is_noop_and_does_not_panic() {
    let mut state = LossyPauliSum::builder().n_qubits(1).build();
    state += ("L", 1.0);

    let result = catch_unwind(AssertUnwindSafe(|| {
        state.rx(0, 0.3);
    }));
    assert!(result.is_ok());

    let l: LossyPauliWord<[u8; 1], fxhash::FxBuildHasher> = "L".into();
    assert_eq!(state.data().len(), 1);
    assert!(state.contains(&l, &1.0));
}

#[test]
fn test_rxx_with_loss_is_noop_and_does_not_panic() {
    let mut state = LossyPauliSum::builder().n_qubits(2).build();
    state += ("LZ", 1.0);
    let mut state2 = state.clone();

    let result = catch_unwind(AssertUnwindSafe(|| {
        state.rxx(0, 1, 0.3);
    }));
    assert!(result.is_ok());

    state2.rx(1, 0.3);

    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(2).build();
    state += ("ZL", 1.0);
    let mut state2 = state.clone();

    let result = catch_unwind(AssertUnwindSafe(|| {
        state.rxx(0, 1, 0.3);
    }));
    assert!(result.is_ok());

    state2.rx(0, 0.3);

    assert_eq!(state, state2);

    let mut state = LossyPauliSum::builder().n_qubits(2).build();
    state += ("LL", 1.0);
    let state2 = state.clone();

    let result = catch_unwind(AssertUnwindSafe(|| {
        state.rxx(0, 1, 0.3);
    }));
    assert!(result.is_ok());

    assert_eq!(state, state2);
}
