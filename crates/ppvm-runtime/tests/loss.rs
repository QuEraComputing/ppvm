use ppvm_runtime::prelude::*;

#[test]
fn test_lossy_sum() {
    let mut state: PauliSum<
        config::indexmap::ByteFxHashF64<
            4,
            NoStrategy,
            LossyPauliWord<[u8; 4], fxhash::FxBuildHasher>,
        >,
    > = PauliSum::<
        config::indexmap::ByteFxHashF64<4, _, LossyPauliWord<[u8; 4], fxhash::FxBuildHasher>>,
    >::builder()
    .n_qubits(2)
    .build();

    state += ("ZI", 1.0);
    state += ("IZ", 1.0);
    state += ("LL", 1.0);

    println!("{}", state);

    assert_eq!(state.data().len(), 3);
}
