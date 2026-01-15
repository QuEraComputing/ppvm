use ppvm_runtime::prelude::*;

#[test]
fn test_loss_channel() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(1).build();
    state += ("I", 1.0);
    state.loss_channel(0, 0.3);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);
    assert!((overlap - 0.7).abs() < 1e-10);
}


#[test]
fn test_ghz_loss_channel() {
    let mut state: PauliSum<config::indexmap::ByteFxHashF64<4>> =
        PauliSum::builder().n_qubits(2).build();

    // build initial state as explicit 00 state
    state += ("ZZ", 1.0);

    // GHZ prep
    state.loss_channel(0, 0.2);
    state.loss_channel(1, 0.2);
    state.cnot(0, 1);

    state.loss_channel(0, 0.3);
    state.h(0);
    // apply loss channel to both qubits

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);

    println!("Calculated overlap: {}", overlap);

    assert!((overlap - 0.44).abs() < 1e-2);
    /* NOTE: number obtained from state-vector simulation in bloqade-circuit
    This is equivalent to counting any shot with loss as zero amplitude state, i.e.,
    zero contribution to the expectation value.

    ```python
    @squin.kernel
    def ghz():
        q = squin.qalloc(2)

        squin.h(q[0])
        pl1 = 0.3
        squin.qubit_loss(pl1, q[0])

        squin.cx(q[0], q[1])
        pl2 = 0.2
        squin.qubit_loss(pl2, q[0])
        squin.qubit_loss(pl2, q[1])

        return q


    sim = StackMemorySimulator(min_qubits=2, loss_m_result=Measurement.Zero)
    nshots = 2000

    Z = np.array([[1, 0], [0, -1]])
    ZZ = np.kron(Z, Z)

    avg = 0.0
    for _ in range(nshots):
        res = sim.run(ghz)

        if any(not q.is_active() for q in res):
            # NOTE: count shots with loss as zero
            continue

        ket = res[0].sim_reg.out_ket()
        exp = np.vdot(ket, ZZ @ ket)
        avg += exp.real / nshots

    print(avg)
    ```
    */
}