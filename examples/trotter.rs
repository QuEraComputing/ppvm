use ppvm_runtime::{prelude::*, strategy::CoefficientThreshold};

/*
Implementation of the XZZ Ising chain Hamiltonian.
See here https://qiskit-community.github.io/qiskit-algorithms/tutorials/13_trotterQRTE.html
and here https://bloqade.quera.com/latest/digital/tutorials/circuits_with_bloqade/#squin-kernel-statements
 */

// type alias
type State = PauliSum<config::indexmap::ByteFxHashF64<4, CoefficientThreshold>>;

fn main() {
    // parameters
    let n = 20;
    let h = 1.0;
    let dt = 0.1 / h;
    let time = 1.0 / h;
    let j = 1.0 / 8.0 * h;

    let strat = CoefficientThreshold(1e-6);
    let mut state: State = PauliSum::builder()
        .n_qubits(n)
        .strategy(strat)
        .capacity(n.pow(2)) // NOTE: the capacity setting has a big impact on performance
        .build();

    // initial state: let's calculate the expectation value of Sum(Z(i))
    let tmp = vec!['I'; n];
    for i in 0..n {
        let mut zi = tmp.clone();
        zi[i] = 'Z';
        let zi_string: String = zi.iter().collect();
        state += (zi_string, 1.0);
    }
    println!("Initial state: {}", state);

    let noise_params = [1e-4; 3];
    let results = trotter(&mut state, n, time, dt, j, h, noise_params);

    // print some output
    println!(
        "Paramters: [n={}, h={}, J={}, dt={}, time={}]",
        n, h, j, dt, time
    );
    // println!("Final state: {}", state);
    println!("Results: {:?}", results);

    // check to see if we should have used a MaxPauliWeight strategy as well
    let max_weight = state.data().iter().map(|(k, _)| k.weight()).max().unwrap();
    println!("Maximum weight encountered: {}", max_weight);
}

fn trotter(
    state: &mut State,
    n: usize,
    total_time: f64,
    dt: f64,
    interaction_strength: f64,
    external_field: f64,
    noise_params: [f64; 3],
) -> Vec<f64> {
    let steps = (total_time / dt) as usize;
    let zero_state_pattern: PauliPattern = "Z?*".into();
    let mut expectation_values = Vec::<f64>::with_capacity(steps);

    let theta_zz = dt * interaction_strength;
    let theta_x = dt * external_field;
    for _ in 0..steps {
        expectation_values.push(state.trace(&zero_state_pattern));

        // perform trotter step
        for i in 0..n {
            state.rx(i, theta_x);
            state.pauli_error(i, noise_params);
        }
        for i in 0..n - 1 {
            state.rzz(i, i + 1, theta_zz);
            state.pauli_error(i, noise_params);
            state.pauli_error(i + 1, noise_params);
        }

        // truncate state in each iteration
        state.truncate();
    }

    // add final state expectation value
    expectation_values.push(state.trace(&zero_state_pattern));

    expectation_values
}
