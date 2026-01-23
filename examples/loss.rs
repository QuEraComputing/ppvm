use ppvm_runtime::prelude::*;

fn test_reset_channel() {
    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("X", 1.0);
    println!("Before {}", state);
    state.reset_loss_channel(0);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("Y", 1.0);
    println!("Before {}", state);
    state.reset_loss_channel(0);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("I", 1.0);
    println!("Before {}", state);
    state.reset_loss_channel(0);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("Z", 1.0);
    println!("Before {}", state);
    state.reset_loss_channel(0);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("L", 1.0);
    println!("Before {}", state);
    state.reset_loss_channel(0);
    println!("After {}", state);
}

fn test_loss_channel() {
    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("X", 1.0);
    println!("Before {}", state);
    state.loss_channel(0, 0.2);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("Y", 1.0);
    println!("Before {}", state);
    state.loss_channel(0, 0.2);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("I", 1.0);
    println!("Before {}", state);
    state.loss_channel(0, 0.2);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("Z", 1.0);
    println!("Before {}", state);
    state.loss_channel(0, 0.2);
    println!("After {}", state);

    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("L", 1.0);
    println!("Before {}", state);
    state.loss_channel(0, 0.2);
    println!("After {}", state);
}

fn test_single_qubit_loss() {
    let mut state = PauliSum::<config::fxhash::Byte<1, f64>>::builder()
        .n_qubits(1)
        .build();

    state += ("Z", 1.0);

    state.reset_loss_channel(0);
    println!("After reset loss channel: {}", state);

    state.x(0);
    state.x(0);

    println!("After two X: {}", state);

    state.loss_channel(0, 0.1);

    state.x(0);

    println!("Final state: {}", state);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);
    println!("Overlap with {}: {}", zero_pattern, overlap);
}

fn test_ghz() {
    let mut state = PauliSum::<config::fxhash::Byte<2, f64>>::builder()
        .n_qubits(2)
        .build();

    let p_l = 0.1;

    state += ("ZZ", 1.0);

    println!("Initial state: {}", state);

    state.reset_loss_channel(0);
    state.reset_loss_channel(1);

    // Applying some identity gates shouldn't affect loss
    state.x(0);
    state.x(1);
    state.x(0);
    state.x(1);

    state.loss_channel(0, p_l);
    state.loss_channel(1, p_l);

    println!("After loss channels: {}", state);
    // state.loss_channel(1, 0.1);

    state.cnot(0, 1);
    state.h(0);

    println!("Final state: {}", state);

    let zero_pattern: PauliPattern = "Z?*".into();
    let overlap = state.trace(&zero_pattern);
    println!("Overlap with {}: {}", zero_pattern, overlap);

    let prob = 0.5 + 0.5 * ((1.0 - p_l) * (1.0 - p_l) - 2.0 * p_l * (1.0 - p_l) + p_l * p_l);
    println!("Expected overlap: {}", prob);

    assert!((overlap - prob).abs() < 1e-10);
}

fn main() {
    // test_reset_channel();
    // test_loss_channel();
    // test_single_qubit_loss();
    test_ghz();
}
