use vihaco::Message;

#[derive(Debug, Clone, Message)]
pub enum CircuitMessage {
    Qubit(usize),                                   // X, Y, Z, ...
    QubitAndFloat(usize, f64),                      // RX, depolarize, ...
    TwoQubit(usize, usize),                         // CX, CZ
    TwoQubitAndFloat(usize, usize, f64),            // RXX, ...
    QubitU3(usize, f64, f64, f64),                  // U3
    QubitAndFloatArr3(usize, [f64; 3]),             // PauliError
    TwoQubitAndFloatArr3(usize, usize, [f64; 3]),   // Correlated loss
    TwoQubitAndFloatArr15(usize, usize, [f64; 15]), // TwoQubitPauliError

    // batched instructions
    QubitBatch(Vec<usize>),                                     // X, Y, Z, ...
    QubitBatchAndFloat(Vec<usize>, f64),                        // RX, depolarize, ...
    TwoQubitBatch(Vec<(usize, usize)>),                         // CX, CZ
    TwoQubitBatchAndFloat(Vec<(usize, usize)>, f64),            // RXX, ...
    QubitBatchU3(Vec<usize>, f64, f64, f64),                    // U3
    QubitBatchAndFloatArr3(Vec<usize>, [f64; 3]),               // PauliError
    TwoQubitBatchAndFloatArr3(Vec<(usize, usize)>, [f64; 3]),   // Correlated loss
    TwoQubitBatchAndFloatArr15(Vec<(usize, usize)>, [f64; 15]), // TwoQubitPauliError
}
