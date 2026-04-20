use vihaco::Message;

#[derive(Debug, Clone, Copy, Message)]
pub enum CircuitMessage {
    Qubit(usize),                                   // X, Y, Z, ...
    QubitAndFloat(usize, f64),                      // RX, depolarize, ...
    TwoQubit(usize, usize),                         // CX, CZ
    TwoQubitAndFloat(usize, usize, f64),            // RXX, ...
    QubitU3(usize, f64, f64, f64),                  // U3
    QubitAndFloatArr3(usize, [f64; 3]),             // PauliError
    TwoQubitAndFloatArr3(usize, usize, [f64; 3]),   // Correlated loss
    TwoQubitAndFloatArr15(usize, usize, [f64; 15]), // TwoQubitPauliError
}
