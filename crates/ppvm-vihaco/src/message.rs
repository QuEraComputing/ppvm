use smallvec::SmallVec;
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
    QubitBatch(SmallVec<[usize; 8]>),              // X, Y, Z, ...
    QubitBatchAndFloat(SmallVec<[usize; 8]>, f64), // RX, depolarize, ...
    TwoQubitBatch(SmallVec<[(usize, usize); 8]>),  // CX, CZ
    TwoQubitBatchAndFloat(SmallVec<[(usize, usize); 8]>, f64), // RXX, ...
    QubitBatchU3(SmallVec<[usize; 8]>, f64, f64, f64), // U3
    QubitBatchAndFloatArr3(SmallVec<[usize; 8]>, [f64; 3]), // PauliError
    TwoQubitBatchAndFloatArr3(SmallVec<[(usize, usize); 8]>, [f64; 3]), // Correlated loss
    TwoQubitBatchAndFloatArr15(SmallVec<[(usize, usize); 8]>, [f64; 15]), // TwoQubitPauliError
}
