use vihaco::Instruction;
use vihaco_parser::Parse;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Instruction, Parse)]
pub enum CircuitInstruction {
    // NOTE: longer tokens need to go first
    TwoQubitPauliError, // needs to go before T

    // Single-Qubit Clifford gates
    X,
    Y,
    Z,
    H,

    #[token = "sqrt_x_adj"]
    SqrtXAdj,

    #[token = "sqrt_x"]
    SqrtX,

    #[token = "sqrt_y_adj"]
    SqrtYAdj,

    #[token = "sqrt_y"]
    SqrtY,

    #[token = "s_adj"]
    SAdj,
    S,

    // Controlled gates
    CNOT,
    CZ,

    // T gate
    TAdj,
    T,

    // Two-qubit rotations
    RXX,
    RYY,
    RZZ,

    // Single-qubit rotations
    RX,
    RY,
    RZ,

    // U3
    U3,

    // Measureme & Reset
    Measure,
    Reset,

    // Loss
    Loss,
    CorrelatedLoss,

    // Noise
    PauliError,
    Depolarize2,
    Depolarize,
}
