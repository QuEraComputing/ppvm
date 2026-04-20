use vihaco::Instruction;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Instruction)]
pub enum CircuitInstruction {
    // Single-Qubit Clifford gates
    X,
    Y,
    Z,
    H,
    S,
    SAdj,
    SqrtX,
    SqrtY,
    SqrtXAdj,
    SqrtYAdj,

    // Controlled gates
    CNOT,
    CZ,

    // T gate
    T,
    TAdj,

    // Single-qubit rotations
    RX,
    RY,
    RZ,

    // Two-qubit rotations
    RXX,
    RYY,
    RZZ,

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
    TwoQubitPauliError,
    Depolarize,
    Depolarize2,
}
