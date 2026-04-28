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
    #[mnemonic = "cnot"]
    CNOT,
    #[mnemonic = "cz"]
    CZ,

    // T gate
    T,
    TAdj,

    // Single-qubit rotations
    #[mnemonic = "rx"]
    RX,
    #[mnemonic = "ry"]
    RY,
    #[mnemonic = "rz"]
    RZ,

    // Two-qubit rotations
    #[mnemonic = "rxx"]
    RXX,
    #[mnemonic = "ryy"]
    RYY,
    #[mnemonic = "rzz"]
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
