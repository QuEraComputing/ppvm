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

    // Measurement & Reset
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

impl std::fmt::Display for CircuitInstruction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use CircuitInstruction::*;
        match self {
            TwoQubitPauliError => write!(f, "TwoQubitPauliError"),

            X => write!(f, "X"),
            Y => write!(f, "Y"),
            Z => write!(f, "Z"),
            H => write!(f, "H"),

            SqrtXAdj => write!(f, "SqrtXAdj"),

            SqrtX => write!(f, "SqrtX"),

            SqrtYAdj => write!(f, "SqrtYAdj"),

            SqrtY => write!(f, "SqrtY"),

            SAdj => write!(f, "SAdj"),
            S => write!(f, "S"),

            CNOT => write!(f, "CNOT"),
            CZ => write!(f, "CZ"),

            TAdj => write!(f, "TAdj"),
            T => write!(f, "T"),

            RXX => write!(f, "RXX"),
            RYY => write!(f, "RYY"),
            RZZ => write!(f, "RZZ"),

            RX => write!(f, "RX"),
            RY => write!(f, "RY"),
            RZ => write!(f, "RZ"),

            U3 => write!(f, "U3"),

            Measure => write!(f, "Measure"),
            Reset => write!(f, "Reset"),

            Loss => write!(f, "Loss"),
            CorrelatedLoss => write!(f, "CorrelatedLoss"),

            PauliError => write!(f, "PauliError"),
            Depolarize2 => write!(f, "Depolarize2"),
            Depolarize => write!(f, "Depolarize"),
        }
    }
}
