//! Normalized AST. Built by `crate::normalize::to_tableau`. Every variant
//! is ready to dispatch directly to a `GeneralizedTableau` method.

#[derive(Debug, Clone, PartialEq)]
pub struct TableauProgram {
    pub instructions: Vec<Instruction>,
    /// Sum over `M` / `MZ` / `MR` target counts, multiplied by enclosing
    /// `REPEAT` counts. Used to pre-size shot result buffers.
    pub expected_measurement_count: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Instruction {
    Gate {
        kind: GateKind,
        targets: Vec<usize>,
        line: usize,
    },
    Noise {
        kind: NoiseKind,
        targets: Vec<usize>,
        args: Vec<f64>,
        line: usize,
    },
    Measure {
        kind: MeasureKind,
        targets: Vec<usize>,
        line: usize,
    },
    /// Phase-1 no-op; preserved so executor can track them for future tooling.
    Annotation {
        line: usize,
    },
    Repeat {
        count: u64,
        body: Vec<Instruction>,
        line: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GateKind {
    Reset,
    X, Y, Z, H, S, SDag, SqrtX, SqrtXDag, SqrtY, SqrtYDag,
    T, TDag,
    RX { theta: f64 },
    RY { theta: f64 },
    RZ { theta: f64 },
    U3 { theta: f64, phi: f64, lambda: f64 },
    CX, CY, CZ,
}

#[derive(Debug, Clone, PartialEq)]
pub enum NoiseKind {
    Depolarize1,
    Depolarize2,
    PauliChannel1,
    PauliChannel2,
    XError,
    YError,
    ZError,
    /// `I_ERROR[loss]` — single-qubit loss with probability `args[0]`.
    Loss,
    /// `I_ERROR[correlated_loss]` — pairwise loss with `args = [p_x, p_y, p_z]`
    /// (or `[p]` shorthand expanded to `[p, 0, 0]` by the normalizer).
    CorrelatedLoss,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeasureKind {
    M,
    MR,
}
