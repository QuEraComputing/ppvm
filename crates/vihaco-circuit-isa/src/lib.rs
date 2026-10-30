// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use vihaco::dialect;

#[dialect]
pub mod circuit {
    #[instruction]
    pub struct TwoQubitPauliError;

    #[instruction]
    pub struct Truncate;

    #[instruction]
    pub struct Trace;

    #[instruction]
    pub struct X;

    #[instruction]
    pub struct Y;

    #[instruction]
    pub struct Z;

    #[instruction]
    pub struct H;

    #[pattern = "'sqrt_x_adj"]
    #[instruction]
    pub struct SqrtXAdj;

    #[pattern = "'sqrt_x"]
    #[instruction]
    pub struct SqrtX;

    #[pattern = "'sqrt_y_adj"]
    #[instruction]
    pub struct SqrtYAdj;

    #[pattern = "'sqrt_y"]
    #[instruction]
    pub struct SqrtY;

    #[pattern = "'s_adj"]
    #[instruction]
    pub struct SAdj;

    #[instruction]
    pub struct S;

    #[instruction]
    pub struct CNOT;

    #[instruction]
    pub struct CZ;

    #[pattern = "'t_adj"]
    #[instruction]
    pub struct TAdj;

    #[instruction]
    pub struct T;

    #[instruction]
    pub struct RXX;

    #[instruction]
    pub struct RYY;

    #[instruction]
    pub struct RZZ;

    #[instruction]
    pub struct RX;

    #[instruction]
    pub struct RY;

    #[instruction]
    pub struct RZ;

    #[instruction]
    pub struct U3;

    #[instruction]
    pub struct Measure;

    #[instruction]
    pub struct Reset;

    #[instruction]
    pub struct R;

    #[instruction]
    pub struct Loss;

    #[instruction]
    pub struct CorrelatedLoss;

    #[instruction]
    pub struct PauliError;

    #[instruction]
    pub struct Depolarize2;

    #[instruction]
    pub struct Depolarize;

    #[instruction]
    pub struct Leakage;
}

#[cfg(test)]
mod tests {
    use super::circuit;
    use chumsky::Parser as _;
    use vihaco::Parse;

    // ─── Parse: tokens are the lowercased variant name ────────────────────

    #[test]
    fn parses_simple_lowercase_tokens() {
        assert_eq!(
            circuit::X::parser().parse("x").into_result(),
            Ok(circuit::X)
        );
        assert_eq!(
            circuit::Y::parser().parse("y").into_result(),
            Ok(circuit::Y)
        );
        assert_eq!(
            circuit::Z::parser().parse("z").into_result(),
            Ok(circuit::Z)
        );
        assert_eq!(
            circuit::H::parser().parse("h").into_result(),
            Ok(circuit::H)
        );
        assert_eq!(
            circuit::CNOT::parser().parse("cnot").into_result(),
            Ok(circuit::CNOT)
        );
        assert_eq!(
            circuit::CZ::parser().parse("cz").into_result(),
            Ok(circuit::CZ)
        );
        assert_eq!(
            circuit::U3::parser().parse("u3").into_result(),
            Ok(circuit::U3)
        );
        assert_eq!(
            circuit::Measure::parser().parse("measure").into_result(),
            Ok(circuit::Measure)
        );
        assert_eq!(
            circuit::Reset::parser().parse("reset").into_result(),
            Ok(circuit::Reset)
        );
        assert_eq!(
            circuit::R::parser().parse("r").into_result(),
            Ok(circuit::R)
        );
        assert_eq!(
            circuit::RX::parser().parse("rx").into_result(),
            Ok(circuit::RX)
        );
        assert_eq!(
            circuit::RY::parser().parse("ry").into_result(),
            Ok(circuit::RY)
        );
        assert_eq!(
            circuit::RZ::parser().parse("rz").into_result(),
            Ok(circuit::RZ)
        );
        assert_eq!(
            circuit::RXX::parser().parse("rxx").into_result(),
            Ok(circuit::RXX)
        );
        assert_eq!(
            circuit::Depolarize2::parser()
                .parse("depolarize2")
                .into_result(),
            Ok(circuit::Depolarize2)
        );
        assert_eq!(
            circuit::Depolarize::parser()
                .parse("depolarize")
                .into_result(),
            Ok(circuit::Depolarize)
        );
        assert_eq!(
            circuit::Loss::parser().parse("loss").into_result(),
            Ok(circuit::Loss)
        );
        assert_eq!(
            circuit::Leakage::parser().parse("leakage").into_result(),
            Ok(circuit::Leakage)
        );
    }
}
