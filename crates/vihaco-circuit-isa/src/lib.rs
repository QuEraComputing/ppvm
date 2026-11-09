// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use vihaco::dialect;

pub mod messages;

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

    // Batched variants: same gate applied to every qubit (or pair) in the batch.

    #[pattern = "'x_many"]
    #[instruction]
    pub struct XMany;

    #[pattern = "'y_many"]
    #[instruction]
    pub struct YMany;

    #[pattern = "'z_many"]
    #[instruction]
    pub struct ZMany;

    #[pattern = "'h_many"]
    #[instruction]
    pub struct HMany;

    #[pattern = "'s_many"]
    #[instruction]
    pub struct SMany;

    #[pattern = "'s_adj_many"]
    #[instruction]
    pub struct SAdjMany;

    #[pattern = "'t_many"]
    #[instruction]
    pub struct TMany;

    #[pattern = "'t_adj_many"]
    #[instruction]
    pub struct TAdjMany;

    #[pattern = "'sqrt_x_many"]
    #[instruction]
    pub struct SqrtXMany;

    #[pattern = "'sqrt_x_adj_many"]
    #[instruction]
    pub struct SqrtXAdjMany;

    #[pattern = "'sqrt_y_many"]
    #[instruction]
    pub struct SqrtYMany;

    #[pattern = "'sqrt_y_adj_many"]
    #[instruction]
    pub struct SqrtYAdjMany;

    #[pattern = "'reset_many"]
    #[instruction]
    pub struct ResetMany;

    #[pattern = "'measure_many"]
    #[instruction]
    pub struct MeasureMany;

    #[pattern = "'cnot_many"]
    #[instruction]
    pub struct CNOTMany;

    #[pattern = "'cz_many"]
    #[instruction]
    pub struct CZMany;

    #[pattern = "'rx_many"]
    #[instruction]
    pub struct RXMany;

    #[pattern = "'ry_many"]
    #[instruction]
    pub struct RYMany;

    #[pattern = "'rz_many"]
    #[instruction]
    pub struct RZMany;

    #[pattern = "'depolarize_many"]
    #[instruction]
    pub struct DepolarizeMany;

    #[pattern = "'loss_many"]
    #[instruction]
    pub struct LossMany;

    #[pattern = "'rxx_many"]
    #[instruction]
    pub struct RXXMany;

    #[pattern = "'ryy_many"]
    #[instruction]
    pub struct RYYMany;

    #[pattern = "'rzz_many"]
    #[instruction]
    pub struct RZZMany;

    #[pattern = "'depolarize2_many"]
    #[instruction]
    pub struct Depolarize2Many;

    #[pattern = "'u3_many"]
    #[instruction]
    pub struct U3Many;

    #[pattern = "'paulierror_many"]
    #[instruction]
    pub struct PauliErrorMany;

    #[pattern = "'correlatedloss_many"]
    #[instruction]
    pub struct CorrelatedLossMany;

    #[pattern = "'twoqubitpaulierror_many"]
    #[instruction]
    pub struct TwoQubitPauliErrorMany;
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

    #[test]
    fn parses_batched_tokens() {
        assert_eq!(
            circuit::XMany::parser().parse("x_many").into_result(),
            Ok(circuit::XMany)
        );
        assert_eq!(
            circuit::YMany::parser().parse("y_many").into_result(),
            Ok(circuit::YMany)
        );
        assert_eq!(
            circuit::ZMany::parser().parse("z_many").into_result(),
            Ok(circuit::ZMany)
        );
        assert_eq!(
            circuit::HMany::parser().parse("h_many").into_result(),
            Ok(circuit::HMany)
        );
        assert_eq!(
            circuit::SMany::parser().parse("s_many").into_result(),
            Ok(circuit::SMany)
        );
        assert_eq!(
            circuit::SAdjMany::parser()
                .parse("s_adj_many")
                .into_result(),
            Ok(circuit::SAdjMany)
        );
        assert_eq!(
            circuit::TMany::parser().parse("t_many").into_result(),
            Ok(circuit::TMany)
        );
        assert_eq!(
            circuit::TAdjMany::parser()
                .parse("t_adj_many")
                .into_result(),
            Ok(circuit::TAdjMany)
        );
        assert_eq!(
            circuit::SqrtXMany::parser()
                .parse("sqrt_x_many")
                .into_result(),
            Ok(circuit::SqrtXMany)
        );
        assert_eq!(
            circuit::SqrtXAdjMany::parser()
                .parse("sqrt_x_adj_many")
                .into_result(),
            Ok(circuit::SqrtXAdjMany)
        );
        assert_eq!(
            circuit::SqrtYMany::parser()
                .parse("sqrt_y_many")
                .into_result(),
            Ok(circuit::SqrtYMany)
        );
        assert_eq!(
            circuit::SqrtYAdjMany::parser()
                .parse("sqrt_y_adj_many")
                .into_result(),
            Ok(circuit::SqrtYAdjMany)
        );
        assert_eq!(
            circuit::ResetMany::parser()
                .parse("reset_many")
                .into_result(),
            Ok(circuit::ResetMany)
        );
        assert_eq!(
            circuit::MeasureMany::parser()
                .parse("measure_many")
                .into_result(),
            Ok(circuit::MeasureMany)
        );
        assert_eq!(
            circuit::CNOTMany::parser().parse("cnot_many").into_result(),
            Ok(circuit::CNOTMany)
        );
        assert_eq!(
            circuit::CZMany::parser().parse("cz_many").into_result(),
            Ok(circuit::CZMany)
        );
        assert_eq!(
            circuit::RXMany::parser().parse("rx_many").into_result(),
            Ok(circuit::RXMany)
        );
        assert_eq!(
            circuit::RYMany::parser().parse("ry_many").into_result(),
            Ok(circuit::RYMany)
        );
        assert_eq!(
            circuit::RZMany::parser().parse("rz_many").into_result(),
            Ok(circuit::RZMany)
        );
        assert_eq!(
            circuit::DepolarizeMany::parser()
                .parse("depolarize_many")
                .into_result(),
            Ok(circuit::DepolarizeMany)
        );
        assert_eq!(
            circuit::LossMany::parser().parse("loss_many").into_result(),
            Ok(circuit::LossMany)
        );
        assert_eq!(
            circuit::RXXMany::parser().parse("rxx_many").into_result(),
            Ok(circuit::RXXMany)
        );
        assert_eq!(
            circuit::RYYMany::parser().parse("ryy_many").into_result(),
            Ok(circuit::RYYMany)
        );
        assert_eq!(
            circuit::RZZMany::parser().parse("rzz_many").into_result(),
            Ok(circuit::RZZMany)
        );
        assert_eq!(
            circuit::Depolarize2Many::parser()
                .parse("depolarize2_many")
                .into_result(),
            Ok(circuit::Depolarize2Many)
        );
        assert_eq!(
            circuit::U3Many::parser().parse("u3_many").into_result(),
            Ok(circuit::U3Many)
        );
        assert_eq!(
            circuit::PauliErrorMany::parser()
                .parse("paulierror_many")
                .into_result(),
            Ok(circuit::PauliErrorMany)
        );
        assert_eq!(
            circuit::CorrelatedLossMany::parser()
                .parse("correlatedloss_many")
                .into_result(),
            Ok(circuit::CorrelatedLossMany)
        );
        assert_eq!(
            circuit::TwoQubitPauliErrorMany::parser()
                .parse("twoqubitpaulierror_many")
                .into_result(),
            Ok(circuit::TwoQubitPauliErrorMany)
        );
    }
}
