use eyre::Result;
use smallvec::SmallVec;
use vihaco::{Effects, observe};

/// Measurement results are represent as an integer enum
/// 0: state |0>
/// 1: state |1>
/// 2: qubit has been lost prior to measurement
/// In byte-code, this is represented as a u32 integer, which is simpler than
/// e.g. two boolean values and matches semantics elsewhere
#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum MeasurementOutcome {
    Zero = 0,
    One = 1,
    Lost = 2,
}

pub type MeasurementResult = SmallVec<[MeasurementOutcome; 8]>;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MeasurementEffect {
    pub measurement_results: MeasurementResult,
}

impl From<Option<bool>> for MeasurementOutcome {
    fn from(m: Option<bool>) -> Self {
        match m {
            Some(false) => Self::Zero,
            Some(true) => Self::One,
            None => Self::Lost,
        }
    }
}

#[derive(Debug, Default)]
pub struct MeasurementObserver {
    pub record: Vec<MeasurementResult>,
}

#[observe(MeasurementEffect)]
impl MeasurementObserver {
    fn observe_measurement_effect(&mut self, effect: &MeasurementEffect) -> Result<Effects<()>> {
        self.record.push(effect.measurement_results.clone());
        Ok(Effects::none())
    }
}
