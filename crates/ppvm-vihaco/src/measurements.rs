use eyre::Result;
use vihaco::{Effects, observe};

#[repr(u8)]
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum MeasurementOutcome {
    Zero = 0,
    One = 1,
    Lost = 2,
}

pub type MeasurementResult = Vec<MeasurementOutcome>;

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
