use eyre::Result;
use vihaco::{Effects, observe};

pub type MeasurementResult = Vec<Option<bool>>;

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MeasurementEffect {
    pub measurement_results: MeasurementResult,
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
