// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use eyre::Result;
use smallvec::SmallVec;
use vihaco::{Effects, observe};

/// Measurement results are represented as an integer enum:
/// 0: state |0>
/// 1: state |1>
/// 2: qubit has been lost prior to measurement
/// In bytecode, this is represented as a u32 integer, which is simpler than
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

/// Per-step trace value emitted by the `Trace` instruction on the PauliSum and
/// LossyPauliSum backends. The plan's Decision 5 keeps trace and measurement
/// records as two parallel streams; this effect feeds the trace stream.
#[derive(Debug, Clone, PartialEq)]
pub struct TraceEffect {
    pub value: f64,
}

#[derive(Debug, Default)]
pub struct TraceObserver {
    pub record: Vec<f64>,
}

#[observe(TraceEffect)]
impl TraceObserver {
    fn observe_trace_effect(&mut self, effect: &TraceEffect) -> Result<Effects<()>> {
        self.record.push(effect.value);
        Ok(Effects::none())
    }
}

/// Union of the two effect types a circuit instruction can produce. The plan's
/// structural note (Task 6) calls for broadening the `#[component(..., effect =
/// ...)]` annotation on `Circuit` and the executors to a union so a single
/// `Trace` (or `Measure`) instruction can fan out to the right observer.
#[derive(Debug, Clone, PartialEq)]
pub enum CircuitOutcomeEffect {
    Measurement(MeasurementEffect),
    Trace(TraceEffect),
}

impl From<MeasurementEffect> for CircuitOutcomeEffect {
    fn from(value: MeasurementEffect) -> Self {
        Self::Measurement(value)
    }
}

impl From<TraceEffect> for CircuitOutcomeEffect {
    fn from(value: TraceEffect) -> Self {
        Self::Trace(value)
    }
}
