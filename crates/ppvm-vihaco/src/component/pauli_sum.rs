// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use eyre::{Result, eyre};
use vihaco::Effects;
use vihaco_circuit_isa::{CircuitInstruction, CircuitMessage};

use super::MAX_QUBITS;
use super::backend::{
    Backend, CorrelatedLossChannel as _, LossChannel as _, Lossy64, Lossy128, Lossy256, Lossy512,
    Lossy1024, Lossy2048, Sum64, Sum128, Sum256, Sum512, Sum1024, Sum2048, draw, executor,
};
use super::dispatch::{batch_draw, batch_pairs_draw, dispatch_pauli_sum};
use crate::device_info::PPVMDeviceInfo;
use crate::measurements::CircuitOutcomeEffect;

pub struct PauliSumExecutor<S> {
    pub state: S,
    initial: S,
    /// See `§ Where the randomness lives` in [`super::backend`].
    ///
    /// The sum channels scale coefficients analytically and never draw; this
    /// stream exists so the executor, not the caller, answers the injected-RNG
    /// gate surface.
    #[cfg(feature = "traits-2")]
    rng: rand::rngs::SmallRng,
}

pub struct LossyPauliSumExecutor<S> {
    pub state: S,
    initial: S,
    /// See `§ Where the randomness lives` in [`super::backend`].
    ///
    /// The sum channels scale coefficients analytically and never draw; this
    /// stream exists so the executor, not the caller, answers the injected-RNG
    /// gate surface.
    #[cfg(feature = "traits-2")]
    rng: rand::rngs::SmallRng,
}

pub enum PauliSumCircuit {
    Bits64(PauliSumExecutor<Sum64>),
    Bits128(PauliSumExecutor<Sum128>),
    Bits256(PauliSumExecutor<Sum256>),
    Bits512(PauliSumExecutor<Sum512>),
    Bits1024(PauliSumExecutor<Sum1024>),
    Bits2048(PauliSumExecutor<Sum2048>),
}

pub enum LossyPauliSumCircuit {
    Bits64(LossyPauliSumExecutor<Lossy64>),
    Bits128(LossyPauliSumExecutor<Lossy128>),
    Bits256(LossyPauliSumExecutor<Lossy256>),
    Bits512(LossyPauliSumExecutor<Lossy512>),
    Bits1024(LossyPauliSumExecutor<Lossy1024>),
    Bits2048(LossyPauliSumExecutor<Lossy2048>),
}

macro_rules! construct {
    ($info:expr, $terms:expr, $wrapper:ident, $variant:ident, $method:ident) => {{
        let state = Backend::$method($info, $terms);
        let initial = state.clone();
        Ok(Self::$variant(executor!(
            $wrapper {
                state: state,
                initial: initial
            },
            None
        )))
    }};
}

impl PauliSumCircuit {
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        match info.n_qubits {
            0..=64 => construct!(info, terms, PauliSumExecutor, Bits64, sum64),
            65..=128 => construct!(info, terms, PauliSumExecutor, Bits128, sum128),
            129..=256 => construct!(info, terms, PauliSumExecutor, Bits256, sum256),
            257..=512 => construct!(info, terms, PauliSumExecutor, Bits512, sum512),
            513..=1024 => construct!(info, terms, PauliSumExecutor, Bits1024, sum1024),
            1025..=MAX_QUBITS => {
                construct!(info, terms, PauliSumExecutor, Bits2048, sum2048)
            }
            n => Err(eyre!("cannot simulate {n} qubits: maximum is {MAX_QUBITS}")),
        }
    }

    pub(super) fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        if matches!(
            inst,
            CircuitInstruction::Loss | CircuitInstruction::CorrelatedLoss
        ) {
            return Err(eyre!(
                "{inst} is not supported on the PauliSum backend; use the LossyPauliSum backend instead"
            ));
        }
        match self {
            Self::Bits64(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
            Self::Bits128(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
            Self::Bits256(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
            Self::Bits512(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
            Self::Bits1024(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
            Self::Bits2048(ex) => dispatch_pauli_sum!(ex, state, inst, msg, "PauliSum"),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => Backend::render(&ex.state),
            Self::Bits128(ex) => Backend::render(&ex.state),
            Self::Bits256(ex) => Backend::render(&ex.state),
            Self::Bits512(ex) => Backend::render(&ex.state),
            Self::Bits1024(ex) => Backend::render(&ex.state),
            Self::Bits2048(ex) => Backend::render(&ex.state),
        }
    }
}

impl LossyPauliSumCircuit {
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        match info.n_qubits {
            0..=64 => construct!(info, terms, LossyPauliSumExecutor, Bits64, lossy64),
            65..=128 => construct!(info, terms, LossyPauliSumExecutor, Bits128, lossy128),
            129..=256 => construct!(info, terms, LossyPauliSumExecutor, Bits256, lossy256),
            257..=512 => construct!(info, terms, LossyPauliSumExecutor, Bits512, lossy512),
            513..=1024 => {
                construct!(info, terms, LossyPauliSumExecutor, Bits1024, lossy1024)
            }
            1025..=MAX_QUBITS => {
                construct!(info, terms, LossyPauliSumExecutor, Bits2048, lossy2048)
            }
            n => Err(eyre!("cannot simulate {n} qubits: maximum is {MAX_QUBITS}")),
        }
    }

    pub(super) fn execute_instruction(
        &mut self,
        inst: &CircuitInstruction,
        msg: &CircuitMessage,
    ) -> Result<Effects<CircuitOutcomeEffect>> {
        use vihaco_circuit_isa::{CircuitInstruction::*, CircuitMessage::*};
        macro_rules! dispatch_lossy {
            ($ex:expr) => {{
                match (inst, msg) {
                    (Loss, &QubitAndFloat(q, p)) => draw!($ex, state, loss_channel(q, p)),
                    (CorrelatedLoss, TwoQubitAndFloatArr3(a, b, p)) => {
                        draw!($ex, state, correlated_loss_channel(*a, *b, *p))
                    }
                    (Loss, QubitBatchAndFloat(qs, p)) => {
                        batch_draw!($ex, state, loss_channel, qs, *p)
                    }
                    (CorrelatedLoss, TwoQubitBatchAndFloatArr3(ps, p)) => {
                        batch_pairs_draw!($ex, state, correlated_loss_channel, ps, *p)
                    }
                    _ => return dispatch_pauli_sum!($ex, state, inst, msg, "LossyPauliSum"),
                }
                Ok(Effects::None)
            }};
        }
        match self {
            Self::Bits64(ex) => dispatch_lossy!(ex),
            Self::Bits128(ex) => dispatch_lossy!(ex),
            Self::Bits256(ex) => dispatch_lossy!(ex),
            Self::Bits512(ex) => dispatch_lossy!(ex),
            Self::Bits1024(ex) => dispatch_lossy!(ex),
            Self::Bits2048(ex) => dispatch_lossy!(ex),
        }
    }

    pub fn state_string(&self) -> String {
        match self {
            Self::Bits64(ex) => Backend::render(&ex.state),
            Self::Bits128(ex) => Backend::render(&ex.state),
            Self::Bits256(ex) => Backend::render(&ex.state),
            Self::Bits512(ex) => Backend::render(&ex.state),
            Self::Bits1024(ex) => Backend::render(&ex.state),
            Self::Bits2048(ex) => Backend::render(&ex.state),
        }
    }
}

macro_rules! reset_circuit {
    ($self:expr) => {
        match $self {
            Self::Bits64(ex) => ex.state = ex.initial.clone(),
            Self::Bits128(ex) => ex.state = ex.initial.clone(),
            Self::Bits256(ex) => ex.state = ex.initial.clone(),
            Self::Bits512(ex) => ex.state = ex.initial.clone(),
            Self::Bits1024(ex) => ex.state = ex.initial.clone(),
            Self::Bits2048(ex) => ex.state = ex.initial.clone(),
        }
    };
}

impl vihaco::Reset for PauliSumCircuit {
    fn reset(&mut self) {
        reset_circuit!(self);
    }
}

impl vihaco::Reset for LossyPauliSumCircuit {
    fn reset(&mut self) {
        reset_circuit!(self);
    }
}
