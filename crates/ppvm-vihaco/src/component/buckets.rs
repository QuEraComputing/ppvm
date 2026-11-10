// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use super::CircuitExecutor;
use super::backend::{
    LossyPauliSumConfig, MAX_QUBITS, PauliSumBackend, PauliSumConfig, paulisum_strategy,
};
use super::ops::CircuitOp;
use crate::device_info::PPVMDeviceInfo;
use bnum::types::{U256, U512, U1024, U2048};
use eyre::{Result, eyre};
use num::complex::Complex64;
use ppvm_pauli_sum::config::fx64hash::Byte8F64;
use ppvm_tableau::prelude::*;

// TODO: enum dispatch can go once we only set the index type in tableau #204
pub(crate) enum TableauBuckets {
    Bits64(CircuitExecutor<GeneralizedTableau<Byte8F64<1>, usize, Vec<(Complex64, usize)>>>),
    Bits128(CircuitExecutor<GeneralizedTableau<Byte8F64<2>, u128, Vec<(Complex64, u128)>>>),
    Bits256(CircuitExecutor<GeneralizedTableau<Byte8F64<4>, U256, Vec<(Complex64, U256)>>>),
    Bits512(CircuitExecutor<GeneralizedTableau<Byte8F64<8>, U512, Vec<(Complex64, U512)>>>),
    Bits1024(CircuitExecutor<GeneralizedTableau<Byte8F64<16>, U1024, Vec<(Complex64, U1024)>>>),
    Bits2048(CircuitExecutor<GeneralizedTableau<Byte8F64<32>, U2048, Vec<(Complex64, U2048)>>>),
}

impl TableauBuckets {
    pub(crate) fn state_string(&self) -> String {
        match self {
            Self::Bits64(c) => c.backend.to_string(),
            Self::Bits128(c) => c.backend.to_string(),
            Self::Bits256(c) => c.backend.to_string(),
            Self::Bits512(c) => c.backend.to_string(),
            Self::Bits1024(c) => c.backend.to_string(),
            Self::Bits2048(c) => c.backend.to_string(),
        }
    }

    pub(crate) fn apply<I: CircuitOp>(
        &mut self,
        op: &I,
        msg: I::Message,
    ) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    pub(crate) fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    pub fn new(n_qubits: usize, coefficient_threshold: f64) -> Result<Self> {
        if n_qubits <= 64 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits64(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 128 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits128(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 256 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits256(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 512 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits512(CircuitExecutor { backend: tab }))
        } else if n_qubits <= 1024 {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits1024(CircuitExecutor { backend: tab }))
        } else if n_qubits <= MAX_QUBITS {
            let tab = GeneralizedTableau::new(n_qubits, coefficient_threshold);
            Ok(Self::Bits2048(CircuitExecutor { backend: tab }))
        } else {
            Err(eyre!(
                "cannot simulate {n_qubits} qubits: maximum is {MAX_QUBITS}"
            ))
        }
    }

    /// Same as [`TableauCircuit::new`], but seed the RNG deterministically so a
    /// shot is reproducible.
    pub fn new_with_seed(n_qubits: usize, coefficient_threshold: f64, seed: u64) -> Result<Self> {
        macro_rules! seeded {
            ($variant:ident) => {{
                let tab = GeneralizedTableau::new_with_seed(n_qubits, coefficient_threshold, seed);
                Ok(Self::$variant(CircuitExecutor { backend: tab }))
            }};
        }
        if n_qubits <= 64 {
            seeded!(Bits64)
        } else if n_qubits <= 128 {
            seeded!(Bits128)
        } else if n_qubits <= 256 {
            seeded!(Bits256)
        } else if n_qubits <= 512 {
            seeded!(Bits512)
        } else if n_qubits <= 1024 {
            seeded!(Bits1024)
        } else if n_qubits <= MAX_QUBITS {
            seeded!(Bits2048)
        } else {
            Err(eyre!(
                "cannot simulate {n_qubits} qubits: maximum is {MAX_QUBITS}"
            ))
        }
    }
}

pub(crate) enum PauliSumBuckets {
    Bits64(CircuitExecutor<PauliSumBackend<PauliSumConfig<8>>>),
    Bits128(CircuitExecutor<PauliSumBackend<PauliSumConfig<16>>>),
    Bits256(CircuitExecutor<PauliSumBackend<PauliSumConfig<32>>>),
    Bits512(CircuitExecutor<PauliSumBackend<PauliSumConfig<64>>>),
    Bits1024(CircuitExecutor<PauliSumBackend<PauliSumConfig<128>>>),
    Bits2048(CircuitExecutor<PauliSumBackend<PauliSumConfig<256>>>),
}

impl PauliSumBuckets {
    pub(crate) fn state_string(&self) -> String {
        match self {
            Self::Bits64(c) => c.backend.state.to_string(),
            Self::Bits128(c) => c.backend.state.to_string(),
            Self::Bits256(c) => c.backend.state.to_string(),
            Self::Bits512(c) => c.backend.state.to_string(),
            Self::Bits1024(c) => c.backend.state.to_string(),
            Self::Bits2048(c) => c.backend.state.to_string(),
        }
    }

    pub(crate) fn apply<I: CircuitOp>(
        &mut self,
        op: &I,
        msg: I::Message,
    ) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    pub(crate) fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    /// Build a PauliSum-backed circuit, seeding the state with every term:
    /// `for (word, coef) in terms { state += (word, coef); }`. Words must
    /// already be validated against `info.n_qubits` by the caller.
    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<PauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                let initial_state = state.clone();
                let psum = PauliSumBackend::<PauliSumConfig<$N>> {
                    state,
                    initial_state,
                };
                Ok(Self::$variant(CircuitExecutor { backend: psum }))
            }};
        }
        if info.n_qubits <= 64 {
            build!(Bits64, 8)
        } else if info.n_qubits <= 128 {
            build!(Bits128, 16)
        } else if info.n_qubits <= 256 {
            build!(Bits256, 32)
        } else if info.n_qubits <= 512 {
            build!(Bits512, 64)
        } else if info.n_qubits <= 1024 {
            build!(Bits1024, 128)
        } else if info.n_qubits <= MAX_QUBITS {
            build!(Bits2048, 256)
        } else {
            Err(eyre!(
                "cannot simulate {} qubits: maximum is {MAX_QUBITS}",
                info.n_qubits
            ))
        }
    }
}

pub(crate) enum LossyPauliSumBuckets {
    Bits64(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<8>>>),
    Bits128(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<16>>>),
    Bits256(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<32>>>),
    Bits512(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<64>>>),
    Bits1024(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<128>>>),
    Bits2048(CircuitExecutor<PauliSumBackend<LossyPauliSumConfig<256>>>),
}

impl LossyPauliSumBuckets {
    pub(crate) fn state_string(&self) -> String {
        match self {
            Self::Bits64(c) => c.backend.state.to_string(),
            Self::Bits128(c) => c.backend.state.to_string(),
            Self::Bits256(c) => c.backend.state.to_string(),
            Self::Bits512(c) => c.backend.state.to_string(),
            Self::Bits1024(c) => c.backend.state.to_string(),
            Self::Bits2048(c) => c.backend.state.to_string(),
        }
    }

    pub(crate) fn apply<I: CircuitOp>(
        &mut self,
        op: &I,
        msg: I::Message,
    ) -> eyre::Result<I::Effect> {
        match self {
            Self::Bits64(c) => op.apply(&mut c.backend, msg),
            Self::Bits128(c) => op.apply(&mut c.backend, msg),
            Self::Bits256(c) => op.apply(&mut c.backend, msg),
            Self::Bits512(c) => op.apply(&mut c.backend, msg),
            Self::Bits1024(c) => op.apply(&mut c.backend, msg),
            Self::Bits2048(c) => op.apply(&mut c.backend, msg),
        }
    }

    pub(crate) fn reset(&mut self) {
        match self {
            Self::Bits64(c) => vihaco::Reset::reset(c),
            Self::Bits128(c) => vihaco::Reset::reset(c),
            Self::Bits256(c) => vihaco::Reset::reset(c),
            Self::Bits512(c) => vihaco::Reset::reset(c),
            Self::Bits1024(c) => vihaco::Reset::reset(c),
            Self::Bits2048(c) => vihaco::Reset::reset(c),
        }
    }

    pub fn new(info: &PPVMDeviceInfo, terms: &[(String, f64)]) -> Result<Self> {
        macro_rules! build {
            ($variant:ident, $N:literal) => {{
                let mut state = PauliSum::<LossyPauliSumConfig<$N>>::builder()
                    .n_qubits(info.n_qubits)
                    .strategy(paulisum_strategy(info))
                    .build();
                for (word, coef) in terms {
                    state += (word.as_str(), *coef);
                }
                let initial_state = state.clone();
                let psum = PauliSumBackend::<LossyPauliSumConfig<$N>> {
                    state,
                    initial_state,
                };
                Ok(Self::$variant(CircuitExecutor { backend: psum }))
            }};
        }
        if info.n_qubits <= 64 {
            build!(Bits64, 8)
        } else if info.n_qubits <= 128 {
            build!(Bits128, 16)
        } else if info.n_qubits <= 256 {
            build!(Bits256, 32)
        } else if info.n_qubits <= 512 {
            build!(Bits512, 64)
        } else if info.n_qubits <= 1024 {
            build!(Bits1024, 128)
        } else if info.n_qubits <= MAX_QUBITS {
            build!(Bits2048, 256)
        } else {
            Err(eyre!(
                "cannot simulate {} qubits: maximum is {MAX_QUBITS}",
                info.n_qubits
            ))
        }
    }
}

pub(crate) enum CircuitKind {
    Tableau(TableauBuckets),
    PauliSum(PauliSumBuckets),
    LossyPauliSum(LossyPauliSumBuckets),
}

impl CircuitKind {
    pub(crate) fn state_string(&self) -> String {
        match self {
            Self::Tableau(b) => b.state_string(),
            Self::PauliSum(b) => b.state_string(),
            Self::LossyPauliSum(b) => b.state_string(),
        }
    }

    pub(crate) fn reset(&mut self) {
        match self {
            Self::Tableau(b) => b.reset(),
            Self::PauliSum(b) => b.reset(),
            Self::LossyPauliSum(b) => b.reset(),
        }
    }
}
