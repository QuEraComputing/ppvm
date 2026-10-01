// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! GPU Pauli propagation for ppvm on NVIDIA cuPauliProp.
//!
//! - [`CudaPauliSum`] (feature `cuda`): a device-resident Pauli sum with
//!   ppvm-named gate methods, propagated in the Heisenberg picture.
//! - [`CpuPauliSum`]: the ppvm `PauliSum` configured exactly like the Python
//!   `ppvm.PauliSum`, used as the CPU baseline.
//! - [`circuits`]: the TFIM and Heisenberg Trotter circuits from
//!   `ppvm-benchmarks/trotter-benchmarks`, generic over [`Propagator`].
//!
//! The `trotter_cpu` / `trotter_gpu` examples run those circuits.

pub mod circuits;
pub mod cpu;
#[cfg(feature = "cuda")]
pub mod cuda;
pub mod encode;

pub use circuits::Propagator;
pub use cpu::CpuPauliSum;
#[cfg(feature = "cuda")]
pub use cuda::{CudaOptions, CudaPauliSum};
