// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

// pub mod bytecode;
pub mod component;
pub mod composite;
pub mod device_info;
pub mod header;
pub mod loader;
pub mod measurements;
pub mod observable;
pub mod shots;

/// Re-exported so consumers (e.g. the TUI) can name gates for
/// [`composite::PPVM::send`] without depending on the ISA crate directly.
pub use vihaco_circuit_isa::circuit;

pub use composite::PPVMModule;
pub use loader::{compile_program, load_module_file, parse_program};

pub mod prelude {
    pub use crate::component::Circuit;
    pub use crate::composite::PPVM;
}
