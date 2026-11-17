// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

pub mod bytecode;
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

/// Compile `.sst` source and write it to a `.ssb` bytecode file.
pub fn dump_program(program: &str, output_path: &str) -> eyre::Result<()> {
    std::fs::write(output_path, bytecode::compile_to_bytes(program)?)?;
    Ok(())
}

/// Read a `.sst` file and dump it to a `.ssb` bytecode file.
pub fn dump_file(input_path: &str, output_path: &str) -> eyre::Result<()> {
    dump_program(&std::fs::read_to_string(input_path)?, output_path)
}

pub mod prelude {
    pub use crate::component::Circuit;
    pub use crate::composite::PPVM;
}
