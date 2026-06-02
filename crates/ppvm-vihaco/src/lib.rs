pub mod bytecode;
pub mod component;
pub mod composite;
pub mod instruction;
pub mod measurements;
pub mod message;
mod syntax;

use chumsky::Parser;
use vihaco::syntax::{ParsedModule, Resolve};
use vihaco::{Type, Value, module::Module};
use vihaco_parser_core::Parse;

use crate::composite::{PPVM, PPVMDeviceInfo, PPVMInstruction};
use crate::syntax::{PPVMHeader, PPVMResolver};

pub fn run_file(path: &str) -> eyre::Result<PPVM> {
    let mut machine = PPVM::default();
    machine.run_file(path)?;
    Ok(machine)
}

pub fn run_program(program: &str) -> eyre::Result<PPVM> {
    let mut machine = PPVM::default();
    machine.run_program(program)?;
    Ok(machine)
}

/// Parse `.sst` source into the unresolved AST.
pub fn parse_program(source: &str) -> eyre::Result<ParsedModule<PPVMInstruction, PPVMHeader>> {
    ParsedModule::<PPVMInstruction, PPVMHeader>::parser()
        .parse(source)
        .into_result()
        .map_err(|errs| eyre::eyre!("parsing failed: {errs:?}"))
}

pub fn compile_program(
    source: &str,
) -> eyre::Result<Module<PPVMInstruction, Value, Type, PPVMDeviceInfo>> {
    PPVMResolver::new().resolve_module(parse_program(source)?)
}

/// Dump `.sst` source to a `.ssb` bytecode file.
pub fn dump_program(program: &str, output_path: &str) -> eyre::Result<()> {
    let bytes = bytecode::compile_to_bytes(program)?;
    std::fs::write(output_path, bytes)?;
    Ok(())
}

/// Read a `.sst` file and dump it to a `.ssb` bytecode file.
pub fn dump_file(input_path: &str, output_path: &str) -> eyre::Result<()> {
    let program = std::fs::read_to_string(input_path)?;
    dump_program(&program, output_path)
}

pub mod prelude {
    pub use crate::component::Circuit;
    pub use crate::composite::PPVM;
    pub use crate::syntax::{PPVMHeader, PPVMResolver};
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dump_program_writes_loadable_bytecode() {
        let src = "device circuit.n_qubits 1;\n\
                   fn @main() { const.u64 0\n gate measure\n ret }\n";
        let path = std::env::temp_dir().join("ppvm_dump_program_test.ssb");
        dump_program(src, path.to_str().unwrap()).unwrap();

        let mut machine = PPVM::default();
        machine.load_bytecode_file(path.to_str().unwrap()).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn dump_file_reads_sst_and_writes_bytecode() {
        let out = std::env::temp_dir().join("ppvm_dump_file_test.ssb");
        dump_file("tests/function_call.sst", out.to_str().unwrap()).unwrap();

        let mut machine = PPVM::default();
        machine.load_bytecode_file(out.to_str().unwrap()).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 1);
        let _ = std::fs::remove_file(&out);
    }
}
