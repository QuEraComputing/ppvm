// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Turning PPVM program source into a [`PPVMModule`] ready to execute.

use crate::composite::{PPVM, PPVMModule};
use crate::header::{HeaderDeclaration, collect_header};
use vihaco::loader::MetadataLoweringContext;
use vihaco::syntax::{HeaderDeclarations, ParsedModule, lower_functions_with_metadata};
use vihaco::{HasInstructionSet, NoContext, SstFile, Type};

/// Name of the SST section holding the PPVM program.
pub const PPVM_SECTION: &str = "ppvm";

type SourceInstruction = <PPVM as HasInstructionSet>::Syntax;

/// Parse and lower SST source into a module whose `extra` is the device info.
pub fn compile_program(source: &str) -> eyre::Result<PPVMModule> {
    let file = SstFile::<NoContext>::from_text(source)?;
    let section = file
        .root()
        .child(PPVM_SECTION)
        .ok_or_else(|| eyre::eyre!("missing `{PPVM_SECTION}` section"))?;

    let ParsedModule { header, functions } = ParsedModule::<
        SourceInstruction,
        Type,
        HeaderDeclarations<HeaderDeclaration>,
    >::parse_section(section)?;
    let (metadata, info) = collect_header(header)?;

    let mut context = MetadataLoweringContext::new(&metadata)?;
    let code = lower_functions_with_metadata(functions, &metadata, &mut context)?;
    metadata.map_extra(|_| info).into_module(code)
}

/// Read a file and compile it with [`compile_program`].
pub fn load_module_file(path: &str) -> eyre::Result<PPVMModule> {
    compile_program(&std::fs::read_to_string(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measurements::MeasurementOutcome;

    /// X on |0> then measure: the outcome is deterministically 1.
    const FLIP_AND_MEASURE: &str = r#"sst v1
.section(root):
    .section(ppvm):
    .header(ppvm):
        circuit.n_qubits 1;
        metadata {
            strings: ["main", "flip.sst"],
            functions: [
                {
                    name: 0,
                    signature: { params: [], ret: [] },
                    local_count: 0,
                    start_address: 0,
                    end_address: 5,
                    file: 1,
                },
            ],
            labels: [],
            constants: [],
            main_function: 0,
            source_symbols: [],
            extra: {},
        };
    .header(ppvm).
    .text(ppvm):
        fn @main() {
            execution::program::const_u64 0
            circuit::circuit::x
            execution::program::const_u64 0
            circuit::circuit::measure
            execution::program::ret 0
        }
    .text(ppvm).
    .section(ppvm).
.section(root).
"#;

    #[test]
    fn compiled_program_runs_on_the_declared_device() -> eyre::Result<()> {
        let module = compile_program(FLIP_AND_MEASURE)?;
        assert_eq!(module.extra.n_qubits, 1);

        let mut machine = PPVM::new(&module.extra, Some(0))?;
        vihaco::execute(&module, &mut machine, &[])?;
        assert_eq!(machine.measurement_record().len(), 1);
        assert_eq!(
            machine.measurement_record()[0].as_slice(),
            [MeasurementOutcome::One]
        );
        Ok(())
    }
}
