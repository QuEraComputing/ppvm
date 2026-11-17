// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Binary bytecode (`.ssb`) round-trip for PPVM modules.
//!
//! Serializes a compiled [`PPVMModule`] to a little-endian container — header
//! (magic + version + device info) → strings → function/label metadata → code —
//! and reads it back. Each instruction is a `u16` opcode followed by its
//! payload fields, as listed in the opcode table below.

use std::io::{Read, Write};
use vihaco::instruction::{FromBytes, WriteBytes};
use vihaco::module::{FunctionInfo, LabelInfo, Parameter, Signature};
use vihaco::{HasInstructionSet, Type};
use vihaco_circuit_isa::circuit;
use vihaco_cpu::cpu;
use vihaco_stdlib::debugger::dialect::debug;
use vihaco_stdlib::program::program;

use crate::PPVMModule;
use crate::composite::{BackendKind, PPVM, PPVM_MAGIC, PPVMDeviceInfo};

type Instruction = <PPVM as HasInstructionSet>::Runtime;

/// Generate the opcode enum and the instruction encoder/decoder from one table.
///
/// An entry's opcode is its position in the table, so the table is append-only:
/// reordering or removing an entry changes the opcodes of every later one.
macro_rules! opcode_table {
    ($( $variant:ident = $($seg:ident)::+ $(( $($field:ident: $ty:ty),* ))? ),* $(,)?) => {
        #[allow(non_camel_case_types, clippy::enum_variant_names)]
        #[repr(u16)]
        enum Opcode {
            $( $variant, )*
        }

        fn write_instruction<W: Write>(inst: &Instruction, w: &mut W) -> eyre::Result<()> {
            match inst {
                $( Instruction::$variant($($seg)::+ $(( $($field),* ))?) => {
                    w.write_all(&(Opcode::$variant as u16).to_le_bytes())?;
                    $( $( $field.write_bytes(w)?; )* )?
                } )*
            }
            Ok(())
        }

        fn read_instruction<R: Read>(r: &mut R) -> eyre::Result<Instruction> {
            let opcode = read_u16(r)?;
            $( if opcode == Opcode::$variant as u16 {
                return Ok(Instruction::$variant($($seg)::+ $(( $( <$ty>::from_bytes(r)? ),* ))?));
            } )*
            Err(eyre::eyre!("unknown opcode {opcode}"))
        }

        /// One instance of every instruction, with distinct payload values.
        #[cfg(test)]
        fn every_instruction() -> Vec<Instruction> {
            let mut next = 0u8;
            let mut value = || {
                next = next.wrapping_add(1);
                next
            };
            vec![$( Instruction::$variant($($seg)::+ $(( $( <$ty>::from({ let _ = stringify!($field); value() }) ),* ))?), )*]
        }
    };
}

// Append new instructions at the end of the table; never reorder or remove.
opcode_table! {
    // execution::program
    V_execution__program__Branch = program::Branch(a: u32),
    V_execution__program__ConditionalBranch = program::ConditionalBranch(a: u32, b: u32),
    V_execution__program__Return = program::Return(a: u32),
    V_execution__program__IndirectCall = program::IndirectCall,
    V_execution__program__Call = program::Call(a: u32, b: u32),
    V_execution__program__Halt = program::Halt,
    V_execution__program__LoadI32 = program::LoadI32(a: u32),
    V_execution__program__LoadI64 = program::LoadI64(a: u32),
    V_execution__program__LoadU32 = program::LoadU32(a: u32),
    V_execution__program__LoadU64 = program::LoadU64(a: u32),
    V_execution__program__LoadF32 = program::LoadF32(a: u32),
    V_execution__program__LoadF64 = program::LoadF64(a: u32),
    V_execution__program__LoadBool = program::LoadBool(a: u32),
    V_execution__program__StoreI32 = program::StoreI32(a: u32),
    V_execution__program__StoreI64 = program::StoreI64(a: u32),
    V_execution__program__StoreU32 = program::StoreU32(a: u32),
    V_execution__program__StoreU64 = program::StoreU64(a: u32),
    V_execution__program__StoreF32 = program::StoreF32(a: u32),
    V_execution__program__StoreF64 = program::StoreF64(a: u32),
    V_execution__program__StoreBool = program::StoreBool(a: u32),
    V_execution__program__Dup = program::Dup,
    V_execution__program__ConstI32 = program::ConstI32(a: u64),
    V_execution__program__ConstI64 = program::ConstI64(a: u64),
    V_execution__program__ConstU32 = program::ConstU32(a: u64),
    V_execution__program__ConstU64 = program::ConstU64(a: u64),
    V_execution__program__ConstF32 = program::ConstF32(a: u64),
    V_execution__program__ConstF64 = program::ConstF64(a: u64),
    V_execution__program__ConstBool = program::ConstBool(a: u64),
    V_execution__program__ConstString = program::ConstString(a: u64),
    V_execution__program__ConstFunctionRef = program::ConstFunctionRef(a: u64),
    V_execution__program__ConstHeapRef = program::ConstHeapRef(a: u64),
    // cpu::cpu
    V_cpu__cpu__HeapAlloc = cpu::HeapAlloc(a: u32),
    V_cpu__cpu__GetItem = cpu::GetItem,
    V_cpu__cpu__HeapDealloc = cpu::HeapDealloc,
    V_cpu__cpu__HeapReserve = cpu::HeapReserve,
    V_cpu__cpu__HeapPush = cpu::HeapPush,
    V_cpu__cpu__HeapLen = cpu::HeapLen,
    V_cpu__cpu__AddI32 = cpu::AddI32,
    V_cpu__cpu__AddF32 = cpu::AddF32,
    V_cpu__cpu__AddI64 = cpu::AddI64,
    V_cpu__cpu__AddU32 = cpu::AddU32,
    V_cpu__cpu__AddU64 = cpu::AddU64,
    V_cpu__cpu__AddF64 = cpu::AddF64,
    V_cpu__cpu__SubI32 = cpu::SubI32,
    V_cpu__cpu__SubI64 = cpu::SubI64,
    V_cpu__cpu__SubU32 = cpu::SubU32,
    V_cpu__cpu__SubU64 = cpu::SubU64,
    V_cpu__cpu__SubF32 = cpu::SubF32,
    V_cpu__cpu__SubF64 = cpu::SubF64,
    V_cpu__cpu__MulI32 = cpu::MulI32,
    V_cpu__cpu__MulI64 = cpu::MulI64,
    V_cpu__cpu__MulU32 = cpu::MulU32,
    V_cpu__cpu__MulU64 = cpu::MulU64,
    V_cpu__cpu__MulF32 = cpu::MulF32,
    V_cpu__cpu__MulF64 = cpu::MulF64,
    V_cpu__cpu__DivI32 = cpu::DivI32,
    V_cpu__cpu__DivI64 = cpu::DivI64,
    V_cpu__cpu__DivU32 = cpu::DivU32,
    V_cpu__cpu__DivU64 = cpu::DivU64,
    V_cpu__cpu__DivF32 = cpu::DivF32,
    V_cpu__cpu__DivF64 = cpu::DivF64,
    V_cpu__cpu__RemI32 = cpu::RemI32,
    V_cpu__cpu__RemI64 = cpu::RemI64,
    V_cpu__cpu__RemU32 = cpu::RemU32,
    V_cpu__cpu__RemU64 = cpu::RemU64,
    V_cpu__cpu__RemF32 = cpu::RemF32,
    V_cpu__cpu__RemF64 = cpu::RemF64,
    V_cpu__cpu__NegI32 = cpu::NegI32,
    V_cpu__cpu__NegI64 = cpu::NegI64,
    V_cpu__cpu__NegF32 = cpu::NegF32,
    V_cpu__cpu__NegF64 = cpu::NegF64,
    V_cpu__cpu__ShlI32 = cpu::ShlI32,
    V_cpu__cpu__ShlI64 = cpu::ShlI64,
    V_cpu__cpu__ShlU32 = cpu::ShlU32,
    V_cpu__cpu__ShlU64 = cpu::ShlU64,
    V_cpu__cpu__ShrI32 = cpu::ShrI32,
    V_cpu__cpu__ShrI64 = cpu::ShrI64,
    V_cpu__cpu__ShrU32 = cpu::ShrU32,
    V_cpu__cpu__ShrU64 = cpu::ShrU64,
    V_cpu__cpu__RolI32 = cpu::RolI32,
    V_cpu__cpu__RolI64 = cpu::RolI64,
    V_cpu__cpu__RolU32 = cpu::RolU32,
    V_cpu__cpu__RolU64 = cpu::RolU64,
    V_cpu__cpu__RorI32 = cpu::RorI32,
    V_cpu__cpu__RorI64 = cpu::RorI64,
    V_cpu__cpu__RorU32 = cpu::RorU32,
    V_cpu__cpu__RorU64 = cpu::RorU64,
    V_cpu__cpu__BitAndI32 = cpu::BitAndI32,
    V_cpu__cpu__BitAndI64 = cpu::BitAndI64,
    V_cpu__cpu__BitAndU32 = cpu::BitAndU32,
    V_cpu__cpu__BitAndU64 = cpu::BitAndU64,
    V_cpu__cpu__BitOrI32 = cpu::BitOrI32,
    V_cpu__cpu__BitOrI64 = cpu::BitOrI64,
    V_cpu__cpu__BitOrU32 = cpu::BitOrU32,
    V_cpu__cpu__BitOrU64 = cpu::BitOrU64,
    V_cpu__cpu__BitXorI32 = cpu::BitXorI32,
    V_cpu__cpu__BitXorI64 = cpu::BitXorI64,
    V_cpu__cpu__BitXorU32 = cpu::BitXorU32,
    V_cpu__cpu__BitXorU64 = cpu::BitXorU64,
    V_cpu__cpu__Not = cpu::Not,
    V_cpu__cpu__And = cpu::And,
    V_cpu__cpu__Or = cpu::Or,
    V_cpu__cpu__Xor = cpu::Xor,
    V_cpu__cpu__EqI32 = cpu::EqI32,
    V_cpu__cpu__EqI64 = cpu::EqI64,
    V_cpu__cpu__EqU32 = cpu::EqU32,
    V_cpu__cpu__EqU64 = cpu::EqU64,
    V_cpu__cpu__EqF32 = cpu::EqF32,
    V_cpu__cpu__EqF64 = cpu::EqF64,
    V_cpu__cpu__NeI32 = cpu::NeI32,
    V_cpu__cpu__NeI64 = cpu::NeI64,
    V_cpu__cpu__NeU32 = cpu::NeU32,
    V_cpu__cpu__NeU64 = cpu::NeU64,
    V_cpu__cpu__NeF32 = cpu::NeF32,
    V_cpu__cpu__NeF64 = cpu::NeF64,
    V_cpu__cpu__LtI32 = cpu::LtI32,
    V_cpu__cpu__LtI64 = cpu::LtI64,
    V_cpu__cpu__LtU32 = cpu::LtU32,
    V_cpu__cpu__LtU64 = cpu::LtU64,
    V_cpu__cpu__LtF32 = cpu::LtF32,
    V_cpu__cpu__LtF64 = cpu::LtF64,
    V_cpu__cpu__GtI32 = cpu::GtI32,
    V_cpu__cpu__GtI64 = cpu::GtI64,
    V_cpu__cpu__GtU32 = cpu::GtU32,
    V_cpu__cpu__GtU64 = cpu::GtU64,
    V_cpu__cpu__GtF32 = cpu::GtF32,
    V_cpu__cpu__GtF64 = cpu::GtF64,
    V_cpu__cpu__LeI32 = cpu::LeI32,
    V_cpu__cpu__LeI64 = cpu::LeI64,
    V_cpu__cpu__LeU32 = cpu::LeU32,
    V_cpu__cpu__LeU64 = cpu::LeU64,
    V_cpu__cpu__LeF32 = cpu::LeF32,
    V_cpu__cpu__LeF64 = cpu::LeF64,
    V_cpu__cpu__GeI32 = cpu::GeI32,
    V_cpu__cpu__GeI64 = cpu::GeI64,
    V_cpu__cpu__GeU32 = cpu::GeU32,
    V_cpu__cpu__GeU64 = cpu::GeU64,
    V_cpu__cpu__GeF32 = cpu::GeF32,
    V_cpu__cpu__GeF64 = cpu::GeF64,
    // circuit::circuit
    V_circuit__circuit__XMany = circuit::XMany,
    V_circuit__circuit__YMany = circuit::YMany,
    V_circuit__circuit__ZMany = circuit::ZMany,
    V_circuit__circuit__HMany = circuit::HMany,
    V_circuit__circuit__SMany = circuit::SMany,
    V_circuit__circuit__SAdjMany = circuit::SAdjMany,
    V_circuit__circuit__TMany = circuit::TMany,
    V_circuit__circuit__TAdjMany = circuit::TAdjMany,
    V_circuit__circuit__SqrtXMany = circuit::SqrtXMany,
    V_circuit__circuit__SqrtXAdjMany = circuit::SqrtXAdjMany,
    V_circuit__circuit__SqrtYMany = circuit::SqrtYMany,
    V_circuit__circuit__SqrtYAdjMany = circuit::SqrtYAdjMany,
    V_circuit__circuit__ResetMany = circuit::ResetMany,
    V_circuit__circuit__MeasureMany = circuit::MeasureMany,
    V_circuit__circuit__CNOTMany = circuit::CNOTMany,
    V_circuit__circuit__CZMany = circuit::CZMany,
    V_circuit__circuit__RXMany = circuit::RXMany,
    V_circuit__circuit__RYMany = circuit::RYMany,
    V_circuit__circuit__RZMany = circuit::RZMany,
    V_circuit__circuit__DepolarizeMany = circuit::DepolarizeMany,
    V_circuit__circuit__LossMany = circuit::LossMany,
    V_circuit__circuit__RXXMany = circuit::RXXMany,
    V_circuit__circuit__RYYMany = circuit::RYYMany,
    V_circuit__circuit__RZZMany = circuit::RZZMany,
    V_circuit__circuit__Depolarize2Many = circuit::Depolarize2Many,
    V_circuit__circuit__U3Many = circuit::U3Many,
    V_circuit__circuit__PauliErrorMany = circuit::PauliErrorMany,
    V_circuit__circuit__CorrelatedLossMany = circuit::CorrelatedLossMany,
    V_circuit__circuit__TwoQubitPauliErrorMany = circuit::TwoQubitPauliErrorMany,
    V_circuit__circuit__TwoQubitPauliError = circuit::TwoQubitPauliError,
    V_circuit__circuit__Truncate = circuit::Truncate,
    V_circuit__circuit__Trace = circuit::Trace,
    V_circuit__circuit__X = circuit::X,
    V_circuit__circuit__Y = circuit::Y,
    V_circuit__circuit__Z = circuit::Z,
    V_circuit__circuit__H = circuit::H,
    V_circuit__circuit__SqrtXAdj = circuit::SqrtXAdj,
    V_circuit__circuit__SqrtX = circuit::SqrtX,
    V_circuit__circuit__SqrtYAdj = circuit::SqrtYAdj,
    V_circuit__circuit__SqrtY = circuit::SqrtY,
    V_circuit__circuit__SAdj = circuit::SAdj,
    V_circuit__circuit__S = circuit::S,
    V_circuit__circuit__CNOT = circuit::CNOT,
    V_circuit__circuit__CZ = circuit::CZ,
    V_circuit__circuit__TAdj = circuit::TAdj,
    V_circuit__circuit__T = circuit::T,
    V_circuit__circuit__RXX = circuit::RXX,
    V_circuit__circuit__RYY = circuit::RYY,
    V_circuit__circuit__RZZ = circuit::RZZ,
    V_circuit__circuit__RX = circuit::RX,
    V_circuit__circuit__RY = circuit::RY,
    V_circuit__circuit__RZ = circuit::RZ,
    V_circuit__circuit__U3 = circuit::U3,
    V_circuit__circuit__Measure = circuit::Measure,
    V_circuit__circuit__Reset = circuit::Reset,
    V_circuit__circuit__R = circuit::R,
    V_circuit__circuit__Loss = circuit::Loss,
    V_circuit__circuit__CorrelatedLoss = circuit::CorrelatedLoss,
    V_circuit__circuit__PauliError = circuit::PauliError,
    V_circuit__circuit__Depolarize2 = circuit::Depolarize2,
    V_circuit__circuit__Depolarize = circuit::Depolarize,
    V_circuit__circuit__Leakage = circuit::Leakage,
    // debugger::debug
    V_debugger__debug__Breakpoint = debug::Breakpoint,
}

/// Current `.ssb` format version. The reader rejects any other version.
pub const PPVM_BYTECODE_VERSION: u16 = 3;

/// Byte length of the fixed portion of the header. The actual `header_size`
/// in the stream may exceed this when the optional `observable` string is
/// populated; the reader uses `header_size` to skip to the strings section.
///
/// Field widths (bytes): magic(4) + version(2) + header_size(4) + n_qubits(4)
/// + coefficient_threshold(8) + backend(1) + max_pauli_weight_present(1)
/// + max_pauli_weight(8) + observable_present(1) = 33.
const FIXED_HEADER_SIZE: u32 = 4 + 2 + 4 + 4 + 8 + 1 + 1 + 8 + 1;

/// Serialize a compiled module to the `.ssb` byte stream.
pub fn write_module<W: Write>(module: &PPVMModule, w: &mut W) -> eyre::Result<()> {
    if !module.constants.is_empty() || !module.source_symbols.is_empty() {
        return Err(eyre::eyre!(
            "bytecode cannot represent constants or source symbols"
        ));
    }

    let info = &module.extra;
    let n_qubits = u32::try_from(info.n_qubits)
        .map_err(|_| eyre::eyre!("n_qubits {} does not fit in u32", info.n_qubits))?;

    // The header is `FIXED_HEADER_SIZE` bytes plus, when an observable is
    // present, a u32 length followed by its UTF-8 bytes.
    let observable_bytes: &[u8] = info
        .observable
        .as_ref()
        .map(String::as_bytes)
        .unwrap_or(&[]);
    let observable_present: u8 = u8::from(info.observable.is_some());
    let observable_len = u32::try_from(observable_bytes.len()).map_err(|_| {
        eyre::eyre!(
            "observable length {} does not fit in u32",
            observable_bytes.len()
        )
    })?;
    let observable_trailer: u32 = if info.observable.is_some() {
        4 + observable_len
    } else {
        0
    };
    let header_size = FIXED_HEADER_SIZE + observable_trailer;

    // Header.
    w.write_all(&PPVM_MAGIC.to_le_bytes())?;
    w.write_all(&PPVM_BYTECODE_VERSION.to_le_bytes())?;
    w.write_all(&header_size.to_le_bytes())?;
    w.write_all(&n_qubits.to_le_bytes())?;
    w.write_all(&info.coefficient_threshold.to_le_bytes())?;
    w.write_all(&[backend_to_u8(info.backend)])?;
    let (mpw_present, mpw_value) = match info.max_pauli_weight {
        Some(w) => (
            1u8,
            u64::try_from(w)
                .map_err(|_| eyre::eyre!("max_pauli_weight {} does not fit in u64", w))?,
        ),
        None => (0u8, 0u64),
    };
    w.write_all(&[mpw_present])?;
    w.write_all(&mpw_value.to_le_bytes())?;
    w.write_all(&[observable_present])?;
    if info.observable.is_some() {
        w.write_all(&observable_len.to_le_bytes())?;
        w.write_all(observable_bytes)?;
    }

    // Strings section: count, then each entry as len-prefixed UTF-8.
    let string_count =
        u32::try_from(module.strings.len()).map_err(|_| eyre::eyre!("string count exceeds u32"))?;
    w.write_all(&string_count.to_le_bytes())?;
    for s in &module.strings {
        let len = u32::try_from(s.len()).map_err(|_| eyre::eyre!("string length exceeds u32"))?;
        w.write_all(&len.to_le_bytes())?;
        w.write_all(s.as_bytes())?;
    }

    // Metadata section: function and label tables are required for calls and
    // branches to retain their resolved targets after a bytecode round trip.
    write_u32(w, u32::try_from(module.functions.len())?)?;
    for function in &module.functions {
        write_u32(w, function.name)?;
        write_u32(w, function.local_count)?;
        write_u32(w, function.start_address)?;
        write_u32(w, function.end_address)?;
        write_u32(w, function.file)?;
        write_u32(w, u32::try_from(function.signature.params.len())?)?;
        for parameter in &function.signature.params {
            write_u32(w, parameter.name)?;
            parameter.ty.write_bytes(w)?;
        }
        write_u32(w, u32::try_from(function.signature.ret.len())?)?;
        for ty in &function.signature.ret {
            ty.write_bytes(w)?;
        }
    }
    write_u32(w, u32::try_from(module.labels.len())?)?;
    for label in &module.labels {
        write_u32(w, label.address)?;
        write_u32(w, label.name)?;
    }
    write_u32(w, module.main_function.unwrap_or(u32::MAX))?;

    // Code section: count, then each instruction's fixed-width frame.
    let code_count =
        u32::try_from(module.code.len()).map_err(|_| eyre::eyre!("code length exceeds u32"))?;
    w.write_all(&code_count.to_le_bytes())?;
    for inst in &module.code {
        write_instruction(inst, w)?;
    }

    Ok(())
}

/// Reconstruct a module from a `.ssb` byte stream.
pub fn read_module<R: Read>(r: &mut R) -> eyre::Result<PPVMModule> {
    // Header.
    let magic = read_u32(r)?;
    if magic != PPVM_MAGIC {
        return Err(eyre::eyre!(
            "not a PPVM bytecode file (magic 0x{magic:08X})"
        ));
    }
    let version = read_u16(r)?;
    if version != PPVM_BYTECODE_VERSION {
        return Err(eyre::eyre!("unsupported bytecode version {version}"));
    }
    let header_size = read_u32(r)?;
    let n_qubits = read_u32(r)? as usize;
    let coefficient_threshold = read_f64(r)?;
    let backend = backend_from_u8(read_u8(r)?)?;
    let mpw_present = read_u8(r)?;
    let mpw_value = read_u64(r)?;
    let max_pauli_weight = match mpw_present {
        0 => None,
        1 => Some(usize::try_from(mpw_value).map_err(|_| {
            eyre::eyre!("max_pauli_weight {mpw_value} does not fit in usize on this platform")
        })?),
        other => {
            return Err(eyre::eyre!(
                "invalid max_pauli_weight presence byte {other}"
            ));
        }
    };
    let observable_present = read_u8(r)?;
    let observable = match observable_present {
        0 => None,
        1 => {
            let len = read_u32(r)? as usize;
            let mut bytes = vec![0u8; len];
            r.read_exact(&mut bytes)?;
            Some(String::from_utf8(bytes)?)
        }
        other => {
            return Err(eyre::eyre!("invalid observable presence byte {other}"));
        }
    };

    // Sections begin at `header_size`; skip any header bytes beyond what this
    // reader knows about (forward compat / self-description).
    let consumed = FIXED_HEADER_SIZE
        + if observable.is_some() {
            4 + u32::try_from(observable.as_deref().unwrap().len())
                .map_err(|_| eyre::eyre!("observable length does not fit in u32"))?
        } else {
            0
        };
    if header_size < consumed {
        return Err(eyre::eyre!(
            "header_size {header_size} smaller than the {consumed} bytes already consumed"
        ));
    }
    skip_bytes(r, u64::from(header_size - consumed))?;

    // Don't pre-allocate from an untrusted count; grow as entries are read.
    let string_count = read_u32(r)?;
    let mut strings = Vec::new();
    for _ in 0..string_count {
        let len = read_u32(r)? as usize;
        let mut bytes = vec![0u8; len];
        r.read_exact(&mut bytes)?;
        strings.push(String::from_utf8(bytes)?);
    }

    let function_count = read_u32(r)?;
    let mut functions = Vec::new();
    for _ in 0..function_count {
        let name = read_u32(r)?;
        let local_count = read_u32(r)?;
        let start_address = read_u32(r)?;
        let end_address = read_u32(r)?;
        let file = read_u32(r)?;
        let parameter_count = read_u32(r)?;
        let mut params = Vec::new();
        for _ in 0..parameter_count {
            params.push(Parameter {
                name: read_u32(r)?,
                ty: Type::from_bytes(r)?,
            });
        }
        let return_count = read_u32(r)?;
        let mut ret = Vec::new();
        for _ in 0..return_count {
            ret.push(Type::from_bytes(r)?);
        }
        functions.push(FunctionInfo {
            name,
            signature: Signature { params, ret },
            local_count,
            start_address,
            end_address,
            file,
        });
    }
    let label_count = read_u32(r)?;
    let mut labels = Vec::new();
    for _ in 0..label_count {
        labels.push(LabelInfo {
            address: read_u32(r)?,
            name: read_u32(r)?,
        });
    }
    let main = read_u32(r)?;
    let main_function = (main != u32::MAX).then_some(main);

    let code_count = read_u32(r)?;
    let mut code = Vec::new();
    for _ in 0..code_count {
        code.push(read_instruction(r)?);
    }

    Ok(PPVMModule {
        extra: PPVMDeviceInfo {
            magic,
            n_qubits,
            coefficient_threshold,
            backend,
            observable,
            max_pauli_weight,
        },
        strings,
        functions,
        labels,
        main_function,
        code,
        ..Default::default()
    })
}

fn backend_to_u8(backend: BackendKind) -> u8 {
    match backend {
        BackendKind::Tableau => 0,
        BackendKind::PauliSum => 1,
        BackendKind::LossyPauliSum => 2,
    }
}

fn backend_from_u8(byte: u8) -> eyre::Result<BackendKind> {
    match byte {
        0 => Ok(BackendKind::Tableau),
        1 => Ok(BackendKind::PauliSum),
        2 => Ok(BackendKind::LossyPauliSum),
        other => Err(eyre::eyre!("invalid backend tag {other}")),
    }
}

/// Serialize a module to an owned byte vector.
pub fn module_to_bytes(module: &PPVMModule) -> eyre::Result<Vec<u8>> {
    let mut buf = Vec::new();
    write_module(module, &mut buf)?;
    Ok(buf)
}

/// Reconstruct a module from a byte slice.
pub fn module_from_bytes(bytes: &[u8]) -> eyre::Result<PPVMModule> {
    read_module(&mut &bytes[..])
}

/// Cheap sniff: does this byte stream begin with the PPVM `.ssb` magic?
///
/// Reads the leading four bytes the same way [`read_module`] does — as a
/// little-endian `u32` — so a positive result here means [`read_module`] will
/// accept the magic. A stream shorter than the magic is not bytecode.
pub fn is_bytecode(bytes: &[u8]) -> bool {
    bytes.len() >= 4 && u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) == PPVM_MAGIC
}

/// "Dump": compile `.sst` source straight to the `.ssb` byte stream.
pub fn compile_to_bytes(source: &str) -> eyre::Result<Vec<u8>> {
    let module = crate::compile_program(source)?;
    module_to_bytes(&module)
}

fn read_u8<R: Read>(r: &mut R) -> eyre::Result<u8> {
    let mut b = [0u8; 1];
    r.read_exact(&mut b)?;
    Ok(b[0])
}

fn write_u32<W: Write>(w: &mut W, value: u32) -> eyre::Result<()> {
    w.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn read_u16<R: Read>(r: &mut R) -> eyre::Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}

fn read_u32<R: Read>(r: &mut R) -> eyre::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: Read>(r: &mut R) -> eyre::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn read_f64<R: Read>(r: &mut R) -> eyre::Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

fn skip_bytes<R: Read>(r: &mut R, n: u64) -> eyre::Result<()> {
    let skipped = std::io::copy(&mut r.take(n), &mut std::io::sink())?;
    if skipped != n {
        return Err(eyre::eyre!("unexpected EOF skipping {n} header bytes"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compile_program;

    fn round_trip(module: &PPVMModule) -> PPVMModule {
        module_from_bytes(&module_to_bytes(module).unwrap()).unwrap()
    }

    fn fixture_paths() -> Vec<std::path::PathBuf> {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let mut paths: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|ext| ext == "sst"))
            .collect();
        paths.sort();
        paths
    }

    #[test]
    fn every_instruction_round_trips() {
        let mut m = PPVMModule::default();
        m.extra.n_qubits = 1;
        m.code = every_instruction();
        assert_eq!(
            m.code.len(),
            Opcode::V_debugger__debug__Breakpoint as usize + 1
        );
        assert_eq!(round_trip(&m), m);
    }

    #[test]
    fn every_fixture_round_trips() {
        let paths = fixture_paths();
        assert!(!paths.is_empty());
        for path in paths {
            let module = compile_program(&std::fs::read_to_string(&path).unwrap())
                .unwrap_or_else(|e| panic!("compile {}: {e:?}", path.display()));
            assert_eq!(round_trip(&module), module, "{}", path.display());
        }
    }

    #[test]
    fn round_trips_device_info() {
        let mut m = PPVMModule::default();
        m.extra.n_qubits = 6;
        m.extra.coefficient_threshold = 1e-9;
        m.extra.backend = BackendKind::PauliSum;
        m.extra.observable = Some("ZZIIII".to_owned());
        m.extra.max_pauli_weight = Some(8);
        assert_eq!(round_trip(&m), m);

        // Absent optional fields take the other encoding path.
        m.extra.backend = BackendKind::LossyPauliSum;
        m.extra.observable = None;
        m.extra.max_pauli_weight = None;
        assert_eq!(round_trip(&m), m);
    }

    #[test]
    fn round_trips_strings() {
        let m = PPVMModule {
            strings: vec![
                String::new(),
                "hello".to_owned(),
                "tab\tnl\nquote\"".to_owned(),
                "üñîçødé ⚛".to_owned(),
            ],
            ..Default::default()
        };
        assert_eq!(round_trip(&m), m);
    }

    #[test]
    fn bytecode_runs_like_source() {
        // Same program and seed, from source and from bytecode: identical records.
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/bell.sst");
        let module = compile_program(&std::fs::read_to_string(path).unwrap()).unwrap();
        let loaded = round_trip(&module);
        let records = |m: &PPVMModule| {
            crate::shots::run_shots(m, 64, Some(3))
                .unwrap()
                .into_iter()
                .map(|shot| shot.measurements)
                .collect::<Vec<_>>()
        };
        assert_eq!(records(&loaded), records(&module));
    }

    #[test]
    fn read_honors_header_size_with_padding() {
        let mut m = PPVMModule::default();
        m.extra.n_qubits = 3;
        m.strings = vec!["hi".to_owned()];
        let mut buf = module_to_bytes(&m).unwrap();
        // Grow the header by 4 padding bytes; the reader must skip past them.
        buf[6..10].copy_from_slice(&(FIXED_HEADER_SIZE + 4).to_le_bytes());
        for i in 0..4 {
            buf.insert(FIXED_HEADER_SIZE as usize + i, 0);
        }
        assert_eq!(module_from_bytes(&buf).unwrap(), m);
    }

    #[test]
    fn read_rejects_bad_magic() {
        let mut buf = module_to_bytes(&PPVMModule::default()).unwrap();
        buf[0] ^= 0xFF;
        assert!(!is_bytecode(&buf));
        let err = module_from_bytes(&buf).unwrap_err();
        assert!(
            err.to_string().contains("not a PPVM bytecode file"),
            "{err}"
        );
    }

    #[test]
    fn read_rejects_unsupported_version() {
        let mut buf = module_to_bytes(&PPVMModule::default()).unwrap();
        buf[4..6].copy_from_slice(&(PPVM_BYTECODE_VERSION + 1).to_le_bytes());
        let err = module_from_bytes(&buf).unwrap_err();
        assert!(
            err.to_string().contains("unsupported bytecode version"),
            "{err}"
        );
    }

    #[test]
    fn read_rejects_unknown_opcode() {
        let m = PPVMModule {
            code: vec![every_instruction().remove(0)],
            ..Default::default()
        };
        let mut buf = module_to_bytes(&m).unwrap();
        // The opcode is the first two bytes of the last instruction (`br` + u32).
        let at = buf.len() - 6;
        buf[at..at + 2].copy_from_slice(&u16::MAX.to_le_bytes());
        let err = module_from_bytes(&buf).unwrap_err();
        assert!(err.to_string().contains("unknown opcode"), "{err}");
    }

    #[test]
    fn read_rejects_truncated_input() {
        let m = PPVMModule {
            code: vec![every_instruction().remove(0)],
            ..Default::default()
        };
        let mut buf = module_to_bytes(&m).unwrap();
        buf.truncate(buf.len() - 3); // cut off mid-instruction
        assert!(module_from_bytes(&buf).is_err());
    }
}
