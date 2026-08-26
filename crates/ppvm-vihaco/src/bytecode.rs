// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

//! Binary bytecode (`.ssb`) round-trip for PPVM modules.
//!
//! Serializes a resolved [`Module`] to a little-endian container — header
//! (magic + version + device info) → strings section → code section — and
//! reads it back. The per-instruction codec is the existing `WriteBytes` /
//! `FromBytes` generated for [`PPVMInstruction`]; this module only frames the
//! container around it.

use std::io::{Read, Write};

use vihaco::instruction::{FromBytes, WriteBytes};
use vihaco::module::{FunctionInfo, LabelInfo, Parameter, Signature};
use vihaco::{Type, Value};
use vihaco_cpu::RuntimeInstruction as CpuInstruction;

use crate::PPVMModule;
use crate::composite::{BackendKind, PPVM_MAGIC, PPVMDeviceInfo, PPVMInstruction};

#[derive(Debug, Clone, vihaco::Instruction)]
enum BytecodeType {
    Undefined,
    String,
    Bool,
    I64,
    U32,
    U64,
    F64,
    FunctionRef,
    HeapRef,
}

#[derive(Debug, Clone, vihaco::Instruction)]
enum BytecodeValue {
    Undefined,
    String(u32),
    Bool(bool),
    I64(i64),
    U32(u32),
    U64(u64),
    F64(f64),
    FunctionRef(u32),
    HeapRef(u32),
}

fn encode_type(value: Type) -> BytecodeType {
    match value {
        Type::Undefined => BytecodeType::Undefined,
        Type::String => BytecodeType::String,
        Type::Bool => BytecodeType::Bool,
        Type::I64 => BytecodeType::I64,
        Type::U32 => BytecodeType::U32,
        Type::U64 => BytecodeType::U64,
        Type::F64 => BytecodeType::F64,
        Type::FunctionRef => BytecodeType::FunctionRef,
        Type::HeapRef => BytecodeType::HeapRef,
    }
}

fn decode_type(value: BytecodeType) -> Type {
    match value {
        BytecodeType::Undefined => Type::Undefined,
        BytecodeType::String => Type::String,
        BytecodeType::Bool => Type::Bool,
        BytecodeType::I64 => Type::I64,
        BytecodeType::U32 => Type::U32,
        BytecodeType::U64 => Type::U64,
        BytecodeType::F64 => Type::F64,
        BytecodeType::FunctionRef => Type::FunctionRef,
        BytecodeType::HeapRef => Type::HeapRef,
    }
}

fn encode_value(value: Value) -> BytecodeValue {
    match value {
        Value::Undefined => BytecodeValue::Undefined,
        Value::String(v) => BytecodeValue::String(v),
        Value::Bool(v) => BytecodeValue::Bool(v),
        Value::I64(v) => BytecodeValue::I64(v),
        Value::U32(v) => BytecodeValue::U32(v),
        Value::U64(v) => BytecodeValue::U64(v),
        Value::F64(v) => BytecodeValue::F64(v),
        Value::FunctionRef(v) => BytecodeValue::FunctionRef(v),
        Value::HeapRef(v) => BytecodeValue::HeapRef(v),
    }
}

fn decode_value(value: BytecodeValue) -> Value {
    match value {
        BytecodeValue::Undefined => Value::Undefined,
        BytecodeValue::String(v) => Value::String(v),
        BytecodeValue::Bool(v) => Value::Bool(v),
        BytecodeValue::I64(v) => Value::I64(v),
        BytecodeValue::U32(v) => Value::U32(v),
        BytecodeValue::U64(v) => Value::U64(v),
        BytecodeValue::F64(v) => Value::F64(v),
        BytecodeValue::FunctionRef(v) => Value::FunctionRef(v),
        BytecodeValue::HeapRef(v) => Value::HeapRef(v),
    }
}

#[derive(Debug, Clone, vihaco::Instruction)]
enum BytecodeCpu {
    Span(u32, u32, u32),
    FunctionStart,
    FunctionEnd,
    Breakpoint,
    Branch(u32),
    ConditionalBranch(u32, u32),
    Return(u32),
    IndirectCall,
    Call(u32, u32),
    Halt,
    Print,
    Load(BytecodeType, u32),
    Store(BytecodeType, u32),
    Dup,
    HeapAlloc(u32),
    GetItem,
    HeapDealloc,
    Const(BytecodeType, BytecodeValue),
    Add(BytecodeType),
    Sub(BytecodeType),
    Mul(BytecodeType),
    Div(BytecodeType),
    Rem(BytecodeType),
    Neg(BytecodeType),
    Shl(BytecodeType),
    Shr(BytecodeType),
    Rol(BytecodeType),
    Ror(BytecodeType),
    BitAnd(BytecodeType),
    BitOr(BytecodeType),
    BitXor(BytecodeType),
    Not,
    And,
    Or,
    Xor,
    Eq(BytecodeType),
    Ne(BytecodeType),
    Lt(BytecodeType),
    Gt(BytecodeType),
    Le(BytecodeType),
    Ge(BytecodeType),
}

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug, Clone, vihaco::Instruction)]
enum BytecodeCircuit {
    TwoQubitPauliError,
    Truncate,
    Trace,
    X,
    Y,
    Z,
    H,
    SqrtXAdj,
    SqrtX,
    SqrtYAdj,
    SqrtY,
    SAdj,
    S,
    CNOT,
    CZ,
    TAdj,
    T,
    RXX,
    RYY,
    RZZ,
    RX,
    RY,
    RZ,
    U3,
    Measure,
    Reset,
    R,
    Loss,
    CorrelatedLoss,
    PauliError,
    Depolarize2,
    Depolarize,
}

#[derive(Debug, Clone, vihaco::Instruction)]
enum BytecodeInstruction {
    Cpu(BytecodeCpu),
    Circuit(BytecodeCircuit),
}

fn encode_instruction(inst: &PPVMInstruction) -> eyre::Result<BytecodeInstruction> {
    let encoded = match inst {
        PPVMInstruction::Cpu(inst) => BytecodeInstruction::Cpu(match inst {
            CpuInstruction::Span(a, b, c) => BytecodeCpu::Span(*a, *b, *c),
            CpuInstruction::FunctionStart => BytecodeCpu::FunctionStart,
            CpuInstruction::FunctionEnd => BytecodeCpu::FunctionEnd,
            CpuInstruction::Breakpoint => BytecodeCpu::Breakpoint,
            CpuInstruction::Branch(v) => BytecodeCpu::Branch(*v),
            CpuInstruction::ConditionalBranch(a, b) => BytecodeCpu::ConditionalBranch(*a, *b),
            CpuInstruction::Return(v) => BytecodeCpu::Return(*v),
            CpuInstruction::IndirectCall => BytecodeCpu::IndirectCall,
            CpuInstruction::Call(a, b) => BytecodeCpu::Call(*a, *b),
            CpuInstruction::Halt => BytecodeCpu::Halt,
            CpuInstruction::Print => BytecodeCpu::Print,
            CpuInstruction::Load(t, v) => BytecodeCpu::Load(encode_type(*t), *v),
            CpuInstruction::Store(t, v) => BytecodeCpu::Store(encode_type(*t), *v),
            CpuInstruction::Dup => BytecodeCpu::Dup,
            CpuInstruction::HeapAlloc(v) => BytecodeCpu::HeapAlloc(*v),
            CpuInstruction::GetItem => BytecodeCpu::GetItem,
            CpuInstruction::HeapDealloc => BytecodeCpu::HeapDealloc,
            CpuInstruction::Const(t, v) => BytecodeCpu::Const(encode_type(*t), encode_value(*v)),
            CpuInstruction::Add(t) => BytecodeCpu::Add(encode_type(*t)),
            CpuInstruction::Sub(t) => BytecodeCpu::Sub(encode_type(*t)),
            CpuInstruction::Mul(t) => BytecodeCpu::Mul(encode_type(*t)),
            CpuInstruction::Div(t) => BytecodeCpu::Div(encode_type(*t)),
            CpuInstruction::Rem(t) => BytecodeCpu::Rem(encode_type(*t)),
            CpuInstruction::Neg(t) => BytecodeCpu::Neg(encode_type(*t)),
            CpuInstruction::Shl(t) => BytecodeCpu::Shl(encode_type(*t)),
            CpuInstruction::Shr(t) => BytecodeCpu::Shr(encode_type(*t)),
            CpuInstruction::Rol(t) => BytecodeCpu::Rol(encode_type(*t)),
            CpuInstruction::Ror(t) => BytecodeCpu::Ror(encode_type(*t)),
            CpuInstruction::BitAnd(t) => BytecodeCpu::BitAnd(encode_type(*t)),
            CpuInstruction::BitOr(t) => BytecodeCpu::BitOr(encode_type(*t)),
            CpuInstruction::BitXor(t) => BytecodeCpu::BitXor(encode_type(*t)),
            CpuInstruction::Not => BytecodeCpu::Not,
            CpuInstruction::And => BytecodeCpu::And,
            CpuInstruction::Or => BytecodeCpu::Or,
            CpuInstruction::Xor => BytecodeCpu::Xor,
            CpuInstruction::Eq(t) => BytecodeCpu::Eq(encode_type(*t)),
            CpuInstruction::Ne(t) => BytecodeCpu::Ne(encode_type(*t)),
            CpuInstruction::Lt(t) => BytecodeCpu::Lt(encode_type(*t)),
            CpuInstruction::Gt(t) => BytecodeCpu::Gt(encode_type(*t)),
            CpuInstruction::Le(t) => BytecodeCpu::Le(encode_type(*t)),
            CpuInstruction::Ge(t) => BytecodeCpu::Ge(encode_type(*t)),
            CpuInstruction::Label(_) => {
                return Err(eyre::eyre!(
                    "runtime labels are not serializable in PPVM bytecode"
                ));
            }
        }),
        PPVMInstruction::Circuit(inst) => BytecodeInstruction::Circuit(match inst {
            vihaco_circuit_isa::CircuitInstruction::TwoQubitPauliError => {
                BytecodeCircuit::TwoQubitPauliError
            }
            vihaco_circuit_isa::CircuitInstruction::Truncate => BytecodeCircuit::Truncate,
            vihaco_circuit_isa::CircuitInstruction::Trace => BytecodeCircuit::Trace,
            vihaco_circuit_isa::CircuitInstruction::X => BytecodeCircuit::X,
            vihaco_circuit_isa::CircuitInstruction::Y => BytecodeCircuit::Y,
            vihaco_circuit_isa::CircuitInstruction::Z => BytecodeCircuit::Z,
            vihaco_circuit_isa::CircuitInstruction::H => BytecodeCircuit::H,
            vihaco_circuit_isa::CircuitInstruction::SqrtXAdj => BytecodeCircuit::SqrtXAdj,
            vihaco_circuit_isa::CircuitInstruction::SqrtX => BytecodeCircuit::SqrtX,
            vihaco_circuit_isa::CircuitInstruction::SqrtYAdj => BytecodeCircuit::SqrtYAdj,
            vihaco_circuit_isa::CircuitInstruction::SqrtY => BytecodeCircuit::SqrtY,
            vihaco_circuit_isa::CircuitInstruction::SAdj => BytecodeCircuit::SAdj,
            vihaco_circuit_isa::CircuitInstruction::S => BytecodeCircuit::S,
            vihaco_circuit_isa::CircuitInstruction::CNOT => BytecodeCircuit::CNOT,
            vihaco_circuit_isa::CircuitInstruction::CZ => BytecodeCircuit::CZ,
            vihaco_circuit_isa::CircuitInstruction::TAdj => BytecodeCircuit::TAdj,
            vihaco_circuit_isa::CircuitInstruction::T => BytecodeCircuit::T,
            vihaco_circuit_isa::CircuitInstruction::RXX => BytecodeCircuit::RXX,
            vihaco_circuit_isa::CircuitInstruction::RYY => BytecodeCircuit::RYY,
            vihaco_circuit_isa::CircuitInstruction::RZZ => BytecodeCircuit::RZZ,
            vihaco_circuit_isa::CircuitInstruction::RX => BytecodeCircuit::RX,
            vihaco_circuit_isa::CircuitInstruction::RY => BytecodeCircuit::RY,
            vihaco_circuit_isa::CircuitInstruction::RZ => BytecodeCircuit::RZ,
            vihaco_circuit_isa::CircuitInstruction::U3 => BytecodeCircuit::U3,
            vihaco_circuit_isa::CircuitInstruction::Measure => BytecodeCircuit::Measure,
            vihaco_circuit_isa::CircuitInstruction::Reset => BytecodeCircuit::Reset,
            vihaco_circuit_isa::CircuitInstruction::R => BytecodeCircuit::R,
            vihaco_circuit_isa::CircuitInstruction::Loss => BytecodeCircuit::Loss,
            vihaco_circuit_isa::CircuitInstruction::CorrelatedLoss => {
                BytecodeCircuit::CorrelatedLoss
            }
            vihaco_circuit_isa::CircuitInstruction::PauliError => BytecodeCircuit::PauliError,
            vihaco_circuit_isa::CircuitInstruction::Depolarize2 => BytecodeCircuit::Depolarize2,
            vihaco_circuit_isa::CircuitInstruction::Depolarize => BytecodeCircuit::Depolarize,
        }),
    };
    Ok(encoded)
}

fn decode_instruction(inst: BytecodeInstruction) -> PPVMInstruction {
    match inst {
        BytecodeInstruction::Cpu(inst) => PPVMInstruction::Cpu(match inst {
            BytecodeCpu::Span(a, b, c) => CpuInstruction::Span(a, b, c),
            BytecodeCpu::FunctionStart => CpuInstruction::FunctionStart,
            BytecodeCpu::FunctionEnd => CpuInstruction::FunctionEnd,
            BytecodeCpu::Breakpoint => CpuInstruction::Breakpoint,
            BytecodeCpu::Branch(v) => CpuInstruction::Branch(v),
            BytecodeCpu::ConditionalBranch(a, b) => CpuInstruction::ConditionalBranch(a, b),
            BytecodeCpu::Return(v) => CpuInstruction::Return(v),
            BytecodeCpu::IndirectCall => CpuInstruction::IndirectCall,
            BytecodeCpu::Call(a, b) => CpuInstruction::Call(a, b),
            BytecodeCpu::Halt => CpuInstruction::Halt,
            BytecodeCpu::Print => CpuInstruction::Print,
            BytecodeCpu::Load(t, v) => CpuInstruction::Load(decode_type(t), v),
            BytecodeCpu::Store(t, v) => CpuInstruction::Store(decode_type(t), v),
            BytecodeCpu::Dup => CpuInstruction::Dup,
            BytecodeCpu::HeapAlloc(v) => CpuInstruction::HeapAlloc(v),
            BytecodeCpu::GetItem => CpuInstruction::GetItem,
            BytecodeCpu::HeapDealloc => CpuInstruction::HeapDealloc,
            BytecodeCpu::Const(t, v) => CpuInstruction::Const(decode_type(t), decode_value(v)),
            BytecodeCpu::Add(t) => CpuInstruction::Add(decode_type(t)),
            BytecodeCpu::Sub(t) => CpuInstruction::Sub(decode_type(t)),
            BytecodeCpu::Mul(t) => CpuInstruction::Mul(decode_type(t)),
            BytecodeCpu::Div(t) => CpuInstruction::Div(decode_type(t)),
            BytecodeCpu::Rem(t) => CpuInstruction::Rem(decode_type(t)),
            BytecodeCpu::Neg(t) => CpuInstruction::Neg(decode_type(t)),
            BytecodeCpu::Shl(t) => CpuInstruction::Shl(decode_type(t)),
            BytecodeCpu::Shr(t) => CpuInstruction::Shr(decode_type(t)),
            BytecodeCpu::Rol(t) => CpuInstruction::Rol(decode_type(t)),
            BytecodeCpu::Ror(t) => CpuInstruction::Ror(decode_type(t)),
            BytecodeCpu::BitAnd(t) => CpuInstruction::BitAnd(decode_type(t)),
            BytecodeCpu::BitOr(t) => CpuInstruction::BitOr(decode_type(t)),
            BytecodeCpu::BitXor(t) => CpuInstruction::BitXor(decode_type(t)),
            BytecodeCpu::Not => CpuInstruction::Not,
            BytecodeCpu::And => CpuInstruction::And,
            BytecodeCpu::Or => CpuInstruction::Or,
            BytecodeCpu::Xor => CpuInstruction::Xor,
            BytecodeCpu::Eq(t) => CpuInstruction::Eq(decode_type(t)),
            BytecodeCpu::Ne(t) => CpuInstruction::Ne(decode_type(t)),
            BytecodeCpu::Lt(t) => CpuInstruction::Lt(decode_type(t)),
            BytecodeCpu::Gt(t) => CpuInstruction::Gt(decode_type(t)),
            BytecodeCpu::Le(t) => CpuInstruction::Le(decode_type(t)),
            BytecodeCpu::Ge(t) => CpuInstruction::Ge(decode_type(t)),
        }),
        BytecodeInstruction::Circuit(inst) => PPVMInstruction::Circuit(match inst {
            BytecodeCircuit::TwoQubitPauliError => {
                vihaco_circuit_isa::CircuitInstruction::TwoQubitPauliError
            }
            BytecodeCircuit::Truncate => vihaco_circuit_isa::CircuitInstruction::Truncate,
            BytecodeCircuit::Trace => vihaco_circuit_isa::CircuitInstruction::Trace,
            BytecodeCircuit::X => vihaco_circuit_isa::CircuitInstruction::X,
            BytecodeCircuit::Y => vihaco_circuit_isa::CircuitInstruction::Y,
            BytecodeCircuit::Z => vihaco_circuit_isa::CircuitInstruction::Z,
            BytecodeCircuit::H => vihaco_circuit_isa::CircuitInstruction::H,
            BytecodeCircuit::SqrtXAdj => vihaco_circuit_isa::CircuitInstruction::SqrtXAdj,
            BytecodeCircuit::SqrtX => vihaco_circuit_isa::CircuitInstruction::SqrtX,
            BytecodeCircuit::SqrtYAdj => vihaco_circuit_isa::CircuitInstruction::SqrtYAdj,
            BytecodeCircuit::SqrtY => vihaco_circuit_isa::CircuitInstruction::SqrtY,
            BytecodeCircuit::SAdj => vihaco_circuit_isa::CircuitInstruction::SAdj,
            BytecodeCircuit::S => vihaco_circuit_isa::CircuitInstruction::S,
            BytecodeCircuit::CNOT => vihaco_circuit_isa::CircuitInstruction::CNOT,
            BytecodeCircuit::CZ => vihaco_circuit_isa::CircuitInstruction::CZ,
            BytecodeCircuit::TAdj => vihaco_circuit_isa::CircuitInstruction::TAdj,
            BytecodeCircuit::T => vihaco_circuit_isa::CircuitInstruction::T,
            BytecodeCircuit::RXX => vihaco_circuit_isa::CircuitInstruction::RXX,
            BytecodeCircuit::RYY => vihaco_circuit_isa::CircuitInstruction::RYY,
            BytecodeCircuit::RZZ => vihaco_circuit_isa::CircuitInstruction::RZZ,
            BytecodeCircuit::RX => vihaco_circuit_isa::CircuitInstruction::RX,
            BytecodeCircuit::RY => vihaco_circuit_isa::CircuitInstruction::RY,
            BytecodeCircuit::RZ => vihaco_circuit_isa::CircuitInstruction::RZ,
            BytecodeCircuit::U3 => vihaco_circuit_isa::CircuitInstruction::U3,
            BytecodeCircuit::Measure => vihaco_circuit_isa::CircuitInstruction::Measure,
            BytecodeCircuit::Reset => vihaco_circuit_isa::CircuitInstruction::Reset,
            BytecodeCircuit::R => vihaco_circuit_isa::CircuitInstruction::R,
            BytecodeCircuit::Loss => vihaco_circuit_isa::CircuitInstruction::Loss,
            BytecodeCircuit::CorrelatedLoss => {
                vihaco_circuit_isa::CircuitInstruction::CorrelatedLoss
            }
            BytecodeCircuit::PauliError => vihaco_circuit_isa::CircuitInstruction::PauliError,
            BytecodeCircuit::Depolarize2 => vihaco_circuit_isa::CircuitInstruction::Depolarize2,
            BytecodeCircuit::Depolarize => vihaco_circuit_isa::CircuitInstruction::Depolarize,
        }),
    }
}

/// Current `.ssb` format version. The reader rejects any other version.
pub const PPVM_BYTECODE_VERSION: u16 = 2;

/// Byte length of the fixed portion of the header. The actual `header_size`
/// in the stream may exceed this when the optional `observable` string is
/// populated; the reader uses `header_size` to skip to the strings section.
///
/// Field widths (bytes): magic(4) + version(2) + header_size(4) + n_qubits(4)
/// + coefficient_threshold(8) + backend(1) + max_pauli_weight_present(1)
/// + max_pauli_weight(8) + observable_present(1) = 33.
const FIXED_HEADER_SIZE: u32 = 4 + 2 + 4 + 4 + 8 + 1 + 1 + 8 + 1;

/// Serialize a resolved module to the v2 `.ssb` byte stream.
pub fn write_module<W: Write>(module: &PPVMModule, w: &mut W) -> eyre::Result<()> {
    if !module.constants.is_empty() || !module.source_symbols.is_empty() {
        return Err(eyre::eyre!(
            "bytecode v2 cannot represent constants or source symbols"
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
        encode_instruction(inst)?.write_bytes(w)?;
    }

    Ok(())
}

/// Reconstruct a module from a v2 `.ssb` byte stream.
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
        code.push(decode_instruction(BytecodeInstruction::from_bytes(r)?));
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
    use vihaco::{Type, Value};

    use super::*;

    fn empty_module() -> PPVMModule {
        PPVMModule::default()
    }

    #[test]
    fn round_trips_device_info() {
        let mut m = empty_module();
        m.extra.n_qubits = 7;
        m.extra.coefficient_threshold = 1e-9;

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn round_trips_paulisum_device_info() {
        let mut m = empty_module();
        m.extra.n_qubits = 6;
        m.extra.backend = BackendKind::PauliSum;
        m.extra.observable = Some("ZZIIII".to_string());
        m.extra.max_pauli_weight = Some(8);

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn round_trips_lossy_backend_without_observable() {
        let mut m = empty_module();
        m.extra.n_qubits = 4;
        m.extra.backend = BackendKind::LossyPauliSum;
        // observable and max_pauli_weight stay None — verifies the absent path.

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
        assert_eq!(back.extra.observable, None);
        assert_eq!(back.extra.max_pauli_weight, None);
    }

    #[test]
    fn round_trips_code() {
        use vihaco_circuit_isa::CircuitInstruction;
        use vihaco_cpu::RuntimeInstruction as Cpu;

        let mut m = empty_module();
        m.extra.n_qubits = 2;
        m.code = vec![
            PPVMInstruction::Cpu(Cpu::Const(Type::U64, Value::U64(0))),
            PPVMInstruction::Circuit(CircuitInstruction::H),
            PPVMInstruction::Circuit(CircuitInstruction::R),
            PPVMInstruction::Cpu(Cpu::Branch(1)),
            PPVMInstruction::Cpu(Cpu::ConditionalBranch(0, 1)),
            PPVMInstruction::Cpu(Cpu::Call(0, 1)),
            PPVMInstruction::Cpu(Cpu::Return(0)),
        ];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn read_honors_header_size_with_padding() {
        let mut m = empty_module();
        m.extra.n_qubits = 3;
        m.strings = vec!["hi".to_string()];
        m.code = vec![PPVMInstruction::Cpu(
            vihaco_cpu::RuntimeInstruction::Return(0),
        )];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();

        // Simulate a larger header: 4 padding bytes after the fixed fields,
        // with header_size bumped to match. The reader must skip to it.
        // (This test uses an empty observable, so the on-disk header size
        // equals FIXED_HEADER_SIZE.)
        buf[6..10].copy_from_slice(&(FIXED_HEADER_SIZE + 4).to_le_bytes());
        for i in 0..4 {
            buf.insert(FIXED_HEADER_SIZE as usize + i, 0x00);
        }

        let back = read_module(&mut buf.as_slice()).unwrap();
        assert_eq!(back, m);
    }

    #[test]
    fn compile_to_bytes_round_trips_through_resolve() {
        let src = "device circuit.n_qubits 2;\n\
                   fn @main() {\n\
                       cpu::cpu.const u64, 0\n\
                       circuit::circuit.h\n\
                       cpu::cpu.const u64, 0\n\
                       cpu::cpu.const u64, 1\n\
                       circuit::circuit.cnot\n\
                       cpu::cpu.ret 0\n\
                   }\n";

        let bytes = compile_to_bytes(src).unwrap();
        let back = module_from_bytes(&bytes).unwrap();
        let expected = crate::compile_program(src).unwrap();

        assert_eq!(back, expected);
    }

    #[test]
    fn loaded_bytecode_executes_like_text() {
        let src = "device circuit.n_qubits 2;\n\
                   fn @main() {\n\
                       cpu::cpu.const u64, 0\n circuit::circuit.h\n\
                       cpu::cpu.const u64, 0\n cpu::cpu.const u64, 1\n circuit::circuit.cnot\n\
                       cpu::cpu.const u64, 0\n circuit::circuit.measure\n\
                       cpu::cpu.const u64, 1\n circuit::circuit.measure\n\
                        cpu::cpu.ret 0\n }\n";
        let bytes = compile_to_bytes(src).unwrap();

        let mut machine = crate::composite::PPVM::default();
        machine.load_bytecode(&bytes).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 2);
    }

    #[test]
    fn load_bytecode_file_reads_from_disk() {
        let src = "device circuit.n_qubits 1;\n\
                   fn @main() { cpu::cpu.const u64, 0\n circuit::circuit.measure\n cpu::cpu.ret 0 }\n";
        let bytes = compile_to_bytes(src).unwrap();
        let path = std::env::temp_dir().join("ppvm_load_bytecode_file_test.ssb");
        std::fs::write(&path, &bytes).unwrap();

        let mut machine = crate::composite::PPVM::default();
        machine.load_bytecode_file(path.to_str().unwrap()).unwrap();
        machine.run().unwrap();

        assert_eq!(machine.measurement_record().len(), 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn read_rejects_truncated_input() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        m.code = vec![PPVMInstruction::Cpu(
            vihaco_cpu::RuntimeInstruction::Return(0),
        )];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        buf.truncate(buf.len() - 3); // cut off mid-instruction

        assert!(read_module(&mut buf.as_slice()).is_err());
    }

    #[test]
    fn round_trips_populated_functions_table() {
        use vihaco::module::{FunctionInfo, Signature};

        let mut m = empty_module();
        m.extra.n_qubits = 1;
        m.functions.push(FunctionInfo {
            name: 0,
            signature: Signature {
                params: vec![],
                ret: vec![],
            },
            local_count: 0,
            start_address: 0,
            end_address: 0,
            file: 0,
        });

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        assert_eq!(read_module(&mut buf.as_slice()).unwrap(), m);
    }

    #[test]
    fn read_rejects_bad_magic() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        // Corrupt the magic (first 4 bytes).
        buf[0] ^= 0xFF;

        let err = read_module(&mut buf.as_slice()).unwrap_err();
        assert!(
            err.to_string().contains("not a PPVM bytecode file"),
            "err: {err}"
        );
    }

    #[test]
    fn round_trips_strings() {
        let mut m = empty_module();
        m.extra.n_qubits = 1;
        m.strings = vec![
            String::new(),
            "hello".to_string(),
            "tab\tnl\nquote\"".to_string(),
            "üñîçødé ⚛".to_string(),
        ];

        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        let back = read_module(&mut buf.as_slice()).unwrap();

        assert_eq!(back, m);
    }

    #[test]
    fn read_rejects_unsupported_version() {
        let mut m = empty_module();
        m.extra.n_qubits = 2;
        let mut buf = Vec::new();
        write_module(&m, &mut buf).unwrap();
        // Bump the version (bytes 4..6) past the supported one.
        let bad = (PPVM_BYTECODE_VERSION + 1).to_le_bytes();
        buf[4] = bad[0];
        buf[5] = bad[1];

        let err = read_module(&mut buf.as_slice()).unwrap_err();
        assert!(
            err.to_string().contains("unsupported bytecode version"),
            "err: {err}"
        );
    }
}
